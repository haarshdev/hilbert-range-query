import json
import math
import pandas as pd


# =========================
# INPUTS (edit if needed)
# =========================
CLUSTER_STATUS_JSON = "cluster-status-15012026.json"
TABLEB_5D = "tableB_nodes_5D.csv"

# numerical tolerance in ms
EPS_MS = 1e-9


# =========================
# HELPERS
# =========================
def short_name(full: str) -> str:
    """
    Converts 'clab-nebula-serf4' -> 'serf4'.
    If already 'serf4', returns it unchanged.
    """
    s = str(full)
    if s.startswith("serf"):
        return s
    if "serf" in s and "-" in s:
        return s.split("-")[-1]
    return s


def euclid(vec_a, vec_b) -> float:
    """Euclidean distance in the same units as Vec (seconds in Serf)."""
    if len(vec_a) != len(vec_b):
        raise ValueError(f"Vec dim mismatch: {len(vec_a)} vs {len(vec_b)}")
    s = 0.0
    for i in range(len(vec_a)):
        d = float(vec_a[i]) - float(vec_b[i])
        s += d * d
    return math.sqrt(s)


def serf_predicted_rtt_seconds(coord_a: dict, coord_b: dict) -> tuple[float, float]:
    """
    Returns (d_vec, d_full) in seconds, following the Serf documentation:
      rtt = sqrt(sum diff^2) + Ha + Hb
      adjusted = rtt + Aa + Ab; if adjusted > 0 use adjusted.
    """
    d_vec = euclid(coord_a["Vec"], coord_b["Vec"])
    rtt = d_vec + float(coord_a["Height"]) + float(coord_b["Height"])
    adjusted = rtt + float(coord_a["Adjustment"]) + float(coord_b["Adjustment"])
    if adjusted > 0.0:
        rtt = adjusted
    return d_vec, rtt


def split_nodes(cell) -> list[str]:
    """TableB stores space-separated node lists; empty/NaN => []"""
    if cell is None or (isinstance(cell, float) and math.isnan(cell)):
        return []
    s = str(cell).strip()
    if not s:
        return []
    return s.split()


# =========================
# LOAD COORDINATES
# =========================
with open(CLUSTER_STATUS_JSON, "r") as f:
    js = json.load(f)

# Map both short and full names to coordinate dict
coord_map = {}
for n in js["nodes"]:
    full = n["name"]
    sn = short_name(full)
    coord = n["coordinate"]  # contains Vec, Height, Adjustment, Error
    coord_map[full] = coord
    coord_map[sn] = coord

# =========================
# LOAD TABLE B
# =========================
dfB = pd.read_csv(TABLEB_5D)

needed = {"query_node", "rtt_ms", "tp_nodes", "fp_nodes", "fn_nodes"}
missing = needed - set(dfB.columns)
if missing:
    raise ValueError(f"{TABLEB_5D} missing columns: {sorted(missing)}")

# normalize names + rtt
dfB["query_node"] = dfB["query_node"].astype(str).map(short_name)
dfB["rtt_ms"] = pd.to_numeric(dfB["rtt_ms"], errors="coerce")


# =========================
# VALIDATE
# =========================
rows = []
summary = []

for _, row in dfB.iterrows():
    q = row["query_node"]
    T_ms = float(row["rtt_ms"])

    if q not in coord_map:
        raise KeyError(f"Query node '{q}' not found in JSON coordinate map.")

    tp = split_nodes(row["tp_nodes"])
    fp = split_nodes(row["fp_nodes"])
    fn = split_nodes(row["fn_nodes"])

    # You can include FN in the validation too, but main focus is TP+FP
    check_nodes = [(n, "TP") for n in tp] + [(n, "FP") for n in fp] + [(n, "FN") for n in fn]

    fp_vec_ok_full_bad = 0
    fp_vec_bad = 0

    for n, cls in check_nodes:
        n = short_name(n)
        if n not in coord_map:
            # keep going but record missing nodes
            rows.append({
                "query": q, "T_ms": T_ms, "node": n, "class": cls,
                "d_vec_ms": None, "d_full_ms": None, "delta_ms": None,
                "reason": "node not found in JSON"
            })
            continue

        d_vec_s, d_full_s = serf_predicted_rtt_seconds(coord_map[q], coord_map[n])
        d_vec_ms = d_vec_s * 1000.0
        d_full_ms = d_full_s * 1000.0
        delta_ms = d_full_ms - d_vec_ms

        vec_in = (d_vec_ms <= T_ms + EPS_MS)
        full_in = (d_full_ms <= T_ms + EPS_MS)

        # Reason tagging (especially for FP)
        reason = ""
        if cls == "FP":
            if vec_in and not full_in:
                reason = "Vec OK, full RTT too high (height/adjust pushed over)"
                fp_vec_ok_full_bad += 1
            elif (not vec_in) and (not full_in):
                reason = "Vec already too high (not caused by height/adjust)"
                fp_vec_bad += 1
            elif vec_in and full_in:
                reason = "FP label inconsistent with computed full RTT (check inputs)"
            else:
                reason = "full in but vec out (rare; check adjustments/heights/units)"
        elif cls == "TP":
            if full_in:
                reason = "OK"
            else:
                reason = "TP label inconsistent with computed full RTT (check inputs)"
        elif cls == "FN":
            if full_in:
                reason = "Missed node (inside full RTT threshold)"
            else:
                reason = "FN label inconsistent with computed full RTT (check inputs)"

        rows.append({
            "query": q,
            "T_ms": T_ms,
            "node": n,
            "class": cls,
            "d_vec_ms": round(d_vec_ms, 6),
            "d_full_ms": round(d_full_ms, 6),
            "delta_ms": round(delta_ms, 6),
            "vec_in_T": vec_in,
            "full_in_T": full_in,
            "reason": reason
        })

    summary.append({
        "query": q,
        "T_ms": T_ms,
        "FP_count": len(fp),
        "FP_vecOK_fullBad": fp_vec_ok_full_bad,
        "FP_vecBad": fp_vec_bad
    })

out_df = pd.DataFrame(rows)
sum_df = pd.DataFrame(summary)

# Save outputs
out_df.to_csv("validate_vec_vs_full_details.csv", index=False)
sum_df.to_csv("validate_vec_vs_full_summary.csv", index=False)

print("Saved:")
print("  validate_vec_vs_full_details.csv  (per-node validation details)")
print("  validate_vec_vs_full_summary.csv  (FP reason breakdown per query)")
print("\nFP reason totals across all queries:")
print(sum_df[["FP_count","FP_vecOK_fullBad","FP_vecBad"]].sum().to_string())
