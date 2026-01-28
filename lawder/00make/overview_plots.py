import os
import math
import pandas as pd
import matplotlib.pyplot as plt


# =========================
# INPUT FILES (edit if needed)
# =========================
SCENARIOS = {
    "5D": {
        "tableA": "tableA_summary_5D.csv",
        "tableB": "tableB_nodes_5D.csv",
    },
    "6D": {
        "tableA": "tableA_summary_6D.csv",
        "tableB": "tableB_nodes_6D.csv",
    },
    "5D+Prune": {
        "tableA": "tableA_summary_5D_pruning.csv",
        "tableB": "tableB_nodes_5D_pruning.csv",
    },
}

OUT_DIR = "plots_overview_order10"
os.makedirs(OUT_DIR, exist_ok=True)

# =========================
# COLORS (user requested)
# =========================
COLOR_TP = "#A8D08D"
COLOR_FP = "#C00000"
COLOR_FN = "#FFC000"

# For metric bars (keep simple and readable)
COLOR_PREC = "#4F81BD"
COLOR_REC = "#F79646"
COLOR_JAC = "#9BBB59"


# =========================
# HELPERS
# =========================
def safe_div(a: float, b: float) -> float:
    return float(a) / float(b) if b else float("nan")


def parse_nodes(cell) -> list[str]:
    """TableB stores space-separated node names; NaN means empty."""
    if cell is None or (isinstance(cell, float) and math.isnan(cell)):
        return []
    s = str(cell).strip()
    if not s:
        return []
    return s.split()


def micro_from_counts(tp: int, fp: int, fn: int) -> dict:
    precision = safe_div(tp, tp + fp)
    recall = safe_div(tp, tp + fn)
    jaccard = safe_div(tp, tp + fp + fn)
    return {"precision": precision, "recall": recall, "jaccard": jaccard}


def read_tableA(path: str) -> pd.DataFrame:
    df = pd.read_csv(path)
    needed = {"query_node", "rtt_ms", "tp", "fp", "fn"}
    missing = needed - set(df.columns)
    if missing:
        raise ValueError(f"{path} missing columns: {sorted(missing)}")

    # normalize rtt_ms to int if possible
    df["rtt_ms"] = pd.to_numeric(df["rtt_ms"], errors="coerce").astype("Int64")
    return df


def read_tableB(path: str) -> pd.DataFrame:
    df = pd.read_csv(path)
    needed = {"query_node", "rtt_ms", "fp_nodes"}
    missing = needed - set(df.columns)
    if missing:
        raise ValueError(f"{path} missing columns: {sorted(missing)}")

    df["rtt_ms"] = pd.to_numeric(df["rtt_ms"], errors="coerce").astype("Int64")
    return df


# =========================
# LOAD + AGGREGATE
# =========================
scenario_totals = []          # per-scenario TP/FP/FN totals
scenario_metrics = []         # per-scenario micro metrics
scenario_fp_by_T = []         # per-scenario FP split by threshold
scenario_fp_node_counts = {}  # per-scenario node -> FP frequency

for scen, files in SCENARIOS.items():
    dfA = read_tableA(files["tableA"])

    tp = int(dfA["tp"].sum())
    fp = int(dfA["fp"].sum())
    fn = int(dfA["fn"].sum())

    scenario_totals.append({"scenario": scen, "tp": tp, "fp": fp, "fn": fn})

    m = micro_from_counts(tp, fp, fn)
    scenario_metrics.append({"scenario": scen, **m})

    # FP split by threshold (5 ms vs 15 ms)
    fp5 = int(dfA.loc[dfA["rtt_ms"] == 5, "fp"].sum())
    fp15 = int(dfA.loc[dfA["rtt_ms"] == 15, "fp"].sum())
    scenario_fp_by_T.append({"scenario": scen, "fp_5ms": fp5, "fp_15ms": fp15})

    # FP per node frequency from TableB
    dfB = read_tableB(files["tableB"])
    counts = {}
    for cell in dfB["fp_nodes"].tolist():
        for n in parse_nodes(cell):
            counts[n] = counts.get(n, 0) + 1
    scenario_fp_node_counts[scen] = counts

totals_df = pd.DataFrame(scenario_totals)
metrics_df = pd.DataFrame(scenario_metrics)
fpT_df = pd.DataFrame(scenario_fp_by_T)


# =========================
# PLOT 1: TP/FP/FN grouped bars per scenario (NOT stacked)
# =========================
plt.figure(figsize=(8, 4.5))

x = list(range(len(totals_df)))
w = 0.25

tp_vals = totals_df["tp"].tolist()
fp_vals = totals_df["fp"].tolist()
fn_vals = totals_df["fn"].tolist()

plt.bar([i - w for i in x], tp_vals, width=w, label="TP", color=COLOR_TP)
plt.bar(x, fp_vals, width=w, label="FP", color=COLOR_FP)
plt.bar([i + w for i in x], fn_vals, width=w, label="FN", color=COLOR_FN)

plt.xticks(x, totals_df["scenario"].tolist())
plt.ylabel("Total count (sum over 18 queries)")
plt.title("Total TP / FP / FN per scenario (Order 10, T=5 and 15 ms)")
plt.legend()
plt.tight_layout()
plt.savefig(os.path.join(OUT_DIR, "plot1_counts_grouped.png"), dpi=200)
plt.close()


# =========================
# PLOT 2: Micro Precision/Recall/Jaccard per scenario
# =========================
plt.figure(figsize=(8, 4.5))

x = list(range(len(metrics_df)))
w = 0.25

prec = metrics_df["precision"].tolist()
rec = metrics_df["recall"].tolist()
jac = metrics_df["jaccard"].tolist()

plt.bar([i - w for i in x], prec, width=w, label="Precision", color=COLOR_PREC)
plt.bar(x, rec, width=w, label="Recall", color=COLOR_REC)
plt.bar([i + w for i in x], jac, width=w, label="Jaccard", color=COLOR_JAC)

plt.xticks(x, metrics_df["scenario"].tolist())
plt.ylim(0, 1.05)
plt.ylabel("Value")
plt.title("Precision, Recall, Jaccard per scenario (Order 10, 18 queries)")
plt.legend()
plt.tight_layout()
plt.savefig(os.path.join(OUT_DIR, "plot2_metrics_micro.png"), dpi=200)
plt.close()


# =========================
# PLOT 3: FP split by RTT threshold per scenario
# =========================
plt.figure(figsize=(8, 4.5))

x = list(range(len(fpT_df)))
w = 0.35

plt.bar([i - w/2 for i in x], fpT_df["fp_5ms"].tolist(), width=w, label="FP @ 5 ms", color=COLOR_FP)
# use a slightly lighter red for contrast but same family
plt.bar([i + w/2 for i in x], fpT_df["fp_15ms"].tolist(), width=w, label="FP @ 15 ms", color="#E06666")

plt.xticks(x, fpT_df["scenario"].tolist())
plt.ylabel("FP count (sum over queries at that threshold)")
plt.title("False positives split by threshold (Order 10)")
plt.legend()
plt.tight_layout()
plt.savefig(os.path.join(OUT_DIR, "plot3_fp_by_threshold.png"), dpi=200)
plt.close()


# =========================
# PLOT 4 (optional): FP frequency per node (one plot per scenario)
# =========================
for scen, counts in scenario_fp_node_counts.items():
    out_path = os.path.join(OUT_DIR, f"plot4_fp_node_frequency_{scen}.png")

    if not counts:
        plt.figure(figsize=(7, 3.5))
        plt.title(f"FP frequency per node ({scen})")
        plt.xlabel("Node")
        plt.ylabel("FP count across all queries")
        plt.text(0.5, 0.5, "No false positives", ha="center", va="center")
        plt.tight_layout()
        plt.savefig(out_path, dpi=200)
        plt.close()
        continue

    items = sorted(counts.items(), key=lambda kv: (-kv[1], kv[0]))
    nodes = [k for k, _ in items]
    vals = [v for _, v in items]

    plt.figure(figsize=(max(8, len(nodes) * 0.35), 4))
    plt.bar(range(len(nodes)), vals, color=COLOR_FP)
    plt.xticks(range(len(nodes)), nodes, rotation=90)
    plt.ylabel("FP count across all queries")
    plt.title(f"FP frequency per node ({scen})")
    plt.tight_layout()
    plt.savefig(out_path, dpi=200)
    plt.close()


print("Done. Plots saved in:", OUT_DIR)
print("\nTotals:")
print(totals_df.to_string(index=False))
print("\nMicro metrics:")
print(metrics_df.to_string(index=False))
print("\nFP split by threshold:")
print(fpT_df.to_string(index=False))
