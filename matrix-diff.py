import pandas as pd

# =========================
# SETTINGS
# =========================
GT_FILE = "rtt-matrix/gt/ground_truth_min_rtt_matrix.csv"
PRED_FILE = "rtt-matrix/predicted/predicted_rtt_matrix.csv"
OUT_FILE = "rtt-matrix/rtt_differences_gt_vs_pred.csv"

# Only consider differences >= 1 ms (ignore anything below 1 ms)
DIFF_THRESHOLD_MS = 1.0

# =========================
# LOAD
# =========================
gt = pd.read_csv(GT_FILE, index_col=0)
pred = pd.read_csv(PRED_FILE, index_col=0)

# Keep only common nodes (in case one file has extra/missing)
common = sorted(set(gt.index) & set(pred.index) & set(gt.columns) & set(pred.columns))
gt = gt.loc[common, common]
pred = pred.loc[common, common]

rows = []
for src in common:
    for dst in common:
        if src == dst:
            continue

        gt_rtt = float(gt.loc[src, dst])
        pr_rtt = float(pred.loc[src, dst])

        diff = pr_rtt - gt_rtt  # + means predicted is higher, - means predicted is lower

        # Ignore small differences (< 1ms), like 25.180 vs 25.883
        if abs(diff) < DIFF_THRESHOLD_MS:
            continue

        rows.append({
            "from": src,
            "to": dst,
            "gt_rtt_ms": round(gt_rtt, 3),
            "pred_rtt_ms": round(pr_rtt, 3),
            "diff_ms": round(diff, 3),  # keep 3 decimals, with sign
            "direction": "OVER" if diff > 0 else "UNDER"
        })

df = pd.DataFrame(rows)

# Sort so biggest differences are easy to spot
if not df.empty:
    df["abs_diff_ms"] = df["diff_ms"].abs()
    df = df.sort_values(["abs_diff_ms"], ascending=False).drop(columns=["abs_diff_ms"])

# Save
df.to_csv(OUT_FILE, index=False, float_format="%.3f")
print("Wrote:", OUT_FILE)

# Show top 20 on screen
print("\nTop 20 biggest differences (>= 1.000 ms only):")
if df.empty:
    print("No differences >= 1.000 ms found.")
else:
    print(df.head(20).to_string(index=False))
