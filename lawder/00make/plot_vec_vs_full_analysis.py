import pandas as pd
import matplotlib.pyplot as plt

# =========================
# INPUT
# =========================
CSV = "validate_vec_vs_full_details.csv"

# Colors (same as your document)
COLOR_TP = "#A8D08D"   # green
COLOR_FP = "#C00000"   # red
COLOR_FN = "#FFC000"   # orange

# =========================
# LOAD
# =========================
df = pd.read_csv(CSV)

# Keep only classified rows
df = df[df["class"].isin(["TP", "FP", "FN"])]

# =========================
# COUNT CAUSES
# =========================
tp_count = len(df[df["class"] == "TP"])

fp_vec_ok_full_bad = len(
    df[
        (df["class"] == "FP") &
        (df["vec_in_T"] == True) &
        (df["full_in_T"] == False)
    ]
)

fn_vec_bad_full_ok = len(
    df[
        (df["class"] == "FN") &
        (df["vec_in_T"] == False) &
        (df["full_in_T"] == True)
    ]
)

labels = [
    "TP\n(full RTT ≤ T)",
    "FP\n(Vec ≤ T,\nfull RTT > T)",
    "FN\n(Vec > T,\nfull RTT ≤ T)",
]

values = [
    tp_count,
    fp_vec_ok_full_bad,
    fn_vec_bad_full_ok,
]

colors = [COLOR_TP, COLOR_FP, COLOR_FN]

# =========================
# PLOT
# =========================
plt.figure(figsize=(6, 4))

plt.bar(labels, values, color=colors)

plt.ylabel("Number of nodes (summed over all queries)")
plt.title("Cause of TP / FP / FN in 5D vector-only queries")

plt.tight_layout()
plt.savefig("plot_fp_cause_summary.png", dpi=200)
plt.close()

print("Saved plot_fp_cause_summary.png")
print("Counts:")
print(" TP:", tp_count)
print(" FP (Vec OK, full RTT too high):", fp_vec_ok_full_bad)
print(" FN (Vec too high, full RTT OK):", fn_vec_bad_full_ok)
