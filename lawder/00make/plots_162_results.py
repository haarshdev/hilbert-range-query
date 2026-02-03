import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# =========================
# SETTINGS
# =========================
BASE_DIR = "."
ORDERS = [10, 11, 12, 13, 14, 15]

FILES = {
    "5D":  "tableA_summary_162_{order}_5D.csv",
    "5DP": "tableA_summary_162_{order}_5DP_new.csv",
}

# Table C (box debug) files (we will focus ONLY on orders 13–15)
TABLEC_FILES = {
    13: "tableC_boxes_162_13_5DP_debug_1.csv",
    14: "tableC_boxes_162_14_5DP_debug_1.csv",
    15: "tableC_boxes_162_15_5DP_debug_1.csv",
}

OUT_DIR = "plots_out"
os.makedirs(OUT_DIR, exist_ok=True)

# =========================
# COLORS (as requested)
# =========================
COLOR_TP   = "#A8D08D"
COLOR_FP   = "#C00000"
COLOR_FN   = "#FFC000"
COLOR_PREC = "#4F81BD"
COLOR_REC  = "#F79646"
COLOR_JAC  = "#9BBB59"

# Dumbbell specific
COLOR_LEFT = "#404040"   # optimistic LB used by pruning (dark gray)
COLOR_LINK = "#7A7A7A"   # connector line

# =========================
# STYLE (as requested)
# =========================
STYLE_5D = dict(
    linestyle=":",
    marker="+",
    markersize=9,
    linewidth=1.3,
    alpha=1.0,
    label="5D",
)

STYLE_5DP = dict(
    linestyle="-",
    marker="*",
    markersize=8,
    linewidth=1.3,
    alpha=0.5,
    label="5D+Pruning",
)

# =========================
# LOAD + AGGREGATE (micro)
# =========================
def load_and_aggregate(path: str) -> dict:
    df = pd.read_csv(path)
    req = ["truth_size", "returned", "tp", "fp", "fn"]

    for c in req:
        if c not in df.columns:
            raise ValueError(f"{path} missing column: {c}")
        df[c] = pd.to_numeric(df[c], errors="coerce").fillna(0.0)

    tp = float(df["tp"].sum())
    fp = float(df["fp"].sum())
    fn = float(df["fn"].sum())
    returned = float(df["returned"].sum())
    truth = float(df["truth_size"].sum())

    precision = tp / returned if returned > 0 else np.nan
    recall    = tp / truth if truth > 0 else np.nan
    jaccard   = tp / (tp + fp + fn) if (tp + fp + fn) > 0 else np.nan

    # FP density for pruning files (will be NaN if column not present)
    if "cover_accepted_leaf_points" in df.columns:
        cap = float(pd.to_numeric(df["cover_accepted_leaf_points"], errors="coerce").fillna(0.0).sum())
        fp_density = fp / cap if cap > 0 else np.nan
    else:
        fp_density = np.nan

    return {
        "tp": tp,
        "fp": fp,
        "fn": fn,
        "precision": precision,
        "recall": recall,
        "jaccard": jaccard,
        "fp_density": fp_density,
    }

def build_table():
    rows = []
    missing = []

    for order in ORDERS:
        for scenario, pat in FILES.items():
            fname = pat.format(order=order)
            path = os.path.join(BASE_DIR, fname)
            if not os.path.exists(path):
                missing.append(fname)
                continue

            agg = load_and_aggregate(path)
            agg.update({"order": order, "scenario": scenario})
            rows.append(agg)

    if not rows:
        raise SystemExit("No input files found.")

    if missing:
        print("Missing files (skipped):")
        for f in missing:
            print("  -", f)

    return pd.DataFrame(rows).sort_values(["order", "scenario"]).reset_index(drop=True)

def get_series(df, scenario, col, orders):
    return np.array(
        [float(df[(df.order == o) & (df.scenario == scenario)].iloc[0][col])
         for o in orders],
        dtype=float
    )

# =========================
# PLOTS (existing)
# =========================
def plot_fp_two_panels(df, outpath):
    orders = sorted(df["order"].unique())
    x = np.array(orders, dtype=int)

    fp_5d  = get_series(df, "5D",  "fp", orders)
    fp_5dp = get_series(df, "5DP", "fp", orders)

    fig, axes = plt.subplots(2, 1, figsize=(11, 7.0), sharex=True)

    axes[0].plot(x, fp_5d, color=COLOR_FP, **STYLE_5D)
    axes[0].set_title("False positives (5D)")
    axes[0].set_ylabel("FP count")
    axes[0].grid(True, axis="y", alpha=0.3)
    axes[0].legend(loc="best")

    axes[1].plot(x, fp_5dp, color=COLOR_FP, **STYLE_5DP)
    axes[1].set_title("False positives (5D+Pruning)")
    axes[1].set_xlabel("Hilbert order")
    axes[1].set_ylabel("FP count")
    axes[1].grid(True, axis="y", alpha=0.3)
    axes[1].legend(loc="best")

    fig.suptitle("False Positives across Hilbert orders", y=0.98)
    plt.tight_layout(rect=[0, 0, 1, 0.96])
    plt.savefig(outpath, dpi=200)
    plt.close()

def plot_fp_reduction_factor(df, outpath):
    orders = sorted(df["order"].unique())
    x = np.array(orders, dtype=int)

    fp_5d  = get_series(df, "5D",  "fp", orders)
    fp_5dp = get_series(df, "5DP", "fp", orders)
    factor = np.where(fp_5dp > 0, fp_5d / fp_5dp, np.nan)

    plt.figure(figsize=(10, 5.2))
    plt.plot(x, factor, linestyle="-", marker="o", linewidth=1.3, markersize=4, alpha=0.95)
    plt.xlabel("Hilbert order")
    plt.ylabel("FP reduction factor (FP_5D / FP_5DP)")
    plt.title("False Positive reduction due to pruning")
    plt.grid(True, axis="y", alpha=0.3)

    for xi, yi in zip(x, factor):
        if np.isfinite(yi):
            plt.text(xi, yi, f"{int(round(yi))}×", ha="center", va="bottom", fontsize=9)

    plt.tight_layout()
    plt.savefig(outpath, dpi=200)
    plt.close()

def plot_invariance_overlap(df, outpath):
    orders = sorted(df["order"].unique())
    x = np.array(orders, dtype=int)

    tp_5d  = get_series(df, "5D",  "tp", orders)
    tp_5dp = get_series(df, "5DP", "tp", orders)
    fn_5d  = get_series(df, "5D",  "fn", orders)
    fn_5dp = get_series(df, "5DP", "fn", orders)

    fig, axes = plt.subplots(2, 1, figsize=(11, 6.5), sharex=True)

    axes[0].plot(x, tp_5d,  color=COLOR_TP, **STYLE_5D)
    axes[0].plot(x, tp_5dp, color=COLOR_TP, **STYLE_5DP)
    axes[0].set_title("True Positive (TP)")
    axes[0].set_ylabel("TP count")
    axes[0].grid(True, axis="y", alpha=0.3)
    axes[0].legend(loc="best")

    axes[1].plot(x, fn_5d,  color=COLOR_FN, **STYLE_5D)
    axes[1].plot(x, fn_5dp, color=COLOR_FN, **STYLE_5DP)
    axes[1].set_title("False Negative (FN)")
    axes[1].set_xlabel("Hilbert order")
    axes[1].set_ylabel("FN count")
    axes[1].grid(True, axis="y", alpha=0.3)
    axes[1].legend(loc="best")

    fig.suptitle("True Positives and False Negatives across Hilbert orders", y=0.98)
    plt.tight_layout(rect=[0, 0, 1, 0.96])
    plt.savefig(outpath, dpi=200)
    plt.close()

def plot_metrics_three_panels_markers(df, outpath):
    orders = sorted(df["order"].unique())
    x = np.array(orders, dtype=int)

    p_5d  = get_series(df, "5D",  "precision", orders)
    p_5dp = get_series(df, "5DP", "precision", orders)

    j_5d  = get_series(df, "5D",  "jaccard", orders)
    j_5dp = get_series(df, "5DP", "jaccard", orders)

    r_5d  = get_series(df, "5D",  "recall", orders)
    r_5dp = get_series(df, "5DP", "recall", orders)

    fig, axes = plt.subplots(3, 1, figsize=(11, 8.2), sharex=True)

    def draw(ax, y5d, y5dp, title, ylabel, color, show_legend=False):
        ax.plot(x, y5d,  color=color, **STYLE_5D)
        ax.plot(x, y5dp, color=color, **STYLE_5DP)
        ax.set_title(title)
        ax.set_ylabel(ylabel)
        ax.set_ylim(0, 1.05)
        ax.grid(True, axis="y", alpha=0.3)
        if show_legend:
            ax.legend(loc="lower left")

    draw(axes[0], p_5d, p_5dp, "Precision", "Precision", COLOR_PREC, show_legend=True)
    draw(axes[1], r_5d, r_5dp, "Recall", "Recall", COLOR_REC, show_legend=True)
    draw(axes[2], j_5d, j_5dp, "Jaccard", "Jaccard", COLOR_JAC, show_legend=True)

    axes[2].set_xlabel("Hilbert order")
    axes[2].set_xticks(x)
    axes[2].set_xticklabels([str(o) for o in orders])

    fig.suptitle("Precision, Recall, Jaccard across Hilbert orders: 5D vs 5D+Pruning", y=0.98)
    plt.tight_layout(rect=[0, 0, 1, 0.96])
    plt.savefig(outpath, dpi=200)
    plt.close()

def plot_fp_density_pruning(df, outpath):
    orders = sorted(df["order"].unique())
    x = np.array(orders, dtype=int)

    fp_density = np.array(
        [float(df[(df.order == o) & (df.scenario == "5DP")].iloc[0]["fp_density"])
         for o in orders],
        dtype=float
    )

    fp_pct = fp_density * 100.0

    plt.figure(figsize=(9, 4.8))
    plt.bar(x, fp_pct, width=0.6, color=COLOR_FP)

    plt.xlabel("Hilbert order")
    plt.ylabel("False positive density (%)")
    plt.title("False positive density among points surviving pruning")
    plt.xticks(x)
    plt.ylim(0, max(fp_pct) * 1.25 if len(fp_pct) else 1)
    plt.grid(True, axis="y", alpha=0.3)

    for xi, yi in zip(x, fp_pct):
        plt.text(xi, yi, f"{yi:.2f}%", ha="center", va="bottom", fontsize=10)

    plt.tight_layout()
    plt.savefig(outpath, dpi=200)
    plt.close()

# =========================
# Q4: DUMBBELL PLOT (orders 13–15 only)
# =========================
def load_tablec_orders_13_15():
    rows = []
    missing = []

    for order, fname in TABLEC_FILES.items():
        path = os.path.join(BASE_DIR, fname)
        if not os.path.exists(path):
            missing.append(fname)
            continue

        df = pd.read_csv(path)

        needed = ["rtt_ms", "fp_node", "adjusted_lb_ms", "serf_effective_lb_ms"]
        for c in needed:
            if c not in df.columns:
                raise ValueError(f"{fname} missing column: {c}")

        # Coerce numerics
        df["rtt_ms"] = pd.to_numeric(df["rtt_ms"], errors="coerce")
        df["adjusted_lb_ms"] = pd.to_numeric(df["adjusted_lb_ms"], errors="coerce")
        df["serf_effective_lb_ms"] = pd.to_numeric(df["serf_effective_lb_ms"], errors="coerce")

        # Determine order robustly
        if "hilbert_order" in df.columns:
            df["hilbert_order"] = pd.to_numeric(df["hilbert_order"], errors="coerce").fillna(order).astype(int)
        else:
            df["hilbert_order"] = int(order)

        df["fp_node"] = df["fp_node"].astype(str)
        df["order_file"] = int(order)

        rows.append(df)

    if not rows:
        raise SystemExit("No Table C files found for orders 13–15.")

    if missing:
        print("Missing Table C files (skipped):")
        for f in missing:
            print("  -", f)

    out = pd.concat(rows, ignore_index=True)
    out = out[out["hilbert_order"].isin([13, 14, 15])].copy()

    # Drop invalid rows
    out = out[np.isfinite(out["rtt_ms"]) &
              np.isfinite(out["adjusted_lb_ms"]) &
              np.isfinite(out["serf_effective_lb_ms"])].copy()

    return out

def plot_q4_guard_mismatch_dumbbell(tablec_all, outpath, prefer_order=13, label_all=True):
    """
    One row per FP node, showing the mismatch:

      left dot  = adjusted_lb_ms (optimistic "best-case" LB used by pruning)
      right dot = serf_effective_lb_ms (effective LB under Serf/Vivaldi guard)

    Interpretation:
      - If left dot is <= T, the box survives pruning.
      - If right dot is  > T, then under Serf semantics it should be prunable.
      - This mismatch explains persistent FPs at orders 13–15.

    We plot a single order (default 13) because you already verified the FP set is the same for 13–15.
    """

    df = tablec_all.copy()
    if prefer_order in set(df["hilbert_order"].unique()):
        df = df[df["hilbert_order"] == prefer_order].copy()
    else:
        df = df[df["hilbert_order"] == int(df["hilbert_order"].min())].copy()

    # Threshold (expected constant, typically 5ms in these debug files)
    Ts = sorted(set(float(x) for x in df["rtt_ms"].dropna().unique()))
    T = float(Ts[0]) if Ts else 5.0

    # One row per node: keep the most optimistic adjusted LB (smallest)
    df_nodes = (
        df.sort_values("adjusted_lb_ms", ascending=True)
          .groupby("fp_node", as_index=False)
          .first()
          .copy()
    )

    # Sort for readability (largest guarded RTT at top)
    df_nodes = df_nodes.sort_values("serf_effective_lb_ms", ascending=False).reset_index(drop=True)

    nodes = df_nodes["fp_node"].tolist()
    y = np.arange(len(nodes))

    x_left = df_nodes["adjusted_lb_ms"].to_numpy(float)
    x_right = df_nodes["serf_effective_lb_ms"].to_numpy(float)

    # Tight x-range
    lo = min(np.min(x_left), T) - 0.8
    hi = max(np.max(x_right), T) + 0.8

    plt.figure(figsize=(12, max(5.5, 0.36 * len(nodes) + 1.6)))

    # Connectors
    for yi, xl, xr in zip(y, x_left, x_right):
        plt.plot([xl, xr], [yi, yi], linewidth=2.0, alpha=0.60, color=COLOR_LINK)

    # Dots
    plt.scatter(
        x_left, y, s=90, color=COLOR_LEFT, alpha=0.95,
        label="Pruning lower bound (most optimistic adjustment available)"
    )
    plt.scatter(
        x_right, y, s=90, color=COLOR_FP, alpha=0.95,
        label="RTT computed by Serf (with adjustment guard)"
    )

    # Threshold line
    plt.axvline(T, linewidth=1.7, alpha=0.90)
    plt.text(T + 0.06, -0.8, f"T = {T:.0f} ms", fontsize=9)

    # Labels and formatting
    plt.yticks(y, nodes if label_all else [""] * len(nodes))
    plt.xlabel("Lower bound RTT (ms)")
    plt.title("Why false positives persist at orders 13–15: adjustment-guard mismatch")
    plt.grid(True, axis="x", alpha=0.25)
    plt.legend(loc="lower right")

    plt.xlim(lo, hi)
    plt.tight_layout()
    plt.savefig(outpath, dpi=200)
    plt.close()

# =========================
# CONSOLE OUTPUT
# =========================
def print_console_summary(df):
    print("\n=== Aggregated results (micro-averaged) ===\n")

    orders_present = sorted(df["order"].unique())
    for order in orders_present:
        print(f"Hilbert order {order}")
        for scenario in ["5D", "5DP"]:
            sub = df[(df.order == order) & (df.scenario == scenario)]
            if sub.empty:
                continue
            row = sub.iloc[0]
            print(
                f"  {scenario:3s} | "
                f"TP={int(row.tp):5d}  "
                f"FP={int(row.fp):5d}  "
                f"FN={int(row.fn):3d}  "
                f"P={row.precision:.4f}  "
                f"R={row.recall:.4f}  "
                f"J={row.jaccard:.4f}"
            )
        print()

    print("=== Compact table ===")
    print(df[["order", "scenario", "tp", "fp", "fn", "precision", "recall", "jaccard"]]
          .to_string(index=False))

# =========================
# MAIN
# =========================
if __name__ == "__main__":
    df = build_table()

    df.to_csv(
        os.path.join(OUT_DIR, "aggregate_summary_orders_10_15_5D_vs_5DP.csv"),
        index=False
    )

    print_console_summary(df)

    plot_fp_two_panels(
        df,
        os.path.join(OUT_DIR, "fig_fp_two_panels_orders_10_15.png")
    )

    plot_fp_reduction_factor(
        df,
        os.path.join(OUT_DIR, "fig_fp_reduction_factor_orders_10_15.png")
    )

    plot_invariance_overlap(
        df,
        os.path.join(OUT_DIR, "fig_invariance_tp_fn_overlap_orders_10_15.png")
    )

    plot_metrics_three_panels_markers(
        df,
        os.path.join(OUT_DIR, "fig_metrics_three_panels_orders_10_15.png")
    )

    plot_fp_density_pruning(
        df,
        os.path.join(OUT_DIR, "fig_fp_density_pruning_orders_10_15.png")
    )

    # Q4: Dumbbell plot from Table C (orders 13–15 only)
    tablec_all = load_tablec_orders_13_15()
    tablec_all.to_csv(os.path.join(OUT_DIR, "q4_tablec_orders_13_15_merged.csv"), index=False)

    plot_q4_guard_mismatch_dumbbell(
        tablec_all,
        os.path.join(OUT_DIR, "fig_q4_dumbbell_guard_mismatch_orders_13_15.png"),
        prefer_order=13,
        label_all=True
    )

    print("\nSaved plots:")
    print("  - fig_fp_two_panels_orders_10_15.png")
    print("  - fig_fp_reduction_factor_orders_10_15.png")
    print("  - fig_invariance_tp_fn_overlap_orders_10_15.png")
    print("  - fig_metrics_three_panels_orders_10_15.png")
    print("  - fig_fp_density_pruning_orders_10_15.png")
    print("  - fig_q4_dumbbell_guard_mismatch_orders_13_15.png")
    print("  - aggregate_summary_orders_10_15_5D_vs_5DP.csv")
    print("  - q4_tablec_orders_13_15_merged.csv")
