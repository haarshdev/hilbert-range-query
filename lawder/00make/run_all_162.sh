#!/usr/bin/env bash
set -euo pipefail

# -------------------------
# SETTINGS (edit these)
# -------------------------
JSON_PATH="cluster-status-2025-07-11-1-rtts_recomputed.json"   # your 162-node JSON cluster-status-2025-07-11-1.json
DRIVER_BIN="./serf_driver.exe"                    # compiled binary path
HORDER="12"
RTTS=(5 15)

OUT_SUMMARY="tableA_summary_162_12_5DP_debug_1.csv"
OUT_NODES="tableB_nodes_162_12_5DP_debug_1.csv"
OUT_BOXES="tableC_boxes_162_12_5DP_debug_1.csv"   # NEW: Table C box debug
FP_COUNTS="fp_counts_h12_5DP_debug_1.json"                # optional but useful

LOG_DIR="run_logs_162_h12_5DP_debug_1"
mkdir -p "$LOG_DIR"

# -------------------------
# SAFETY: start fresh outputs
# -------------------------
rm -f "$OUT_SUMMARY" "$OUT_NODES" "$FP_COUNTS"

# -------------------------
# Extract node names from JSON (supports root array OR {"nodes":[...]} )
# -------------------------
mapfile -t NODES < <(python3 - <<'PY'
import json

path = "cluster-status-2025-07-11-1.json"
with open(path, "r", encoding="utf-8") as f:
    root = json.load(f)

if isinstance(root, list):
    nodes = root
elif isinstance(root, dict) and isinstance(root.get("nodes"), list):
    nodes = root["nodes"]
else:
    raise SystemExit("ERROR: JSON must be an array OR an object with a 'nodes' array.")

names = []
for n in nodes:
    name = n.get("name")
    if isinstance(name, str) and name:
        names.append(name)

# stable order, but do not sort if you want original order
for name in names:
    print(name)
PY
)

echo "Found ${#NODES[@]} nodes in JSON."
echo "Running: horder=$HORDER, rtts={${RTTS[*]}} for every node."

# -------------------------
# Run driver for each node and RTT
# -------------------------
for q in "${NODES[@]}"; do
  for T in "${RTTS[@]}"; do
    tag="$(echo "$q" | sed 's/[^a-zA-Z0-9._-]/_/g')"
    log="${LOG_DIR}/${tag}_T${T}_H${HORDER}.log"

    echo "QNODE=$q  RTT=$T  HORDER=$HORDER"
    "$DRIVER_BIN" \
      --json "$JSON_PATH" \
      --qnode "$q" \
      --rtt "$T" \
      --horder "$HORDER" \
      --out_summary_csv "$OUT_SUMMARY" \
      --out_nodes_csv "$OUT_NODES" \
      --out_boxes_csv "$OUT_BOXES" \
      >"$log" 2>&1
  done
done

echo "Done."
echo "Saved:"
echo "  $OUT_SUMMARY"
echo "  $OUT_NODES"
echo "  $OUT_BOXES"
echo "  $FP_COUNTS"
echo "  Logs in: $LOG_DIR/"
