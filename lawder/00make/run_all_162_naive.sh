#!/usr/bin/env bash
set -euo pipefail

# -------------------------
# SETTINGS (edit these)
# -------------------------
JSON_PATH="cluster-status-2025-07-11-1-rtts_recomputed.json"   # 162-node JSON
DRIVER_BIN="./serf_driver_naive.exe"                          # compiled naive binary
RTTS=(5 15)

OUT_SUMMARY="tableA_naive_summary_162.csv"
OUT_NODES="tableB_naive_nodes_162.csv"   # optional, comment out if you don't want it

# -------------------------
# SAFETY: start fresh outputs
# -------------------------
rm -f "$OUT_SUMMARY" "$OUT_NODES"

# -------------------------
# Extract node names from JSON (supports root array OR {"nodes":[...]} )
# -------------------------
mapfile -t NODES < <(python3 - <<PY
import json

path = "$JSON_PATH"
with open(path, "r", encoding="utf-8") as f:
    root = json.load(f)

if isinstance(root, list):
    nodes = root
elif isinstance(root, dict) and isinstance(root.get("nodes"), list):
    nodes = root["nodes"]
else:
    raise SystemExit("ERROR: JSON must be an array OR an object with a 'nodes' array.")

for n in nodes:
    name = n.get("name")
    if isinstance(name, str) and name:
        print(name)
PY
)

echo "Found ${#NODES[@]} nodes in JSON."
echo "Running naive: rtts={${RTTS[*]}} for every node."

# -------------------------
# Run driver for each node and RTT
# -------------------------
for q in "${NODES[@]}"; do
  for T in "${RTTS[@]}"; do
    tag="$(echo "$q" | sed 's/[^a-zA-Z0-9._-]/_/g')"
    log="${LOG_DIR}/${tag}_T${T}.log"

    echo "QNODE=$q  RTT=$T"

    # If you do NOT want the nodes CSV, remove the --out_naive_nodes_csv line
    "$DRIVER_BIN" \
      --json "$JSON_PATH" \
      --qnode "$q" \
      --rtt "$T" \
      --out_naive_summary_csv "$OUT_SUMMARY" \
      #--out_naive_nodes_csv "$OUT_NODES" \
      >"$log" 2>&1
  done
done

echo "Done."
echo "Saved:"
echo "  $OUT_SUMMARY"
echo "  $OUT_NODES"
