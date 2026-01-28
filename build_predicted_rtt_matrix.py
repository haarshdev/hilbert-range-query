import json
import math
import pandas as pd

# =========================
# INPUT / OUTPUT
# =========================
INPUT_JSON = "cluster-status/cluster-status15012026.json"
OUTPUT_CSV = "rtt-matrix/predicted/predicted_rtt_matrix_15012026.csv"

# =========================
# HELPERS
# =========================
def trunc3(x: float) -> float:
    # keep first 3 decimals, no rounding
    return math.floor(x * 1000.0) / 1000.0

def short_name(full_name: str) -> str:
    # clab-nebula-serf12 -> serf12
    if "serf" in full_name:
        return "serf" + full_name.split("serf")[-1]
    return full_name

# =========================
# LOAD CLUSTER STATUS
# =========================
with open(INPUT_JSON, "r") as f:
    data = json.load(f)

# Build name mapping
name_map = {}
for node in data["nodes"]:
    name_map[node["name"]] = short_name(node["name"])

# Sorted short names (serf1, serf2, ..., serf30)
short_nodes = sorted(
    name_map.values(),
    key=lambda x: int(x.replace("serf", ""))
)

# Initialize matrix
matrix = pd.DataFrame(index=short_nodes, columns=short_nodes, dtype=float)

# Diagonal
for n in short_nodes:
    matrix.loc[n, n] = 0.0

# Fill RTTs
for node in data["nodes"]:
    src_short = name_map[node["name"]]
    for dst_full, rtt in node.get("rtts", {}).items():
        dst_short = name_map[dst_full]
        matrix.loc[src_short, dst_short] = trunc3(float(rtt))

# Write CSV
matrix.to_csv(OUTPUT_CSV)
print("RTT matrix written to:", OUTPUT_CSV)
