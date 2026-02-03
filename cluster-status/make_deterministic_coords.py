import json, math, re

IN_FILE  = "cluster-status-15012026.json"
OUT_FILE = "cluster-status-deterministic.json"

def serf_index(name: str) -> int:
    m = re.search(r"serf(\d+)$", name)
    if not m:
        raise ValueError(f"Cannot parse serf index from name: {name}")
    return int(m.group(1))

def coord_for_idx(idx: int, cols: int = 6):
    # 6x5 grid for 30 nodes (idx 1..30)
    i = idx - 1
    row = i // cols   # 0..4
    col = i % cols    # 0..5

    # Step sizes in seconds
    x_step_s = 0.002000        # 2.000 ms
    y_step_s = 0.0043301270189 # 4.330 ms (chosen so max corner distance ~ 20 ms)

    x = col * x_step_s
    y = row * y_step_s
    return [x, y, 0.0, 0.0, 0.0]  # 5D


def dist_ms(vec_a, vec_b, ha=0.0, hb=0.0, aa=0.0, ab=0.0) -> float:
    # Same idea as Serf: Euclidean(Vec) + heights, then apply adjustments.
    sumsq = 0.0
    for va, vb in zip(vec_a, vec_b):
        d = va - vb
        sumsq += d * d
    rtt_s = math.sqrt(sumsq) + ha + hb
    adjusted = rtt_s + aa + ab
    if adjusted > 0.0:
        rtt_s = adjusted
    return round(rtt_s * 1000.0, 3)

with open(IN_FILE, "r") as f:
    data = json.load(f)

nodes = data["nodes"]

# Build deterministic Vec for each node
name_to_vec = {}
for n in nodes:
    idx = serf_index(n["name"])
    name_to_vec[n["name"]] = coord_for_idx(idx)

# Build output JSON
out = {
    "timestamp": data.get("timestamp"),
    "nodes": []
}

for n in nodes:
    name = n["name"]
    vec = name_to_vec[name]

    coord = {
        "Vec": vec,
        "Error": 0.0,
        "Adjustment": 0.0,
        "Height": 0.0
    }

    rtts = {}
    for m in nodes:
        other = m["name"]
        if other == name:
            continue
        rtts[other] = dist_ms(vec, name_to_vec[other])

    out_node = {
        "name": n["name"],
        "addr": n["addr"],
        "port": n["port"],
        "status": n["status"],
        "tags": n.get("tags", {}),
        "coordinate": coord,
        "rtts": rtts
    }
    out["nodes"].append(out_node)

with open(OUT_FILE, "w") as f:
    json.dump(out, f, indent=2)

print("Wrote:", OUT_FILE)
