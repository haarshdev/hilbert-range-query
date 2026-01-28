import argparse
import json
import math
import csv
import statistics
from typing import Dict, Any, List, Tuple


def short_name(full: str) -> str:
    if "-" not in full:
        return full
    return full.rsplit("-", 1)[-1]


def guess_scale(nodes: List[Dict[str, Any]]) -> float:
    """
    Guess whether coords are in seconds (need *1000) or already ms.
    Heuristic: if median abs(Vec component) < 1.0 => seconds.
    """
    vals: List[float] = []
    for n in nodes:
        coord = n.get("coordinate") or n.get("Coord") or {}
        vec = coord.get("Vec")
        if not isinstance(vec, list) or not vec:
            continue
        for x in vec[:5]:
            try:
                vals.append(abs(float(x)))
            except Exception:
                pass
    if not vals:
        return 1000.0
    med = statistics.median(vals)
    return 1000.0 if med < 1.0 else 1.0


def load_nodes(json_path: str) -> Tuple[List[Dict[str, Any]], Dict[str, str]]:
    with open(json_path, "r", encoding="utf-8") as f:
        root = json.load(f)

    if isinstance(root, list):
        nodes = root
    elif isinstance(root, dict) and isinstance(root.get("nodes"), list):
        nodes = root["nodes"]
    else:
        raise ValueError("JSON must be an array or an object with a 'nodes' array")

    short_to_full: Dict[str, str] = {}
    for n in nodes:
        name = n.get("name")
        if isinstance(name, str):
            short_to_full[short_name(name)] = name

    return nodes, short_to_full


def get_coord_fields_ms(n: Dict[str, Any], scale: float) -> Tuple[List[float], float, float]:
    coord = n.get("coordinate") or n.get("Coord") or {}
    vec = coord.get("Vec") or []
    height = coord.get("Height", 0.0)
    adj = coord.get("Adjustment", 0.0)

    vec5 = []
    for i in range(min(5, len(vec))):
        vec5.append(float(vec[i]) * scale)

    height_ms = float(height) * scale
    adj_ms = float(adj) * scale
    return vec5, height_ms, adj_ms


def vec_distance_ms(a_vec5: List[float], b_vec5: List[float]) -> float:
    m = min(len(a_vec5), len(b_vec5))
    s = 0.0
    for i in range(m):
        d = a_vec5[i] - b_vec5[i]
        s += d * d
    return math.sqrt(s)


def serf_full_rtt_ms(vec_only: float, base_vec_plus_h: float, adj_sum: float) -> float:
    """
    Serf logic:
      rtt = vec_dist + height_a + height_b
      adjusted = rtt + adj_a + adj_b
      if adjusted > 0: rtt = adjusted
    """
    adjusted = base_vec_plus_h + adj_sum
    return adjusted if adjusted > 0.0 else base_vec_plus_h


def parse_nodes_list(cell: str) -> List[str]:
    cell = (cell or "").strip()
    if not cell:
        return []
    return [x for x in cell.split() if x]


def classify(node_type: str, vec_only: float, full: float, T: float) -> str:
    # node_type is "TP" or "FP" just for readable labels
    if node_type == "FP":
        if vec_only <= T and full > T:
            return "PUSHED_OVER_BY_HEIGHT_ADJ"
        if vec_only > T:
            return "ALREADY_FAR_IN_VEC"
        return "OTHER_OR_EDGE"
    else:
        # TP
        if vec_only <= T and full <= T:
            return "VEC_ALONE_SUFFICIENT"
        if vec_only > T and full <= T:
            return "PULLED_IN_BY_HEIGHT_ADJ"
        if vec_only <= T and full > T:
            # should not happen if truth and full align, but keep for sanity
            return "INCONSISTENT"
        return "OTHER_OR_EDGE"


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--json", default="cluster-status-15012026.json")
    ap.add_argument("--tableB", default="tableB_nodes.csv")
    ap.add_argument("--rtt", type=float, default=None, help="Filter to a specific RTT threshold (e.g., 15)")
    ap.add_argument("--query", default=None, help="Filter to one query node short name (e.g., serf1)")
    ap.add_argument("--scale", type=float, default=None, help="Force scale (1000 for seconds->ms, 1 for already-ms)")
    ap.add_argument("--show", choices=["fp", "tp", "both"], default="both")
    args = ap.parse_args()

    nodes, short_to_full = load_nodes(args.json)
    full_to_node: Dict[str, Dict[str, Any]] = {n["name"]: n for n in nodes if isinstance(n.get("name"), str)}

    scale = args.scale if args.scale is not None else guess_scale(nodes)
    print(f"Scale used: {scale} (1000 means coords are seconds; 1 means coords are ms)\n")

    with open(args.tableB, "r", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        rows = list(reader)

    for row in rows:
        q_short = row.get("query_node", "")
        rtt_ms = float(row.get("rtt_ms", "nan"))

        if args.rtt is not None and abs(rtt_ms - args.rtt) > 1e-9:
            continue
        if args.query is not None and q_short != args.query:
            continue

        q_full = short_to_full.get(q_short, q_short)
        if q_full not in full_to_node:
            print(f"Skip {q_short}: not found in JSON names")
            continue

        fp_list = parse_nodes_list(row.get("fp_nodes", ""))
        tp_list = parse_nodes_list(row.get("tp_nodes", ""))
        T = rtt_ms

        qnode = full_to_node[q_full]
        q_vec, q_h, q_adj = get_coord_fields_ms(qnode, scale)

        rtts_obj = qnode.get("rtts") or {}
        if not isinstance(rtts_obj, dict):
            rtts_obj = {}

        print(f"=== Query {q_short} RTT {T:.1f} ms | TP {len(tp_list)} | FP {len(fp_list)} ===")

        def analyze_list(node_type: str, lst: List[str]) -> None:
            if not lst:
                return

            counts: Dict[str, int] = {}
            for s in lst:
                full = short_to_full.get(s, s)
                if full not in full_to_node:
                    counts["MISSING_IN_JSON"] = counts.get("MISSING_IN_JSON", 0) + 1
                    continue

                n = full_to_node[full]
                v, h, adj = get_coord_fields_ms(n, scale)

                vec_only = vec_distance_ms(q_vec, v)
                base = vec_only + q_h + h
                full_rtt = serf_full_rtt_ms(vec_only, base, q_adj + adj)

                truth = None
                if full in rtts_obj:
                    try:
                        truth = float(rtts_obj[full])
                    except Exception:
                        truth = None

                label = classify(node_type, vec_only, full_rtt, T)
                counts[label] = counts.get(label, 0) + 1

                truth_s = f"{truth:.3f}" if truth is not None else "NA"
                print(
                    f"  {node_type} {s:>6} | truth={truth_s:>8} | "
                    f"vec_only={vec_only:>7.3f} | vec+H={base:>7.3f} | full={full_rtt:>7.3f} | {label}"
                )

            print(f"{node_type} summary:")
            for k in sorted(counts.keys()):
                print(f"  {k}: {counts[k]}")
            print()

        if args.show in ("fp", "both"):
            analyze_list("FP", sorted(fp_list))
        if args.show in ("tp", "both"):
            analyze_list("TP", sorted(tp_list))

        print()

if __name__ == "__main__":
    main()
