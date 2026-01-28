#!/usr/bin/env python3
import subprocess
import re
import csv
import time
from pathlib import Path

# ----------------------------
# CONFIG: set these
# ----------------------------
LAB = "nebula"  # must match containerlab topology name (e.g., "century")
PING_COUNT = 3         # number of ICMP probes per pair
PING_TIMEOUT_S = 2     # per-probe timeout
SLEEP_BETWEEN_PINGS_S = 0.0  # set small value if you want to reduce load
OUT_CSV = "rtt-matrix/gt/ground_truth_rtt_matrix.csv"

# Option A: define nodes + IPs explicitly
# (edit to match your IP plan)
NODES = [
    ("serf1",  "10.0.1.10"),
    ("serf2",  "10.0.1.11"),
    ("serf3",  "10.0.1.12"),
    ("serf4",  "10.0.1.13"),
    ("serf5",  "10.0.1.14"),
    ("serf6",  "10.0.1.15"),
    ("serf7",  "10.0.1.16"),
    ("serf8",  "10.0.1.17"),
    ("serf9",  "10.0.1.18"),
    ("serf10", "10.0.1.19"),

    ("serf11", "10.0.2.10"),
    ("serf12", "10.0.2.11"),
    ("serf13", "10.0.2.12"),
    ("serf14", "10.0.2.13"),
    ("serf15", "10.0.2.14"),
    ("serf16", "10.0.2.15"),
    ("serf17", "10.0.2.16"),
    ("serf18", "10.0.2.17"),
    ("serf19", "10.0.2.18"),
    ("serf20", "10.0.2.19"),

    ("serf21", "10.0.3.10"),
    ("serf22", "10.0.3.11"),
    ("serf23", "10.0.3.12"),
    ("serf24", "10.0.3.13"),
    ("serf25", "10.0.3.14"),
    ("serf26", "10.0.3.15"),
    ("serf27", "10.0.3.16"),
    ("serf28", "10.0.3.17"),
    ("serf29", "10.0.3.18"),
    ("serf30", "10.0.3.19"),
]

# ----------------------------
# Helpers
# ----------------------------

PING_RE = re.compile(r"=\s*([\d\.]+)/([\d\.]+)/([\d\.]+)/([\d\.]+)\s*ms")  # min/avg/max/mdev

def docker_exec(container: str, cmd: list[str]) -> tuple[int, str, str]:
    full = ["sudo", "docker", "exec", container] + cmd
    p = subprocess.run(full, capture_output=True, text=True)
    return p.returncode, p.stdout, p.stderr

def ping_avg_ms(from_container: str, dst_ip: str) -> float | None:
    # Busybox ping and iputils ping formats differ a bit, but most include "min/avg/max"
    cmd = [
        "ping",
        "-c", str(PING_COUNT),
        "-W", str(PING_TIMEOUT_S),
        dst_ip
    ]
    rc, out, err = docker_exec(from_container, cmd)
    text = out + "\n" + err

    # Look for the min/avg/max line
    m = PING_RE.search(text)
    if not m:
        return None

    # group(2) is avg
    return float(m.group(2))

def ping_min_ms(from_container: str, dst_ip: str) -> float | None:
    cmd = ["ping", "-c", str(PING_COUNT), "-W", str(PING_TIMEOUT_S), dst_ip]
    rc, out, err = docker_exec(from_container, cmd)
    text = out + "\n" + err

    m = PING_RE.search(text)
    if not m:
        return None

    return float(m.group(1))   # MIN RTT

def main() -> None:
    names = [n for (n, _ip) in NODES]
    ip_map = {n: ip for (n, ip) in NODES}

    # Initialize matrix with None
    matrix: dict[str, dict[str, float | None]] = {src: {dst: None for dst in names} for src in names}

    for i, src in enumerate(names, start=1):
        src_container = f"clab-{LAB}-{src}"
        print(f"[{i}/{len(names)}] From {src_container}")

        for j, dst in enumerate(names, start=1):
            if src == dst:
                matrix[src][dst] = 0.0
                continue

            dst_ip = ip_map[dst]
            #avg = ping_avg_ms(src_container, dst_ip)
            rtt = ping_min_ms(src_container, dst_ip)

            #matrix[src][dst] = avg
            matrix[src][dst] = rtt
            status = f"{rtt:.2f} ms" if rtt is not None else "NA"
            print(f"  -> {dst:6s} {dst_ip:12s} = {status}")

            if SLEEP_BETWEEN_PINGS_S > 0:
                time.sleep(SLEEP_BETWEEN_PINGS_S)

    # Write matrix CSV (header row + each row)
    out_path = Path(OUT_CSV)
    with out_path.open("w", newline="") as f:
        w = csv.writer(f)
        w.writerow([""] + names)  # top-left blank cell, then column headers
        for src in names:
            row = [src]
            for dst in names:
                v = matrix[src][dst]
                row.append("" if v is None else f"{v:.3f}")
            w.writerow(row)

    print(f"\n[OK] Wrote RTT matrix to: {out_path.resolve()}")

if __name__ == "__main__":
    main()
