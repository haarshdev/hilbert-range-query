#!/usr/bin/env bash
set -euo pipefail

LAB="nebula"          # topology name
QUEUE_LENGTH=5000     # like before (packets)

# Helper: set netem delay + queue sizes on a container interface (idempotent)
set_delay() {
  local node="$1" iface="$2" delay_ms="$3"
  local c="clab-${LAB}-${node}"

  echo "[+] ${c}:${iface} <= ${delay_ms}ms (txqueuelen=${QUEUE_LENGTH}, netem_limit=${QUEUE_LENGTH})"
  sudo docker exec "$c" sh -lc "
    ip link set ${iface} up || true
    ip link set ${iface} txqueuelen ${QUEUE_LENGTH} || true
    tc qdisc replace dev ${iface} root netem delay ${delay_ms}ms limit ${QUEUE_LENGTH}
  "
}

# -------------------------
# Intra-cluster "arms"
# -------------------------
set_delay "serf1"  "eth1" 0
set_delay "serf2"  "eth1" 1
set_delay "serf3"  "eth1" 2
set_delay "serf4"  "eth1" 3
set_delay "serf5"  "eth1" 4
set_delay "serf6"  "eth1" 5
set_delay "serf7"  "eth1" 6
set_delay "serf8"  "eth1" 7
set_delay "serf9"  "eth1" 8
set_delay "serf10" "eth1" 9

set_delay "serf11" "eth1" 0
set_delay "serf12" "eth1" 1
set_delay "serf13" "eth1" 2
set_delay "serf14" "eth1" 3
set_delay "serf15" "eth1" 4
set_delay "serf16" "eth1" 5
set_delay "serf17" "eth1" 6
set_delay "serf18" "eth1" 7
set_delay "serf19" "eth1" 8
set_delay "serf20" "eth1" 9

set_delay "serf21" "eth1" 0
set_delay "serf22" "eth1" 1
set_delay "serf23" "eth1" 2
set_delay "serf24" "eth1" 3
set_delay "serf25" "eth1" 4
set_delay "serf26" "eth1" 5
set_delay "serf27" "eth1" 6
set_delay "serf28" "eth1" 7
set_delay "serf29" "eth1" 8
set_delay "serf30" "eth1" 9

# -------------------------
# Inter-cluster uplink offsets on the router
# -------------------------
set_delay "router" "eth1" 5
set_delay "router" "eth2" 7
set_delay "router" "eth3" 9

echo "[OK] Delays applied."
