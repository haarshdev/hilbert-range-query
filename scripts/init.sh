#!/usr/bin/env bash
set -Eeuo pipefail

HOSTNAME="$(hostname)"
IFACE_DATA="eth1"
CFG="/ip-mapping.txt"    # your mapping file inside the container

log() { echo "[INIT] $*"; }

log "Waiting for $IFACE_DATA to appear"
for _ in {1..40}; do
  ip link show "$IFACE_DATA" &>/dev/null && break
  sleep 0.5
done

NODE_IP_DATA=""
GATEWAY=""

# File format:
# net_name subnet hostname cidr gw
while read -r net_name subnet name cidr gw || [[ -n "${name:-}" ]]; do
  # skip empty lines and comments
  [[ -z "${name:-}" || "$net_name" =~ ^# ]] && continue

  if [[ "$HOSTNAME" == "$name" ]]; then
    log "Configuring $IFACE_DATA for $HOSTNAME with $cidr gw $gw"

    ip link set "$IFACE_DATA" up
    ip addr flush dev "$IFACE_DATA" || true
    ip addr add "$cidr" dev "$IFACE_DATA"

    # replace default route to use eth1 gateway
    ip route del default || true
    ip route add default via "$gw" dev "$IFACE_DATA" || true

    NODE_IP_DATA="${cidr%%/*}"
    GATEWAY="$gw"

    log "Configured $IFACE_DATA: IP=$cidr, Gateway=$gw"
    break
  fi
done < "$CFG"

if [[ -z "$NODE_IP_DATA" ]]; then
  log "ERROR: could not find $HOSTNAME in $CFG"
  exit 1
fi

log "Initialization complete for $HOSTNAME"
