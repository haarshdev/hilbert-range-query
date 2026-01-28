#!/bin/bash

# The node that will run the 'serf join' command
joining_node="clab-nebula-serf1"

# List of containers (Ubuntu nodes from 2 to 162)
containers=()
for i in {2..30}; do
  containers+=(clab-nebula-serf$i)
done

# Gather IP addresses of all nodes and map them to their container names
declare -A ip_to_node  # Associative array to store IP -> Node name mapping
ip_addresses=()

for container in "${containers[@]}"; do
  # Attempt to get the IP address of eth1
  ip_address=$(docker exec "$container" ip -4 addr show eth1 | grep -oP '(?<=inet\s)\d+(\.\d+){3}')
  
  # Check if the IP was found
  if [[ -z "$ip_address" ]]; then
    echo "Warning: Could not find IP address for eth1 on container $container"
  else
    ip_addresses+=("$ip_address")
    ip_to_node["$ip_address"]="$container"
  fi
done

# Check if any IP addresses were collected
if [[ ${#ip_addresses[@]} -eq 0 ]]; then
  echo "Error: No IP addresses were found on eth1 for the specified containers."
  exit 1
fi

# Join each node separately with a sleep timer in between
for ip in "${ip_addresses[@]}"; do
  target_node="${ip_to_node[$ip]}"  # Get the node name for the IP
  join_command="/opt/serfapp/serf join $ip"

  # Print formatted message
  echo "Joining target node $target_node ($ip) to $joining_node with the following command:"
  echo "$join_command"

  # Execute the join command on the first node
  docker exec "$joining_node" bash -c "$join_command"
  
  # Sleep for a specified duration (e.g., 5 seconds)
  #sleep 1
done

echo "Cluster joined successfully."
