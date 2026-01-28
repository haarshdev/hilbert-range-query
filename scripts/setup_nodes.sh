#!/bin/bash

# List of containers (Ubuntu nodes for 100 nodes)
containers=()
for i in {1..30}; do
  containers+=(clab-nebula-serf$i)
done

# Paths and file names
json_file="node.json"
binary_file_loc="../serf_binaries/5D"  # Use the full path to the serf binary
destination_dir="/opt/serfapp"

: << 'EOF'
This is a multiline comment.
bind_addr: The address that Serf will bind to for communication with other Serf nodes. By default this is "0.0.0.0:7946". 
rpc_addr: The address that Serf will bind to for the agent's RPC server. By default this is "127.0.0.1:7373", allowing only loopback connections.
advertise: The advertise flag is used to change the address that we advertise to other nodes in the cluster. By default, the bind address is advertised.
EOF

# Function to set up Ubuntu nodes
setup_ubuntu_nodes() {
  for i in "${!containers[@]}"; do
    container="${containers[$i]}"
    
    # Check if container is running
    if ! docker ps --format '{{.Names}}' | grep -q "$container"; then
      echo "Container $container is not running, skipping..."
      continue
    fi
    
    echo "Setting up $container..."

    # Get the IP address of the eth1 interface directly from within the container
    ip_address=$(docker exec $container ip -4 addr show eth1 | grep -oP '(?<=inet\s)\d+(\.\d+){3}')
    if [ -z "$ip_address" ]; then
      echo "Failed to retrieve IP address for $container"
      continue
    fi
    
    echo "IP address for $container (eth1): $ip_address"
    
    # Generate the JSON configuration file dynamically
    json_content=$(cat <<EOF
{
  "node_name": "$container",
  "bind": "0.0.0.0:7946",
  "rpc_addr": "0.0.0.0:7373",
  "advertise": "$ip_address:7946",
  "log_level": "debug"
}
EOF
)
    
    # Create a temporary JSON file on the host
    temp_json_file=$(mktemp)
    echo "$json_content" > "$temp_json_file"

    # Create the destination directory inside the container
    docker exec "$container" mkdir -p "$destination_dir"

    # Copy the generated JSON file and serf binary into the /opt/serfapp/ directory
    docker cp "$temp_json_file" "$container":"$destination_dir"/node.json || { echo "Failed to copy node.json to $container"; exit 1; }
    docker cp "$binary_file_loc/serf" "$container":"$destination_dir"/ || { echo "Failed to copy serf binary to $container"; exit 1; }

    # Make the binary executable
    docker exec "$container" chmod +x "$destination_dir"/"serf" || { echo "Failed to make serf executable on $container"; exit 1; }

    # Remove the temporary JSON file
    rm "$temp_json_file"

    echo "$container setup complete."
  done
}

# Main script execution
setup_ubuntu_nodes

