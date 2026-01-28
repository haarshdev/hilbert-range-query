#!/bin/bash

# List of containers (Ubuntu nodes for 162 nodes)
containers=()
for i in {1..30}; do
  containers+=(clab-nebula-serf$i)
done

# Start Serf agent on each node
start_serf_agents() {
  for container in "${containers[@]}"; do
    echo "Starting Serf agent on $container..."
    
    # Ensure the serf binary is executable
    docker exec "$container" chmod +x /opt/serfapp/serf
    
    # Start the serf agent with the specified config file
    # docker exec -d "$container" /opt/serfapp/serf_2D agent -profile wan -config-file=/opt/serfapp/node.json -log-level debug
    # docker exec -d "$container" bash -c '/opt/serfapp/serf agent -log-level=debug -profile=lan -config-file=/opt/serfapp/node.json > /opt/serfapp/serf.log 2>&1'
    docker exec -d "$container" bash -c '/opt/serfapp/serf agent -log-level=debug -profile=lan -config-file=/opt/serfapp/node.json'
    
    echo "Serf agent started on $container."
    echo ""
    #sleep 1
  done
}

# Main script execution
start_serf_agents

