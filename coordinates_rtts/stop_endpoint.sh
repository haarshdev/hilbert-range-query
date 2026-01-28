#!/bin/bash

# Container prefix
PREFIX="clab-nebula-serf"

# Number of containers
NUM_CONTAINERS=162

# Loop over all containers
for i in $(seq 1 $NUM_CONTAINERS); do
    CONTAINER="${PREFIX}${i}"
    echo "Stopping api in container $CONTAINER..."
    
    # Use pkill to stop the binary by name inside the container
    docker exec "$CONTAINER" pkill -f "/opt/serfapp/endpoint" 2>/dev/null
    
    # Optional: check if it’s stopped
    docker exec "$CONTAINER" pgrep -f "/opt/serfapp/endpoint" >/dev/null
    if [ $? -eq 0 ]; then
        echo "Failed to stop api in $CONTAINER"
    else
        echo "Stopped successfully"
    fi
done

echo "Done."
