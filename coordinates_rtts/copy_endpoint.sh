#!/bin/bash

# Target path inside container
TARGET_PATH="/opt/serfapp/rtts-api"

# Loop through all containers 1 to 162
for i in $(seq 1 162); do
    CONTAINER="clab-nebula-serf$i"
    
    echo "Copying api to $CONTAINER..."
    docker cp ./rtts-api "$CONTAINER:$TARGET_PATH"
    
    echo "Starting api in background on $CONTAINER..."
    docker exec "$CONTAINER" sh -c "nohup $TARGET_PATH > /dev/null 2>&1 &"
    
    echo "$CONTAINER is running api in background."
done

echo "All 162 containers updated."
