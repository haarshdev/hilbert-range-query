#!/bin/bash

for i in $(seq 1 162); do
    CONTAINER="clab-nebula-serf$i"
    echo "Cleaning on $CONTAINER..."
    docker exec "$CONTAINER" sh -c "rm -f /opt/serfapp/endpoint /opt/serfapp/endpoint.log /opt/serfapp/rqserf /opt/serfapp/range-query-serf"
done

echo "Cleanup complete on all containers."
