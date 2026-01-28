#!/bin/bash

# Get a list of all OVS bridges
BRIDGES=$(sudo ovs-vsctl list-br)

if [[ -z "$BRIDGES" ]]; then
    echo "No OVS bridges found."
    exit 0
fi

# Delete each bridge
for br in $BRIDGES; do
    sudo ovs-vsctl del-br $br
    echo "Deleted OVS bridge: $br"
done

echo "All Open vSwitch bridges deleted successfully."
