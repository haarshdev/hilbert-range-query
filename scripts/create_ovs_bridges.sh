#!/bin/bash

# Create 37 OVS bridges (switch1 to switch37)
for i in {a..c}; do
    sudo ovs-vsctl add-br switch_$i
    echo "Created OVS bridge: switch_$i"
done

echo "All Open vSwitch bridges created successfully."
