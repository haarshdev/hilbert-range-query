FROM ubuntu:24.04

ENV DEBIAN_FRONTEND=noninteractive
WORKDIR /opt/serfapp

# Update the package list and install essential packages
RUN apt update && apt install -y \
    iproute2 \
    iputils-ping \
    net-tools \
    curl \
    && apt clean && rm -rf /var/lib/apt/lists/*
    
# Set the default command to bash
CMD ["bash"]
