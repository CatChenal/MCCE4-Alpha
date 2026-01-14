FROM continuumio/miniconda3

# Install system tools
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    curl \
    git \
    iproute2 \
    net-tools \
    openssh-client \
    sudo \
    tree \
    && apt-get clean && rm -rf /var/lib/apt/lists/*

# Install Node via Conda
RUN conda install -c conda-forge nodejs

# Create user 'builder', add to sudo group, and enable passwordless sudo
RUN useradd -ms /bin/bash builder && \
    usermod -aG sudo builder && \
    echo '%sudo ALL=(ALL) NOPASSWD:ALL' >> /etc/sudoers

# 4. Switch to the new user
USER builder
WORKDIR /home/builder

# Optional: Ensure bash is used for RUN commands later (if you add any)
SHELL ["/bin/bash", "-l", "-c"]