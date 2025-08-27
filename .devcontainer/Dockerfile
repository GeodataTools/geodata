FROM mcr.microsoft.com/devcontainers/base:bookworm

USER root
RUN apt-get update && export DEBIAN_FRONTEND=noninteractive \
    && apt-get -y install --no-install-recommends libgdal-dev \
    && rm -rf /var/lib/apt/lists/*

USER vscode
RUN curl -LsSf https://astral.sh/uv/install.sh | sh
