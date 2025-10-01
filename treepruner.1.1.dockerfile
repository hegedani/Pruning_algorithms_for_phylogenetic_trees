FROM python:3.12-slim AS base

# Avoid interactive prompts
ARG DEBIAN_FRONTEND=noninteractive
RUN apt-get update && apt-get install -y --no-install-recommends python3-ete3 python3-numpy python3-matplotlib && rm -rf /var/lib/apt/lists/*

COPY bin/treepruner.py /usr/local/bin/treepruner
RUN chmod +x /usr/local/bin/treepruner

# Default entrypoint
ENTRYPOINT ["treepruner"]
