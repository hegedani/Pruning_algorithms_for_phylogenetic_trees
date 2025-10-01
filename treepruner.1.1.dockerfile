FROM python:3.12-slim AS base

# Avoid interactive prompts
ARG DEBIAN_FRONTEND=noninteractive

RUN pip install ete3==3.1.3 numpy==2.3.3 matplotlib==3.10.6

COPY bin/treepruner.py /usr/local/bin/treepruner.py
RUN cat << 'EOF' > /usr/local/bin/treepruner
#!/bin/sh
python3 /usr/local/bin/treepruner.py "$@"
EOF
RUN chmod +x /usr/local/bin/treepruner


# Default entrypoint
ENTRYPOINT ["/bin/bash"]
