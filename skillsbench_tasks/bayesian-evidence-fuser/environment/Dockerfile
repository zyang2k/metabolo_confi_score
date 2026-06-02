# Minimal Ubuntu-based environment for the bayesian-evidence-fuser task.
FROM ubuntu:22.04

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update \
    && apt-get install -y --no-install-recommends \
        python3 \
        python3-pip \
    && pip3 install --no-cache-dir pandas==2.2.2 scipy==1.13.1 pytest==8.4.1 pytest-json-ctrf==0.3.5 \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /root

# The input feature table the agent must read.
COPY data /root/environment/data

# NOTE: Skills (environment/skills/bayesian-evidence-fuser/SKILL.md) are
# deployed at runtime by the benchmark harness, NOT baked into the image.
# Baking them in would leak the skill into the no-skills experimental arm and
# (because bench excludes skills/ from the build context) break the build.
