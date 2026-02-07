# Multi-stage Dockerfile for SCCmecExtractor with Bakta
# Build stage for compiling dependencies
FROM mambaorg/micromamba:1.5.10

USER root
RUN apt-get update && apt-get install -y procps git curl wget && rm -rf /var/lib/apt/lists/*

USER $MAMBA_USER

# Copy environment and install dependencies
COPY --chown=$MAMBA_USER:$MAMBA_USER environment.yml /tmp/environment.yml
RUN micromamba install -y -n base -f /tmp/environment.yml && micromamba clean --all --yes

# Set working directory
WORKDIR /app
COPY --chown=$MAMBA_USER:$MAMBA_USER . /app

# Install your package in editable mode using micromamba run
RUN micromamba run -n base pip install --no-cache-dir -e .

# Set environment variables for Bakta
ENV BAKTA_DB=/data/bakta_db

# Temporarily switch to root to create mount points
USER root
RUN mkdir -p /data /work && chown -R $MAMBA_USER:$MAMBA_USER /data /work

# Switch back to non-root user
USER $MAMBA_USER

# Set working directory for user data
WORKDIR /work

# Add entrypoint script
COPY --chown=$MAMBA_USER:$MAMBA_USER docker-entrypoint.sh /usr/local/bin/
RUN chmod +x /usr/local/bin/docker-entrypoint.sh

ENTRYPOINT ["/usr/local/bin/docker-entrypoint.sh"]
CMD ["/bin/bash"]

# Labels
LABEL maintainer="Alison MacFadyen <alison.macfadyen86@gmail.com>"
LABEL description="SCCmecExtractor: Pipeline for extracting SCCmec sequences from Staphylococcus genomes"
LABEL version="1.2.2"