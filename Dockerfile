FROM ghcr.io/prefix-dev/pixi:latest

WORKDIR /app

# Copy configuration files
COPY pyproject.toml pixi.lock* ./

# Copy source code
COPY src/ ./src/
COPY README.md LICENSE ./

# Install dependencies (including STAR from bioconda via pixi)
# Install dependencies and the package itself
RUN pixi install && \
    pixi run python -m pip install .

# Add the pixi environment to PATH so Nextflow can find the executable
ENV PATH="/app/.pixi/envs/default/bin:$PATH"
ENV PYTHONPATH="/app/src:$PYTHONPATH"

# Default command
CMD ["getRPF", "--help"]
