FROM ghcr.io/prefix-dev/pixi:latest

WORKDIR /app

# Copy configuration files
COPY pyproject.toml pixi.lock* ./

# Copy source code
COPY src/ ./src/
COPY README.md LICENSE ./

# Install dependencies (including STAR from bioconda via pixi)
RUN pixi install

# Set entrypoint to run getRPF
# Note: we run via pixi run to ensure environment is activated
ENTRYPOINT ["pixi", "run", "getRPF"]
