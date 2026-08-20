FROM python:3.13-slim

# System dependencies required to build scientific Python packages
RUN apt-get update \
    && apt-get install -y --no-install-recommends \
        build-essential \
    && rm -rf /var/lib/apt/lists/*

# Install uv
COPY --from=ghcr.io/astral-sh/uv:latest /uv /uvx /bin/

WORKDIR /app

# Install Python dependencies into /app/.venv.
# Keeping this before COPY . gives us good Docker layer caching.
COPY pyproject.toml uv.lock ./

RUN uv sync --frozen --no-install-project

# Copy application source
COPY . .

# Install the project
RUN uv sync --frozen

EXPOSE 8000

CMD ["uv", "run", "uvicorn", "main:app", "--host", "0.0.0.0", "--port", "8000"]
