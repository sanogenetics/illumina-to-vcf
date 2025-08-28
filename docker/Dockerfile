# ---- build the wheel ----
FROM python:3.11-slim AS buildenv
WORKDIR /app
RUN pip install --no-cache-dir hatch
COPY pyproject.toml README.md ./
COPY illumina2vcf ./illumina2vcf
RUN hatch build -t wheel

# ---- runtime (Sano base) ----
FROM sanogenetics/htslib-bcftools-samtools:latest AS runtime
WORKDIR /app

# Python for runtime (keep AWS CLI as-is)
RUN apt-get update \
  && apt-get install -y --no-install-recommends python3 python3-pip ca-certificates \
  && rm -rf /var/lib/apt/lists/*

# Install your wheel
COPY --from=buildenv /app/dist/*.whl /tmp/
RUN pip install --no-cache-dir /tmp/*.whl --break-system-packages && rm -f /tmp/*.whl

ENV AWS_CLI_ARGS=
COPY docker/run.sh /app/run.sh
CMD ["sh", "/app/run.sh"]