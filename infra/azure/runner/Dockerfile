# Runner image for the WhestBench experiment grid.  Same harness pins as the grader
# (whestbench 0.16.x / flopscope 0.12.x); evaluation uses the official scoring path.
FROM python:3.11-slim
ENV PYTHONUNBUFFERED=1 PIP_NO_CACHE_DIR=1 HF_HOME=/data/hf \
    OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 MKL_NUM_THREADS=4
RUN apt-get update && apt-get install -y --no-install-recommends git curl && rm -rf /var/lib/apt/lists/*
RUN pip install "flopscope>=0.12.1,<0.13.0" "whestbench>=0.16.1,<0.17.0" \
    "datasets>=3" pyarrow huggingface_hub azure-storage-blob azure-identity numpy
WORKDIR /app
COPY stage_dataset.py run_shard.py /app/
# estimators are shipped per task (downloaded from the submissions container), not baked in
ENTRYPOINT []
