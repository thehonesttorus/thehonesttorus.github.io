#!/usr/bin/env bash
# Build the runner images in Azure Container Registry (cloud build; no local docker needed):
#   $IMAGE      from runner/Dockerfile      (CPU grid tasks, dataset staging)
#   $GPU_IMAGE  from runner/Dockerfile.gpu  (bakes and moment atlases on the GPU pool)
# ACR Tasks (az acr build) are refused on free-trial / free-credit subscriptions with
# TasksOperationsNotAllowed; the script then falls back to a local docker build + push.
# SKIP_GPU_IMAGE=1 builds only the CPU image.
source "$(dirname "$0")/env.sh"
cd "$(dirname "$0")/runner"
build() {  # $1 Dockerfile, $2 repository:tag
  if ! az acr build -r "$ACR" -t "$2" -f "$1" --platform linux .; then
    command -v docker >/dev/null || { echo "az acr build failed and docker is not installed" >&2; return 1; }
    echo "az acr build failed; falling back to local docker build + push" >&2
    az acr login -n "$ACR"
    docker build -f "$1" -t "$ACR.azurecr.io/$2" .
    docker push "$ACR.azurecr.io/$2"
  fi
  echo "image: $ACR.azurecr.io/$2"
}
build Dockerfile "$IMAGE"
if [ "${SKIP_GPU_IMAGE:-0}" != "1" ]; then build Dockerfile.gpu "$GPU_IMAGE"; fi
