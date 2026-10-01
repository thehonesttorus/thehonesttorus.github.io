#!/usr/bin/env bash
# Build the runner image in Azure Container Registry (cloud build; no local docker needed).
source "$(dirname "$0")/env.sh"
cd "$(dirname "$0")/runner"
az acr build -r "$ACR" -t "$IMAGE" -f Dockerfile . 
echo "image: $ACR.azurecr.io/$IMAGE"
