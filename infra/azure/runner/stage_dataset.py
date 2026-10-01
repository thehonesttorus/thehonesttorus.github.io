"""Download the public WhestBench dataset shards from Hugging Face and upload them to an
Azure blob container (run once, inside Azure).  Also uploads metadata.json and README.md.

    python stage_dataset.py --dataset aicrowd/arc-whestbench-public-2026 --revision v2-phase2 \
        --dest 'https://<account>.blob.core.windows.net/dataset?<sas>' [--splits mini full]
"""
import argparse, os, sys
from huggingface_hub import HfApi, hf_hub_download
from azure.storage.blob import ContainerClient


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dataset", required=True)
    ap.add_argument("--revision", required=True)
    ap.add_argument("--dest", required=True, help="container URL with SAS")
    ap.add_argument("--splits", nargs="+", default=["mini", "full"])
    ap.add_argument("--workdir", default="/data/stage")
    args = ap.parse_args()
    os.makedirs(args.workdir, exist_ok=True)
    api = HfApi()
    files = api.list_repo_files(args.dataset, repo_type="dataset", revision=args.revision)
    wanted = [f for f in files if f in ("metadata.json", "README.md")
              or any(f.startswith(f"data/{s}-") and f.endswith(".parquet") for s in args.splits)]
    cc = ContainerClient.from_container_url(args.dest)
    existing = {b.name: b.size for b in cc.list_blobs()}
    for f in sorted(wanted):
        name = f"{args.revision}/{f}"
        local = hf_hub_download(args.dataset, f, repo_type="dataset", revision=args.revision, local_dir=args.workdir)
        size = os.path.getsize(local)
        if existing.get(name) == size:
            print("skip (present)", name, flush=True)
        else:
            print("upload", name, size, flush=True)
            with open(local, "rb") as fh:
                cc.upload_blob(name, fh, overwrite=True, max_concurrency=8)
        os.remove(local)
    print("done", len(wanted), "files", flush=True)


if __name__ == "__main__":
    sys.exit(main())
