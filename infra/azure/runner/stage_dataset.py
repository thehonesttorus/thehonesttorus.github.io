"""Download the public WhestBench dataset shards from Hugging Face and upload them to an
Azure blob container (run once, inside Azure).  Also uploads metadata.json and README.md.

    python stage_dataset.py --dataset aicrowd/arc-whestbench-public-2026 --revision v2-phase2 \
        --dest 'https://<account>.blob.core.windows.net/dataset?<sas>' [--splits mini full]
"""
import argparse, os, sys, tempfile


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dataset", required=True)
    ap.add_argument("--revision", required=True)
    ap.add_argument("--dest", required=True, help="container URL with SAS")
    ap.add_argument("--splits", nargs="+", default=["mini", "full"])
    ap.add_argument("--all-files", action="store_true", help="stage every file of the repo (community atlases)")
    ap.add_argument("--prefix", default=None, help="blob prefix (default: the revision)")
    ap.add_argument("--workdir", default=os.path.join(os.environ.get("AZ_BATCH_TASK_WORKING_DIR") or tempfile.gettempdir(), "stage"),
                    help="download scratch (default: the Batch task working directory, on the node resource disk)")
    args = ap.parse_args()
    os.makedirs(args.workdir, exist_ok=True)
    # keep the Hugging Face cache (incl. the hf_xet chunk cache) next to the scratch, not on the OS disk
    if os.environ.get("AZ_BATCH_TASK_WORKING_DIR"):
        os.environ["HF_HOME"] = os.path.join(args.workdir, "hf")
    from huggingface_hub import HfApi, hf_hub_download
    from azure.storage.blob import ContainerClient
    api = HfApi()
    files = api.list_repo_files(args.dataset, repo_type="dataset", revision=args.revision)
    wanted = files if args.all_files else [f for f in files if f in ("metadata.json", "README.md")
              or any(f.startswith(f"data/{s}-") and f.endswith(".parquet") for s in args.splits)]
    prefix = args.prefix or args.revision
    cc = ContainerClient.from_container_url(args.dest)
    existing = {b.name: b.size for b in cc.list_blobs()}
    for f in sorted(wanted):
        name = f"{prefix}/{f}"
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
