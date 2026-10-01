"""One GPU task: bake slice K of M of a dataset of fresh He-init MLPs with the official recipe,
then upload the baked directory to blob.  Seeds come from a shared seeds.json so that slices
are disjoint and the merge (`whest dataset merge`) reproduces one dataset.

    python bake_shard.py --name fresh-A --seeds <url|path to seeds.json> --n-mlps 1000 \
        --n-samples 100000000 --slice 3/16 --dest 'https://<acct>.blob.core.windows.net/dataset?<sas>' \
        [--split dev] [--chunk-size 131072]
"""
import argparse, json, os, shutil, subprocess, sys, tempfile, urllib.parse


def fetch(src, dest):
    if src.startswith("http"):
        from azure.storage.blob import BlobClient
        with open(dest, "wb") as f:
            BlobClient.from_blob_url(src).download_blob().readinto(f)
    else:
        shutil.copyfile(src, dest)


def upload_dir(local, container_url, prefix):
    from azure.storage.blob import ContainerClient
    cc = ContainerClient.from_container_url(container_url)
    for root, _, files in os.walk(local):
        for fn in files:
            p = os.path.join(root, fn)
            rel = os.path.relpath(p, local)
            with open(p, "rb") as f:
                cc.upload_blob(f"{prefix}/{rel}", f, overwrite=True, max_concurrency=8)
            print("uploaded", rel, flush=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--name", required=True)
    ap.add_argument("--seeds", required=True)
    ap.add_argument("--n-mlps", type=int, required=True)
    ap.add_argument("--n-samples", type=int, required=True)
    ap.add_argument("--slice", required=True, help="K/M")
    ap.add_argument("--dest", required=True, help="container URL with SAS")
    ap.add_argument("--split", default="dev")
    ap.add_argument("--width", type=int, default=1024)
    ap.add_argument("--depth", type=int, default=16)
    ap.add_argument("--chunk-size", type=int, default=131072)
    ap.add_argument("--device", default="cuda")
    args = ap.parse_args()
    work = tempfile.mkdtemp(prefix="bake-")
    seeds = os.path.join(work, "seeds.json")
    fetch(args.seeds, seeds)
    out = os.path.join(work, "out")
    k, m = args.slice.split("/")
    cmd = ["whest", "dataset", "bake", "--n-mlps", str(args.n_mlps), "--n-samples", str(args.n_samples),
           "--width", str(args.width), "--depth", str(args.depth), "--split", args.split, "--mlp-seeds", seeds,
           "--output", out, "--slice", args.slice, "--chunk-size", str(args.chunk_size)]
    if args.device != "cpu":
        cmd += ["--torch", "--device", args.device]
    print("running:", " ".join(cmd), flush=True)
    rc = subprocess.call(cmd)
    if rc != 0:
        sys.exit(rc)
    upload_dir(out, args.dest, f"bakes/{args.name}/slice-{int(k):03d}-of-{int(m):03d}")
    print("done", flush=True)


if __name__ == "__main__":
    main()
