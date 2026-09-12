# Build a LOCAL copy of the full 100-MLP mini split so the flopscope-metered harness can score all of
# it.  Only shard 0 (16 MLPs) was downloaded; the weights of the other 84 are regenerable from
# mlp_seed -- verified bit-exact against the stored weights for MLPs 0, 1, 7, 15 -- and every other
# column is small enough to fetch by HTTP range read.  Written in batches so peak RAM stays low.
import numpy as np, pyarrow as pa, pyarrow.parquet as pq, fsspec, json, os, sys, shutil, time
SP=os.popen("ls -d /tmp/claude-0/-home-user-thehonesttorus-github-io/*/scratchpad").read().strip()
BASE=("https://huggingface.co/datasets/aicrowd/arc-whestbench-public-2026/resolve/v2-phase2/data/"
      "mini-{:05d}-of-00007.parquet")
SMALL=['mlp_id','mlp_name','mlp_seed','all_layer_means','avg_variance','sampling_budget_breakdown']
CACHE=SP+'/est/split_small.npz'
n=1024; L=16
def flat(c):
    a=c.combine_chunks() if isinstance(c,pa.ChunkedArray) else c
    while pa.types.is_list(a.type) or pa.types.is_fixed_size_list(a.type) or pa.types.is_large_list(a.type):
        a=a.flatten()
    return a.to_numpy(zero_copy_only=False)

def fetch_small():
    if os.path.exists(CACHE):
        d=np.load(CACHE,allow_pickle=True); return {k:d[k] for k in d.files}
    fs=fsspec.filesystem("https")
    ids=[];names=[];seeds=[];alm=[];avgv=[];sbb=[]
    for i in range(7):
        for att in range(4):
            try:
                if i==0:
                    tb=pq.ParquetFile(SP+'/data/mini-00000-of-00007.parquet').read_row_group(0,columns=SMALL)
                else:
                    with fs.open(BASE.format(i), block_size=2**20) as f:
                        tb=pq.ParquetFile(f).read_row_group(0,columns=SMALL)
                break
            except Exception as e:
                print(f"  shard {i} attempt {att}: {e}",flush=True); time.sleep(2**att)
        else: raise SystemExit(f"shard {i} failed")
        m=len(tb.column('mlp_seed'))
        ids+=tb.column('mlp_id').to_pylist(); names+=tb.column('mlp_name').to_pylist()
        seeds+=tb.column('mlp_seed').to_pylist()
        alm.append(flat(tb.column('all_layer_means')).astype(np.float32).reshape(m,L,n))
        avgv+=tb.column('avg_variance').to_pylist()
        sbb+=tb.column('sampling_budget_breakdown').to_pylist()
        print(f"  shard {i}: {m} MLPs",flush=True)
    d=dict(ids=np.array(ids),names=np.array(names),seeds=np.array(seeds),
           alm=np.concatenate(alm,0),avgv=np.array(avgv),sbb=np.array(sbb))
    np.savez_compressed(CACHE,**d); return d

def regen(seed):
    ss=np.random.SeedSequence(int(seed)).spawn(3); rng=np.random.default_rng(ss[0]); sc=float(np.sqrt(2.0/n))
    return np.stack([(rng.standard_normal((n,n))*sc).astype(np.float32) for _ in range(L)])

SCHEMA=pq.ParquetFile(SP+'/data/mini-00000-of-00007.parquet').schema_arrow

def nested(vals, dims):
    arr=pa.array(vals)                       # flat float32
    for d in reversed(dims[1:]):
        off=pa.array(np.arange(0,len(arr)//d+1,dtype=np.int32)*d)
        arr=pa.ListArray.from_arrays(off,arr)
    off=pa.array(np.arange(0,len(arr)//dims[0]+1,dtype=np.int32)*dims[0])
    return pa.ListArray.from_arrays(off,arr)

def build(idx, out_dir, D):
    if os.path.exists(out_dir): shutil.rmtree(out_dir)
    os.makedirs(out_dir+'/data')
    W=np.concatenate([regen(D['seeds'][i]) for i in idx]).reshape(-1)
    cols={
      'mlp_id': pa.array([int(D['ids'][i]) for i in idx], pa.int32()),
      'mlp_name': pa.array([str(D['names'][i]) for i in idx], pa.string()),
      'mlp_seed': pa.array([int(D['seeds'][i]) for i in idx], pa.int64()),
      'weights': nested(pa.array(W, pa.float32()), [L,n,n]),
      'all_layer_means': nested(pa.array(D['alm'][idx].reshape(-1), pa.float32()), [L,n]),
      'final_means': pa.FixedSizeListArray.from_arrays(
            pa.array(D['alm'][idx][:,-1,:].reshape(-1), pa.float32()), n),
      'avg_variance': pa.array([float(D['avgv'][i]) for i in idx], pa.float64()),
      'sampling_budget_breakdown': pa.array([str(D['sbb'][i]) for i in idx], pa.string()),
    }
    tb=pa.Table.from_arrays([cols[f.name] for f in SCHEMA], schema=SCHEMA)
    pq.write_table(tb, out_dir+f'/data/mini-00000-of-00001.parquet')
    md=json.load(open(SP+'/localds/metadata.json'))
    md['splits']={'mini':{**md['splits']['mini'],'n_mlps':len(idx)}}
    md['default_split']='mini'
    md.pop('prepared_splits',None); md.pop('partials_count',None)
    json.dump(md, open(out_dir+'/metadata.json','w'))
    del W, tb
    return len(idx)

if __name__=="__main__":
    D=fetch_small()
    print("total MLPs with small columns:", len(D['seeds']))
    which=sys.argv[1]
    if which=='verify':
        m=build(list(range(16)), SP+'/ds_verify', D)
        print(f"rebuilt {m} MLPs into ds_verify (should reproduce run20 exactly)")
    else:
        lo,hi=[int(x) for x in which.split('-')]
        m=build(list(range(lo,hi)), SP+f'/ds_{lo}_{hi}', D)
        print(f"built MLPs {lo}..{hi-1} ({m}) into ds_{lo}_{hi}")
