export PYTHONPATH=../num12
python3 eval_official.py 0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15
for k in 1 2 3 4 5 6; do
  curl -sSL -o mini$k.parquet "https://huggingface.co/datasets/aicrowd/arc-whestbench-public-2026/resolve/v2-phase2/data/mini-0000$k-of-00007.parquet" || curl -sSL -o mini$k.parquet "https://huggingface.co/datasets/aicrowd/arc-whestbench-public-2026/resolve/v2-phase2/data/mini-0000$k-of-00007.parquet"
  ids=$(python3 extract.py mini$k.parquet | tr ' ' ',' | sed 's/,$//')
  rm -f mini$k.parquet
  python3 eval_official.py $ids
done
echo DONE
