#!/bin/bash
cd /home/user/thehonesttorus.github.io/notes/frontier
while [ $(ls results/sw_fbs150_lam080.json results/sw_fbs150_lam110.json 2>/dev/null | wc -l) -lt 2 ]; do sleep 15; done
sleep 5
python3 sweep.py cfg_batch9.json 0,1,2 2 2 > results/sweep_batch9.log 2>&1
