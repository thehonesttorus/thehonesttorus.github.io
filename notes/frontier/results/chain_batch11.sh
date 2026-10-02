#!/bin/bash
cd /home/user/thehonesttorus.github.io/notes/frontier
while [ $(ls results/sw_fbs150_lam065.json results/sw_fbs200_lam065.json 2>/dev/null | wc -l) -lt 2 ]; do sleep 15; done
sleep 5
python3 sweep.py cfg_batch11.json 0,1,2 2 2 > results/sweep_batch11.log 2>&1
