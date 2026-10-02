#!/bin/bash
cd /home/user/thehonesttorus.github.io/notes/frontier
until grep -q SWEEP_DONE results/sweep_batch2.log; do sleep 15; done
python3 sweep.py cfg_batch3.json 0,1,2 1 4 > results/sweep_batch3.log 2>&1
