#!/bin/bash
cd /home/user/thehonesttorus.github.io/notes/frontier
until grep -q SWEEP_DONE results/sweep_batch3.log 2>/dev/null; do sleep 15; done
python3 sweep.py cfg_batch4.json 0,1,2 2 2 > results/sweep_batch4.log 2>&1
