#!/bin/bash
cd /home/user/thehonesttorus.github.io/notes/frontier
until grep -q SWEEP_DONE results/sweep_batch5.log 2>/dev/null; do sleep 15; done
python3 sweep.py cfg_batch6.json 0,1,2 2 2 > results/sweep_batch6.log 2>&1
