#!/bin/bash
cd /home/user/thehonesttorus.github.io/notes/frontier
while pgrep -f "sweep.py cfg_batch1" > /dev/null; do sleep 15; done
python3 sweep.py cfg_batch2.json 0,1,2 2 2 > results/sweep_batch2.log 2>&1
