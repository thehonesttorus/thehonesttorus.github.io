#!/bin/bash
cd /home/user/thehonesttorus.github.io/notes/essence/H-rational-memory
python3 run_sliceage.py 0 all m0_all
for K in 1 2 4; do python3 run_sliceage.py 0,1,2 $K m012_K$K; done
echo DRIVER_DONE
