#!/bin/bash
cd /home/user/thehonesttorus.github.io/notes/fresh-slate/breakthrough/costate
while kill -0 1577 2>/dev/null; do sleep 10; done
export PRED_DIR=/tmp/claude-0/-home-user-thehonesttorus-github-io/73934192-bb11-53f5-adb1-9561b059184f/scratchpad/preds OMP_NUM_THREADS=3 OPENBLAS_NUM_THREADS=3
(OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 /root/whest/bin/python anatomy.py w1024_d16 0 1 > anatomy_w1024_mlp0_A1.txt 2>&1 &)
RES_SUFFIX=_filt /root/whest/bin/python run.py --set w1024_d16 --variants A1oldR1,A1oldR8,A1gsm,A1hfitd,A1oldR32,A1oldD > logs_w1024_filt.txt 2>&1
