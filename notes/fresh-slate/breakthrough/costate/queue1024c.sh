#!/bin/bash
cd /home/user/thehonesttorus.github.io/notes/fresh-slate/breakthrough/costate
while kill -0 1820 2>/dev/null; do sleep 10; done
export PRED_DIR=/tmp/claude-0/-home-user-thehonesttorus-github-io/73934192-bb11-53f5-adb1-9561b059184f/scratchpad/preds OMP_NUM_THREADS=3 OPENBLAS_NUM_THREADS=3
RES_SUFFIX=_filt /root/whest/bin/python run.py --set w1024_d16 --mlps 0,1,2 --variants A1oldR32,A1hfitd,A1oldD,A3oldR8,A3gsm > logs_w1024_filt2.txt 2>&1
RES_SUFFIX=_filt /root/whest/bin/python run.py --set w1024_d16 --mlps 1,2 --variants A1oldR1,A1oldR8,A1gsm >> logs_w1024_filt2.txt 2>&1
