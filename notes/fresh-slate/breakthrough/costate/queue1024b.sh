#!/bin/bash
cd /home/user/thehonesttorus.github.io/notes/fresh-slate/breakthrough/costate
while kill -0 1775 2>/dev/null; do sleep 10; done
export PRED_DIR=/tmp/claude-0/-home-user-thehonesttorus-github-io/73934192-bb11-53f5-adb1-9561b059184f/scratchpad/preds OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4
RES_SUFFIX=_nc /root/whest/bin/python run.py --set w1024_d16 --variants A1gl_nc,A2gl_nc,A3gl_nc,A3gsl_nc > logs_w1024_nc.txt 2>&1
