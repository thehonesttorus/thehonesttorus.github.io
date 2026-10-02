#!/bin/bash
# wait for a running sweep to finish, then start the next: chain.sh PREVLOG CFG TAG
cd "$(dirname "$0")"
until grep -q SWEEP_DONE "$1"; do sleep 20; done
/root/whest/bin/python sweepF1.py "$2" ${4:-0,1,2} 2 2 > results/sweep$3.log 2>&1
