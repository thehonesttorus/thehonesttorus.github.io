#!/bin/bash
# F2 screening: exact check, then (AGE3, K3) grid on MLPs 0,1
cd /home/user/thehonesttorus.github.io/notes/fresh-slate/bench
O=../../frontier/F2/results
mkdir -p $O
F2_AGE3=13 F2_K3=288 F2_EXACT=1 /root/whest/bin/python run_p_inproc.py --estimator ../../frontier/F2/estimator_f2.py --set w1024_d16 --mlps 0 --json $O/exact_a13.json > $O/exact_a13.log 2>&1
for cfg in "8 96" "8 128" "10 64" "10 96" "11 64" "8 64"; do
  set -- $cfg
  F2_AGE3=$1 F2_K3=$2 /root/whest/bin/python run_p_inproc.py --estimator ../../frontier/F2/estimator_f2.py --set w1024_d16 --mlps 0,1 --json $O/a$1_k$2.json > $O/a$1_k$2.log 2>&1
done
echo ALLDONE > $O/sweep1.done
