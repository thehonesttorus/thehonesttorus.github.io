#!/bin/bash
cd /home/user/thehonesttorus.github.io/notes/fresh-slate/bench
O=../../frontier/F2/results
for a in 8 11 13; do
  F2_RAISE=1 F2_DROP=1 F2_AGE3=$a /root/whest/bin/python run_p_inproc.py --estimator ../../frontier/F2/estimator_f2.py --set w1024_d16 --mlps 0,1 --json $O/drop_a$a.json > $O/drop_a$a.log 2>&1
done
echo ALLDONE > $O/sweep2.done
