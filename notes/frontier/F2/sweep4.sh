#!/bin/bash
cd /home/user/thehonesttorus.github.io/notes/fresh-slate/bench
O=../../frontier/F2/results
/root/whest/bin/python run_p_inproc.py --estimator ../../frontier/estimator_v29r3.py --set w1024_d16 --mlps 3,4,5 --json $O/base_345.json > $O/base_345.log 2>&1
F2_RAISE=1 F2_POSTT=1 /root/whest/bin/python run_p_inproc.py --estimator ../../frontier/F2/estimator_f2.py --set w1024_d16 --mlps 3,4,5 --json $O/postt_345.json > $O/postt_345.log 2>&1
echo ALLDONE > $O/sweep4.done
