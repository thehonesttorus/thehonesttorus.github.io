#!/bin/bash
cd /home/user/thehonesttorus.github.io/notes/fresh-slate/bench
O=../../frontier/F2/results
F2_RAISE=1 F2_POSTT=1 V21_R_OLD=320 /root/whest/bin/python run_p_inproc.py --estimator ../../frontier/F2/estimator_f2.py --set w1024_d16 --mlps 0,1,2 --json $O/postt_r320.json > $O/postt_r320.log 2>&1
F2_RAISE=1 V21_R_OLD=320 /root/whest/bin/python run_p_inproc.py --estimator ../../frontier/F2/estimator_f2.py --set w1024_d16 --mlps 0,1,2 --json $O/base_r320.json > $O/base_r320.log 2>&1
echo ALLDONE > $O/sweep3.done
