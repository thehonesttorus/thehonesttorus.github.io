#!/bin/bash
cd /home/user/thehonesttorus.github.io/notes/frontier/d21
V17_DEBUG=1 H_FEED_SCALE=1.3 /root/whest/bin/python rec_d21.py est_rec_F1.py 0 rec2_F1_feed13_m0.npz > rec31.log 2>&1
V17_DEBUG=1 V17_LAM_SCALE=0.95 /root/whest/bin/python rec_d21.py est_rec_F1.py 0 rec2_F1_lam95_m0.npz > rec32.log 2>&1
echo REC_DONE >> rec32.log
