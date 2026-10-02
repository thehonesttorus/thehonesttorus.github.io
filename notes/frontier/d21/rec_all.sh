#!/bin/bash
cd /home/user/thehonesttorus.github.io/notes/frontier/d21
V17_DEBUG=1 /root/whest/bin/python rec_d21.py ../F1/est_F1_final.py 0 rec_F1final_m0.npz > rec1.log 2>&1
V17_DEBUG=1 H_FB_SCALE=1.0 /root/whest/bin/python rec_d21.py ../F1/est_F1_final.py 0 rec_F1final_fb1_m0.npz > rec2.log 2>&1
V17_DEBUG=1 /root/whest/bin/python rec_d21.py ../estimator_v29r3.py 0 rec_v29_m0.npz > rec3.log 2>&1
echo REC_DONE >> rec3.log
