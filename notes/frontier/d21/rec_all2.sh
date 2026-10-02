#!/bin/bash
cd /home/user/thehonesttorus.github.io/notes/frontier/d21
V17_DEBUG=1 /root/whest/bin/python rec_d21.py est_rec_F1.py 0 rec2_F1final_m0.npz > rec21.log 2>&1
V17_DEBUG=1 /root/whest/bin/python rec_d21.py est_rec_v29.py 0 rec2_v29_m0.npz > rec23.log 2>&1
V17_DEBUG=1 H_FB_SCALE=1.0 /root/whest/bin/python rec_d21.py est_rec_F1.py 0 rec2_F1final_fb1_m0.npz > rec22.log 2>&1
echo REC_DONE >> rec22.log
