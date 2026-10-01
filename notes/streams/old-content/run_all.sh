#!/bin/bash
# full width-128 batch on the N = 6e5 atlases (A1, A2: same MLP seed 770100, independent samples; B1: seed 770101)
export PYTHONDONTWRITEBYTECODE=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
P=/root/whest/bin/python; A=/root/atl; R=results
cd "$(dirname "$0")"
RK=1,2,4,8,16,32,64
(
  $P tracker.py $A/A1/mlp_00000.npz > $R/tracker_A1.txt 2>&1
  $P tracker.py $A/A1/mlp_00000.npz --noad > $R/tracker_A1_noad.txt 2>&1
  $P tracker.py $A/A1/mlp_00000.npz --model --noad > $R/tracker_A1_model_noad.txt 2>&1
  $P tracker.py $A/B1/mlp_00000.npz --noad > $R/tracker_B1_noad.txt 2>&1
  $P tracker.py $A/A2/mlp_00000.npz --noad > $R/tracker_A2_noad.txt 2>&1
) &
(
  for w in 4 1; do $P carriers.py $A/A1/mlp_00000.npz --w $w --noad --which modes,modesC,modesdyn,ageregress,regress --ranks $RK > $R/carriers_A1_noad_w$w.txt 2>&1; done
  $P carriers.py $A/A1/mlp_00000.npz --w 4 --which modes,modesdyn,ageregress,regress --ranks $RK > $R/carriers_A1_ad_w4.txt 2>&1
  $P carriers.py $A/A1/mlp_00000.npz --w 4 --noad --model --which modes,modesdyn,ageregress --ranks $RK > $R/carriers_A1model_noad_w4.txt 2>&1
) &
(
  $P propproj.py $A/A1/mlp_00000.npz --noad --w 4 --ks 4,8,16,32,64 > $R/propproj_A1_noad_w4.txt 2>&1
  $P propproj.py $A/A1/mlp_00000.npz --model --noad --w 4 --ks 4,8,16,32,64 > $R/propproj_A1model_noad_w4.txt 2>&1
  $P propproj.py $A/B1/mlp_00000.npz --noad --w 4 --ks 4,8,16,32,64 > $R/propproj_B1_noad_w4.txt 2>&1
  $P propproj.py $A/A1/mlp_00000.npz --w 4 --ks 4,8,16,32,64 > $R/propproj_A1_ad_w4.txt 2>&1
  $P carriers.py $A/B1/mlp_00000.npz --w 4 --noad --which modes,modesdyn,ageregress --ranks $RK > $R/carriers_B1_noad_w4.txt 2>&1
) &
(
  $P carriers.py $A/A1/mlp_00000.npz --w 4 --noad --which cp --cpranks 32,64,128,256 --cplayers 6,10,14 --iters 40 > $R/cp_A1_noad_w4.txt 2>&1
  $P carriers.py $A/A1/mlp_00000.npz --w 4 --noad --which cpdyn --cpranks 64,128 --iters 25 > $R/cpdyn_A1_noad_w4.txt 2>&1
) &
wait
echo ALLDONE
