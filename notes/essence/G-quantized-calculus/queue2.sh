while pgrep -f sweep.sh > /dev/null; do sleep 20; done
export OMP_NUM_THREADS=4
P=pow:3.46:1.5,pow:6:2,lin:2:0.11:0.015625,lin:2:0.07:0.015625
python run_prof.py w256_d16 0,1,2,3 $P prof_w256 > results/prof_w256.log 2>&1
python run_prof.py w512_d16 0,1,2,3 $P prof_w512 > results/prof_w512.log 2>&1
python run_prof.py w1024_d16 0,1 $P prof_w1024 > results/prof_w1024.log 2>&1
