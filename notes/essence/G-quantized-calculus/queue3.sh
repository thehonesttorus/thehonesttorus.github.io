while pgrep -f queue2.sh > /dev/null; do sleep 20; done
export OMP_NUM_THREADS=4
P=pow:1.15:0.5,pow:1.2:0.75
python run_prof.py w256_d16 0,1,2,3 $P prof2_w256 > results/prof2_w256.log 2>&1
python run_prof.py w512_d16 0,1 $P prof2_w512 > results/prof2_w512.log 2>&1
