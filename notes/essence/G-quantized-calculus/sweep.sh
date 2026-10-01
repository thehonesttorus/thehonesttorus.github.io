export OMP_NUM_THREADS=4
python run_g.py w256_d16 0,1,2,3,4,5,6,7 0,0.5,1,1.5,2,3 t3_w256 diag > results/t3_w256.log 2>&1
python run_g.py w512_d16 0,1,2,3 0,0.5,1,1.5,2,3 t3_w512 diag > results/t3_w512.log 2>&1
python run_g.py w1024_d16 0,1 0,1,1.5,2,3 t3_w1024a diag > results/t3_w1024a.log 2>&1
python run_g.py w1024_d16 2,3,4,5 0,1,2 t3_w1024b > results/t3_w1024b.log 2>&1
