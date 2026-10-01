for f in predict_gauss predict_win6 predict_hyb21_w16 predict_mem21; do
  /root/whest/bin/python run_q.py $f w1024_d16 > q1024_$f.log 2>&1
done
echo ALLDONE > q1024_done.flag
