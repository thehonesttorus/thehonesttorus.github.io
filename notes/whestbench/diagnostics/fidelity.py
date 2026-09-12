# How faithful is the numpy twin to the flopscope-metered harness?  Same 16 MLPs (the local shard).
import numpy as np, cov_rho as CR
v=np.array([CR.predict(k,S3=0.22) for k in range(16)])
print("twin mean over MLPs 0-15: %.4e" % v.mean())
print("harness run20           : 8.0987e-07")
print("relative difference     : %.3f%%" % (100*(v.mean()-8.0987e-07)/8.0987e-07))
