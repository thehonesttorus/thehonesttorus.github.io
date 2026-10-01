import sys; sys.path.insert(0, "../../bench"); sys.path.insert(0, ".")
import numpy as np, eval_q
from tct0 import tct0
eval_q.evaluate(lambda W: tct0([w for w in W]), sys.argv[1].split(","), units_1024=19)
