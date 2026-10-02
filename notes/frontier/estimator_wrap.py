"""V29r3 (504aldo, MIT) with the thin-leg ranks R_FB / R_RES settable from the environment (H_R_FB, H_R_RES)."""
import os
import sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from estimator_v29r3 import Estimator as _Base  # noqa: E402


class Estimator(_Base):
    R_FB = int(os.environ.get("H_R_FB", "16"))
    R_RES = int(os.environ.get("H_R_RES", "16"))
