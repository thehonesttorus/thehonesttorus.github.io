"""F1 wrapper over the public V29 chain (504aldo, MIT): thin-leg ranks from the environment.
F1_RFB (default 16) = D21 feedback rank, F1_RRES (default 16) = S21 residual leg rank."""
import importlib.util as _iu
import os as _os
_p = _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), '..', 'estimator_v29r3.py')
_s = _iu.spec_from_file_location('estimator_v29r3_base', _p)
_m = _iu.module_from_spec(_s)
_s.loader.exec_module(_m)


class Estimator(_m.Estimator):
    R_FB = int(_os.environ.get('F1_RFB', '16'))
    R_RES = int(_os.environ.get('F1_RRES', '16'))
