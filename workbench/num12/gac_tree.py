# GAC + fresh one-point tree corrections (Edgeworth mean shift from the previous layer's kinks), excluding the coherent
# pair fourth cumulant PP4, which is the gain injection GAC already carries.
import numpy as np, sys
from ledger import gstep
from gac import inject, EG
from edgeworth import cumulants_next, edgeworth_shift
TERMS = ("D3", "P3", "T3", "D4", "PD4")
def gac_tree(Ws, J=1, terms=TERMS, lmax=99):
    n = Ws[0].shape[0]; out = []; gam = 0.0; g1 = 1.0
    for l, W in enumerate(Ws):
        if l == 0: mu = np.zeros(n); S = W @ W.T
        else:
            nu = W @ mbar; Q = W @ Sbar @ W.T
            mu0 = nu/g1; S0 = Q - np.outer(mu0, mu0); gam = gam + inject(W, mu0, S0, prev)
            g1 = EG(gam); mu = nu/g1; S = Q - np.outer(mu, mu)
        M, C = gstep(mu, S); sig = np.sqrt(np.diag(S))
        if 0 < l <= lmax:
            pm, ps, pR = prev
            _, _, cu = cumulants_next(W, pm, ps, pR, J=J, terms=terms)
            k3 = sum(cu[t] for t in terms if t.endswith("3")); k4 = sum(cu[t] for t in terms if t.endswith("4"))
            M = M + edgeworth_shift(mu, sig, k3, k4)
        mbar = g1*M; Sbar = C + np.outer(M, M)
        prev = (mu, sig, S/np.outer(sig, sig))
        out.append(dict(m=mbar, gamma=gam))
    return out
if __name__ == "__main__":
    from gac import gac
    tag = sys.argv[1]; Ws = list(np.load(f"W_{tag}.npy").astype(np.float64)); mt = np.load(f"truth_{tag}.npz")["m"][-1].astype(np.float64)
    mse = lambda a: np.mean((a-mt)**2); sc = lambda v: 8*(mt @ (v-mt))/(mt @ mt)
    g = gac(Ws)[-1]["m"]; gt = gac_tree(Ws)[-1]["m"]; g3 = gac_tree(Ws, terms=("D3", "P3", "T3"))[-1]["m"]; g4 = gac_tree(Ws, terms=("D4", "PD4"))[-1]["m"]
    print(f"{tag}: GAC {mse(g):.3e} ({sc(g):+.4f}) | +trees {mse(gt):.3e} ({sc(gt):+.4f}) | +k3 only {mse(g3):.3e} ({sc(g3):+.4f}) | +k4 (D4,PD4) only {mse(g4):.3e} ({sc(g4):+.4f})", flush=True)
