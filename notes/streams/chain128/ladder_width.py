"""Interface ladder vs width: eps of D21(z_{l+1}) when the ReLU step is fed an atlas's TRUE z_l cumulants and the
all-distinct kappa3(a_l) is built by each rule (slices of kappa3(a) from the step itself), noise-corrected with an
independent atlas of the same MLP when given.
    python ladder_width.py A.npz [B.npz] > results/ladder_w<n>.txt
Rules: none (slices only), wick (Gaussian rho^2 + Phi^3 kappa3), closure (oracle leg-partition coefficients),
engine (all first-order diagrams), engine_reg211 / engine_zero211 (kappa4 (2,1,1) slice regenerated u_i C_jk / zeroed),
true (the atlas's own kappa3(a): measures the step's slice error + noise).
"""
import sys
import numpy as np
import chain as ch
import oracle_k3 as ok


def rel(a, b):
    return float(np.sqrt(np.sum((a - b) ** 2) / np.sum(b ** 2)))


def main():
    A = np.load(sys.argv[1]); B = np.load(sys.argv[2]) if len(sys.argv) > 2 else None
    W = A["weights"].astype(np.float64)
    L, n, _ = W.shape
    idx = np.arange(n)
    rules = ["none", "wick", "closure", "engine", "engine_reg211", "engine_zero211", "true"]
    print(f"# {sys.argv[1]} width {n} depth {L} N {int(A['n_samples'])}" + (" (noise-corrected with pair)" if B is not None else ""))
    print("l  noise " + " ".join(f"{r:>14}" for r in rules))
    for l in range(L - 1):
        st = ch.atlas_state(A, l)
        T = ok.offdiag(ch.atlas_state(A, l + 1, with_k4=False).K3[idx, idx, :])
        nz = 0.0
        if B is not None:
            TB = ok.offdiag(ch.atlas_state(B, l + 1, with_k4=False).K3[idx, idx, :])
            nz = rel(T, TB) / np.sqrt(2)
            T = 0.5 * (T + TB); nz = nz / np.sqrt(2)   # target = mean of the two atlases
        out = []
        for r in rules:
            s = ch.State(st.mu, st.C, st.K3, st.X)
            mode = r
            if r.startswith("engine_"):
                mode = "engine"
                K211 = ok.all_distinct(st.X)
                if r == "engine_zero211":
                    s.X = st.X - K211
                else:
                    Co = ok.offdiag(st.C)
                    u = np.einsum("ijk,jk->i", K211, Co) / float(np.sum(Co * Co))
                    s.X = st.X - K211 + ok.all_distinct(np.einsum("i,jk->ijk", u, Co))
            if r == "true":
                _, _, K3a = ch.atlas_post(A, l)
            else:
                _, _, K3a, _, _ = ch.relu_step(s, mode, want_k4=False)
            e = rel(ok.offdiag(ok.transport_d21(K3a, W[l + 1])), T)
            out.append(np.sqrt(max(e * e - nz * nz, 0.0)))
        print(f"{l:<2} {nz:5.3f} " + " ".join(f"{v:14.3f}" for v in out), flush=True)


if __name__ == "__main__":
    main()
