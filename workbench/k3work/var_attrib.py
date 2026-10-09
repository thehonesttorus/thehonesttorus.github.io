# Which of the chain's two-site slices carries its per-layer variance injection? (note XLIII section 3)
#   python var_attrib.py LOCDIR MC2DIR CD2DIR NET "LAYERS"
# Per target layer t (source s = t - 1), every first-order Edgeworth term of the post-activation covariance, read
# through the next layer's rows, diag W_t D W_t^T / var_true_t, computed twice: from the chain's own slices at its own
# state (chaindump2: D3, D21, g4row, wk4m, wk431; C_off, var; means W_s out_(s-1)) and from the Monte Carlo cumulants at
# the true state (halves A, B for noise-free statistics). Terms: 'D21' kappa(a,a,b); 'k3' kappa3(a) (diagonal and (3,0));
# 'K22' kappa(a,a,b,b); 'K31' kappa(a,a,a,b); 'k4' kappa4(a) (diagonal and (4,0)).
# Checks first that the chain's own correction x = var_chain_t - diag W KG(chain state) W^T equals the sum of its own
# terms (the dumps are what its Wick stage reads), then splits the injection x - g term by term.
import sys, numpy as np
sys.argv, _argv = sys.argv[:1] + [".", ".", "0", "1"], sys.argv
exec(open(__file__.replace("var_attrib.py", "var_ladder.py")).read().split("cd = np.load")[0])
sys.argv = _argv
loc, mcd, cdd, net = sys.argv[1], sys.argv[2], sys.argv[3], int(sys.argv[4])
layers = [int(v) for v in sys.argv[5].split(",")]
TERMS = ("D21", "k3", "K22", "K31", "k4")


def terms(mu, C, k3, k4, D21, K22, K31):
    v = np.diag(C).copy(); s = np.sqrt(v); al = mu / s
    R = C / np.outer(s, s); np.fill_diagonal(R, 0.0)
    c = ccoef(al, M + 5)
    sa = s[:, None]; sb = s[None, :]
    _, m1 = relu_var(mu, v)
    out = {}
    T = 0.5 * D21 / sa * series(c, c, R, 2, 1, 0); T = T + T.T; np.fill_diagonal(T, 0.0); out["D21"] = T
    T = (k3[:, None] / 6.0) * sa ** -2 * sb * series(c, c, R, 3, 0, 1); T = T + T.T
    T[np.diag_indices_from(T)] = k3 / 6.0 * 2.0 * c[2] / s - 2.0 * m1 * k3 / 6.0 * c[3] / s ** 2; out["k3"] = T
    T = 0.25 * K22 / (sa * sb) * series(c, c, R, 2, 2, 0); np.fill_diagonal(T, 0.0); out["K22"] = T
    T = (K31 / 6.0) * sa ** -2 * series(c, c, R, 3, 1, 0); T = T + T.T; np.fill_diagonal(T, 0.0); out["K31"] = T
    T = (k4[:, None] / 24.0) * sa ** -3 * sb * series(c, c, R, 4, 0, 1); T = T + T.T
    T[np.diag_indices_from(T)] = k4 / 24.0 * 2.0 * c[3] / s ** 2 - 2.0 * m1 * k4 / 24.0 * c[4] / s ** 3; out["k4"] = T
    return out


cd = np.load(f"{loc}/chaindump_{net}.npz"); c2 = np.load(f"{cdd}/chaindump2_{net}.npz")
W = np.load(f"{loc}/W_off{net}.npy").astype(np.float64)
TT = {h: np.load(f"{mcd}/mc2_off{net}_{h}.npz") for h in ("full", "h0", "h1")}
out = c2["out"]
print(f"=== network {net}: diag W_t (term) W_t^T / var_true_t; chain = own slices at own state; true = MC (nf = noise-free)",
      flush=True)
for t in layers:
    s_ = t - 1
    mu_c = W[s_] @ out[s_ - 1]
    Cc = c2[f"C_off_{s_}"].astype(np.float64); Cc = 0.5 * (Cc + Cc.T); np.fill_diagonal(Cc, c2[f"var_{s_}"])
    z = lambda k: c2[f"{k}_{s_}"].astype(np.float64)
    D21c = z("D21"); np.fill_diagonal(D21c, 0.0)
    K22c = z("wk4m"); K22c = 0.5 * (K22c + K22c.T); np.fill_diagonal(K22c, 0.0)
    K31c = z("wk431").T.copy(); np.fill_diagonal(K31c, 0.0)
    Tc = terms(mu_c, Cc, z("D3"), z("g4row"), D21c, K22c, K31c)
    x = cd["var"][t] - qdiag(W[t], KG(mu_c, Cc))
    ec = {k: qdiag(W[t], v) for k, v in Tc.items()}
    res = {}
    for h in ("full", "h0", "h1"):
        F = TT[h]
        mu = F["mu"][s_].astype(np.float64); C = F["cov"][s_].astype(np.float64); C = 0.5 * (C + C.T)
        D21 = F["D21"][s_].astype(np.float64).copy(); np.fill_diagonal(D21, 0.0)
        K22 = F["K22"][s_].astype(np.float64); K22 = 0.5 * (K22 + K22.T); np.fill_diagonal(K22, 0.0)
        K31 = F["K31"][s_].astype(np.float64).copy(); np.fill_diagonal(K31, 0.0)
        Tt = terms(mu, C, F["k3"][s_].astype(np.float64), F["k4"][s_].astype(np.float64), D21, K22, K31)
        g = F["var"][t].astype(np.float64) - qdiag(W[t], KG(mu, C))
        res[h] = dict(g=g, vt=F["var"][t].astype(np.float64), **{k: qdiag(W[t], v) for k, v in Tt.items()})
    ref = res["full"]["vt"]; A, B = res["h0"], res["h1"]
    sumc = sum(ec.values())
    print(f" layer {t}: chain check: rms(x - sum of own terms) {np.sqrt(np.mean(((x - sumc) / ref) ** 2)):.2e} "
          f"(rms x {np.sqrt(np.mean((x / ref) ** 2)):.2e})", flush=True)
    yA, yB = x - A["g"], x - B["g"]
    e_y = np.mean(yA * yB / ref ** 2)
    line = []
    for k in TERMS:
        dA, dB = ec[k] - A[k], ec[k] - B[k]          # chain minus truth, this term
        nf = np.sqrt(max(np.mean(dA * dB / ref ** 2), 0.0))
        sh = 0.5 * np.mean((dA * yB + dB * yA) / ref ** 2) / e_y     # share of the injection energy (projection)
        line.append(f"{k}: chain {np.sqrt(np.mean((ec[k] / ref) ** 2)):.1e} true {np.sqrt(max(np.mean(A[k] * B[k] / ref ** 2), 0)):.1e} "
                    f"diff nf {nf:.1e} share {sh:+.2f}")
    print("   " + "\n   ".join(line), flush=True)
    rA = yA - sum(ec[k] - A[k] for k in TERMS); rB = yB - sum(ec[k] - B[k] for k in TERMS)
    print(f"   injection nf {np.sqrt(max(e_y, 0)):.2e}; left after all term differences nf "
          f"{np.sqrt(max(np.mean(rA * rB / ref ** 2), 0)):.2e}", flush=True)
