# The chain's kappa3 readouts split by class against the exact transport law (note XXXVI section 3c).
#   python k3chain.py NET DUMP MCPREFIX
# DUMP: one-step dump with every oracle on and KEEP_EXTRA=K21,K3v (the chain's own D3/D21 recorded before replacement).
# At each l -> l+1:  chain readout = T3_pair(chain's y slices at l) + chain legs part.  Truth = T3_pair(y) + R3, and R3 =
# legs-1st-order(true T) + facet(true C) to 99% (mclegs.py). Reported for the D3 slice | D21 off-diagonal slice:
#   total readout error, pair-part error, legs-part error (all relative to the truth slice), and the legs-part error
#   regressed on the facet term and on the first-order legs term.
import sys, math, numpy as np
net, dmp, pre = int(sys.argv[1]), sys.argv[2], sys.argv[3]
ch = np.load(dmp); F = np.load(f"{pre}_full.npz"); G = np.load(f"{pre}_legs.npz")
Wcol = np.load(f"../official/W_off{net}.npy").astype(np.float64); L, n, _ = Wcol.shape; off = ~np.eye(n, dtype=bool)
mu = F["mu"].astype(np.float64); sd = np.sqrt(F["var"].astype(np.float64)); al = mu / sd
Phi = 0.5 * (1 + np.vectorize(math.erf)(al / math.sqrt(2))); phi = np.exp(-al * al / 2) / math.sqrt(2 * math.pi)
def offd(A):
    A = np.array(A, dtype=np.float64); np.fill_diagonal(A, 0.0); return A
def T3_pair(W, H, k3v, Dm):
    G3 = Dm @ W.T
    d3 = (W * H) @ k3v + 3 * np.einsum("ia,ai->i", H, G3)
    d21 = offd((H * k3v[None, :]) @ W.T + H @ Dm @ W.T + 2 * (W * G3.T) @ W.T)
    return d3, d21
def rel(a, b): return float(np.linalg.norm(a) / np.linalg.norm(b))
def fit(r, x):
    b = float(r.ravel() @ x.ravel() / (x.ravel() @ x.ravel())); return b, 1 - rel(r - b * x, r) ** 2
print(f"net {net}: chain kappa3 readouts (one-step, oracle inputs) by class; slices D3 | D21")
acc = []
for l in range(1, L - 1):
    need = [f"D3own_{l + 1}", f"D21own_{l + 1}", f"K3v_{l}", f"K21_{l}"]
    if any(k not in ch.files for k in need):
        continue
    W = Wcol[l + 1]; H = W * W; X = W * Phi[l][None, :]; XX = X * X
    k3z = F["k3"][l].astype(np.float64); D21z = offd(F["D21"][l]); Cz = offd(0.5 * (F["cov"][l] + F["cov"][l].T))
    t3 = (F["k3"][l + 1].astype(np.float64), offd(F["D21"][l + 1]))
    pair_true = T3_pair(W, H, F["k3_y"][l].astype(np.float64), offd(F["D21_y"][l]))
    R3 = tuple(a - b for a, b in zip(t3, pair_true))
    own = (ch[f"D3own_{l + 1}"].astype(np.float64), offd(ch[f"D21own_{l + 1}"]))
    pair_ch = T3_pair(W, H, ch[f"K3v_{l}"].astype(np.float64), offd(ch[f"K21_{l}"]))
    legs_ch = tuple(a - b for a, b in zip(own, pair_ch))
    # exact-law pieces from truth
    u3, P = G["u3"][l].astype(np.float64), G["P"][l].astype(np.float64)
    leg1 = (u3 - 3 * np.einsum("ia,ai->i", XX, D21z @ X.T) - (XX * X) @ k3z,
            offd(P - XX @ D21z @ X.T - 2 * (X * (D21z @ X.T).T) @ X.T - (XX * k3z[None, :]) @ X.T))
    rho = phi[l] / sd[l]; Wr = W * rho[None, :]; Gm = Cz @ X.T; C2 = Cz * Cz
    fac = (3 * np.einsum("ia,ia->i", Wr, Gm.T ** 2 - XX @ C2),
           offd((Wr @ (Gm * Gm - C2 @ XX.T)).T + 2 * ((Wr * Gm.T) @ Gm - (X * (Wr @ C2)) @ X.T)))
    out = []; rowacc = []
    for s in range(2):
        e_tot = own[s] - t3[s]; e_pair = pair_ch[s] - pair_true[s]; e_legs = legs_ch[s] - R3[s]
        bf, ef = fit(e_legs, -fac[s]); bl, el = fit(e_legs, leg1[s]); rowacc += [rel(e_tot, t3[s]), rel(e_legs, t3[s])]
        out.append(f"total {rel(e_tot, t3[s]):.4f} pair {rel(e_pair, t3[s]):.4f} legs {rel(e_legs, t3[s]):.4f} "
                   f"(legs part/R3 err {rel(e_legs, R3[s]):.3f}); legs err ~ -facet: coef {bf:+.2f} expl {ef:+.2f}, ~ leg1: coef {bl:+.2f} expl {el:+.2f}")
    print(f"layer {l:2d}->{l + 1:2d}: D3: {out[0]}\n             D21: {out[1]}", flush=True)
    acc.append(rowacc)
a = np.array(acc)
print("SUMMARY mean over layers: D3 total %.4f legs %.4f | D21 total %.4f legs %.4f; layers 10-14: D3 legs %.4f D21 legs %.4f"
      % (*a.mean(0), a[-5:, 1].mean(), a[-5:, 3].mean()))
