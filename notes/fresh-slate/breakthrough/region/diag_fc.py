"""Compare an FC run's per-layer state with the n = 1024 MC atlas: eps of D21 (cross-half noise-corrected), covariance
error of the chain, truth-noise floor per layer."""
import sys, json, numpy as np
sys.path.insert(0, '../../bench')
import bench, fc, mc_analyse as ma

def main(i, path, **kw):
    S = bench.load_set('w1024_d16'); W = bench.weights(S, i); T = S['means'][i]
    if kw.pop('k4tf', False):
        DD = np.load(path); k4f = {}
        for l in range(1, W.shape[0]):
            cA = ma.central(DD, 'A', l); cB = ma.central(DD, 'B', l)
            k4f[l] = (0.5 * (cA['K22'] + cB['K22']), 0.5 * (cA['K31'] + cB['K31']), 0.5 * (cA['k4'] + cB['k4']))
        mode = kw.pop('k4mode', 'full')
        for l in list(k4f):
            K22, K31, k4 = k4f[l]
            if mode == 'nooff':
                K22 = 0 * K22; K31 = 0 * K31
            elif mode.startswith('colmean') and mode != 'colmean':   # colmean*x / colmean_no31 / colmean_no22
                K22 = np.broadcast_to(K22.mean(0, keepdims=True), K22.shape).copy(); K31 = np.broadcast_to(K31.mean(0, keepdims=True), K31.shape).copy()
                if '*' in mode:
                    x = float(mode.split('*')[1]); K22 = x * K22; K31 = x * K31
                if mode.endswith('no31'):
                    K31 = 0 * K31
                if mode.endswith('no22'):
                    K22 = 0 * K22
            elif mode == 'colmean':      # keep only the coherent column sums: K_ik -> mean_i' K_i'k
                K22 = np.broadcast_to(K22.mean(0, keepdims=True), K22.shape).copy(); K31 = np.broadcast_to(K31.mean(0, keepdims=True), K31.shape).copy()
            elif mode == 'diagonly':
                K22 = 0 * K22; K31 = 0 * K31
            elif mode.startswith('scale'):
                x = float(mode[5:]); K22 = x * K22; K31 = x * K31
            elif mode == 'nodiag':
                k4 = 0 * k4
            k4f[l] = (K22, K31, k4)
        kw['k4f'] = k4f
    tr = []; p = fc.run(W, trace=tr, **kw)
    D = np.load(path)
    out = []
    for l in range(W.shape[0]):
        cA = ma.central(D, 'A', l); cB = ma.central(D, 'B', l)
        varA = np.diag(D['A_s2t'][l]) - D['A_s1t'][l] ** 2
        noise = varA.mean() / S['n_samples']
        mse = ((p[l] - T[l]) ** 2).mean()
        r = dict(l=l, mse_excess=mse - noise, noise=noise)
        if tr[l]['D21'] is not None:
            M = tr[l]['D21']
            KA, KB = cA['K21'], cB['K21']
            num = np.sum((M - KA) * (M - KB)); den = np.sum(KA * KB)
            r['eps_D21'] = float(np.sqrt(max(num, 0) / den))
            r['eps_D3'] = float(np.sqrt(max(np.sum((np.diag(M) - cA['k3']) * (np.diag(M) - cB['k3'])), 0) / np.sum(cA['k3'] * cB['k3'])))
        # pre-activation covariance error vs MC (noise-corrected)
        Sm = tr[l]['S']
        num = np.sum((Sm - cA['C']) * (Sm - cB['C'])); den = np.sum(cA['C'] * cB['C'])
        r['eps_S'] = float(np.sqrt(max(num, 0) / den))
        dv = np.diag(Sm) - 0.5 * (cA['var'] + cB['var'])
        r['rel_var_rms'] = float(np.sqrt(np.mean((dv / np.diag(Sm)) ** 2)))
        out.append(r)
        print(json.dumps({k: (f'{v:.3e}' if isinstance(v, float) else v) for k, v in r.items()}), flush=True)
    return out

if __name__ == '__main__':
    import ast
    kw = ast.literal_eval(sys.argv[3]) if len(sys.argv) > 3 else {}
    main(int(sys.argv[1]), sys.argv[2], **kw)
