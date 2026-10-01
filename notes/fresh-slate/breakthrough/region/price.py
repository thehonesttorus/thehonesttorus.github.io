"""Price the local errors of the exact Gaussian closure with the transfer coefficients (REPORT §3) for one MLP."""
import json, sys, numpy as np
sys.path.insert(0, '../../bench')
import bench, gclose

def main(i):
    rows = json.load(open(f'results/mc_mc1024_mlp{i}_N262144.json')); lr = json.load(open(f'results/lr_mlp{i}.json'))
    D = np.load(f'data/mc1024_mlp{i}_N262144.npz')
    S = bench.load_set('w1024_d16'); W = bench.weights(S, i)
    meas = float(((gclose.predict_exact(W)[-1] - S['means'][i][-1]) ** 2).mean() - S['noise'][i])
    tot = dict(off=0., diag=0., mean=0., off_after21=0., off_after21_22=0., mean_after_edgeworth=0.)
    for r in rows:
        l = r['layer']
        var = np.diag(D['A_s2t'][l]) - D['A_s1t'][l] ** 2; noise = var.mean() / S['n_samples']
        K = lambda k: (lr[k][l] if l < 15 else (1.0 if k == 'mean' else 0.0))
        tot['off'] += K('off') * r['fro2_dCoff']; tot['diag'] += K('diag') * r['rms2_dvar']
        tot['mean'] += K('mean') * (r['rms2_dm'] - noise)
        tot['off_after21'] += K('off') * r['fro2_res_G21']; tot['off_after21_22'] += K('off') * r['fro2_res_G21+31+G22']
        tot['mean_after_edgeworth'] += K('mean') * (r['rms2_dm_res'] - noise)
    tot['sum'] = tot['off'] + tot['diag'] + tot['mean']; tot['measured_gauss_raw'] = meas
    print(i, json.dumps({k: f'{v:.3e}' for k, v in tot.items()}))
    return tot

if __name__ == '__main__':
    out = {int(a): main(int(a)) for a in sys.argv[1:]}
    json.dump(out, open('results/price_' + '_'.join(sys.argv[1:]) + '.json', 'w'), indent=1)
