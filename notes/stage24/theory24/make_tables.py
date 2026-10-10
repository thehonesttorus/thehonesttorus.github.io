"""Collect the Stage 24 width-1024 measurements into LaTeX tables (notes/stage24/s24/tables.tex).

  python notes/stage24/theory24/make_tables.py KIKUCHI_DIR ELLIPSOID_DIR"""
import sys, os, re, glob, numpy as np

KD, ED = sys.argv[1], sys.argv[2]; OUT = os.path.join(os.path.dirname(__file__), "..", "s24", "tables.tex")
num = lambda s: [float(x) for x in re.findall(r"[-+]?\d+\.\d+(?:e[-+]?\d+)?|[-+]?\d+", s)]


def parse(path):
    """Per-layer dicts from a s24_kikuchi text file (fields split on ' | ')."""
    t1, t2, cur = {}, {}, None
    for line in open(path):
        if line.startswith(" l | level"): cur = 1; continue
        if line.startswith(" l | excitation"): cur = 2; continue
        if cur and re.match(r"^\s*\d+ \|", line):
            f = [x.strip() for x in line.split(" | ")]; l = int(f[0]) - 1
            if cur == 1:
                lv = num(f[1]); r2 = num(f[3]); fl = num(f[4]); gpr = num(f[7]); et = num(f[8]); gp = num(f[10])
                lg = re.findall(r"gap ([\d.]+) \(16/128 PCs out: ([\d.]+) ([\d.]+)\)", line)[0]
                t1[l] = dict(r=lv[0], rsd=lv[1], r2s=r2[0] if r2 else np.nan, r2v=r2[1] if r2 else np.nan, rate=fl[0], het=fl[1], ud=fl[2],
                             rice=num(f[6])[0], eta=gpr[0], pr=gpr[1], eta16=et[5], eta128=et[8], g=gp[0], prod=gp[1],
                             lin=float(lg[0]), lin16=float(lg[1]), lin128=float(lg[2]))
            else:
                ex = num(f[1]); rel = num(f[2]); xs = num(f[4].split("(edge")[0])
                t2[l] = dict(ex=ex[0], exsd=ex[1], expred=ex[2], relax=rel[0], x1=xs[0] if xs else np.nan)
    return t1, t2


files = sorted(f for f in glob.glob(f"{KD}/s24_kikuchi_net*.txt") if "_N17" not in f)
P = {int(re.findall(r"net(\d+)", f)[0]): parse(f) for f in files}; nets = sorted(P)
rg = lambda v: f"{min(v):.3f}--{max(v):.3f}"
L = []
L += [r"\begin{table}[ht]\centering\footnotesize",
      r"\caption{Levels and walls per layer, width $1024$ ($N=2^{16}$ inputs and rotation partners per net). Net 0 values;"
      r" net 0 only (the multi-net run was stopped).}\label{tab:levels}",
      r"\resizebox{\textwidth}{!}{\begin{tabular}{@{}rccccccc@{}}\toprule",
      r"$\ell$ & $r/n$ & sd$(r)/n$ & excitation$/n$ & $R^2$ next level: signs / values & wall rate & rate sd/mean & Rice--Gauss / measured\\\midrule"]
for l in range(16):
    a = P[nets[0]][0][l]; b = P[nets[0]][1][l]
    rate = [P[k][0][l]["rate"] for k in nets]; rice = [P[k][0][l]["rice"] for k in nets]; ex = [P[k][1][l]["ex"] for k in nets]
    r2 = f"{a['r2s']:.2f} / {a['r2v']:.2f}" if not np.isnan(a["r2s"]) else "--"
    L.append(f"{l + 1} & {a['r']:.3f} & {a['rsd']:.4f} & {b['ex']:.3f} [{rg(ex)}] & {r2} & {a['rate']:.3f} [{rg(rate)}] & {a['het']:.2f} & {a['rice']:.3f} [{rg(rice)}]\\\\")
L += [r"\bottomrule\end{tabular}}\end{table}"]
L += [r"\begin{table}[ht]\centering\footnotesize",
      r"\caption{Dependence of the gate law, single-layer wall-chain bounds (per unit rotation angle) and depth gains, net 0.}\label{tab:gap}",
      r"\resizebox{\textwidth}{!}{\begin{tabular}{@{}rcccccccc@{}}\toprule",
      r"$\ell$ & $\eta$ & $\eta$, 16 / 128 PCs out & PR$(\corr g)$ & $\lambda_{\rm lin}$ & $\lambda_{\rm lin}$, 16 / 128 PCs out & level relax. & $\|\corr(g_\ell,g_{\ell+1})\|$ & $\prod_{k>\ell}g_k$\\\midrule"]
for l in range(16):
    m = lambda key, t=0: np.nanmean([P[k][t][l][key] for k in nets])
    eta = [P[k][0][l]["eta"] for k in nets]; lin = [P[k][0][l]["lin"] for k in nets]
    x1 = f"{m('x1', 1):.2f}" if l < 15 else "--"
    L.append(f"{l + 1} & {m('eta'):.2f} [{min(eta):.1f}--{max(eta):.1f}] & {m('eta16'):.2f} / {m('eta128'):.2f} & {m('pr'):.0f} & "
             f"{m('lin'):.3f} [{rg(lin)}] & {m('lin16'):.2f} / {m('lin128'):.2f} & {m('relax', 1):.2f} & {x1} & {m('prod'):.3f}\\\\")
L += [r"\bottomrule\end{tabular}}\end{table}"]
ell = sorted(glob.glob(f"{ED}/s24_ellipsoid_net*.txt"))
if ell:
    L += [r"\begin{table}[ht]\centering\footnotesize",
          r"\caption{The terminal ellipsoid: visible spectrum of $S=\cov F$ and the share of first-pass errors in the eigenvector bands"
          r" $[0,16),[16,64),[64,256),[256,512),[512,1024)$; R = Rayleigh ratio $e\T Se/(|e|^2\overline\lambda)$.}\label{tab:ell}",
          r"\resizebox{\textwidth}{!}{\begin{tabular}{@{}rcccll@{}}\toprule",
          r"net & PR & trace, top 16/64/256 & bottom half & closure: bands; R & CC1: bands; R\\\midrule"]
    for f in ell:
        t = open(f).read(); k = int(re.findall(r"net (\d+)", t)[0])
        pr = float(re.findall(r"PR ([\d.]+)", t)[0]); top = re.findall(r"top 1/4/16/64/256: ([\d./]+)", t)[0].split("/")
        bh = float(re.findall(r"bottom half ([\d.]+)", t)[0])
        rows = {m_.group(1): (m_.group(2).split(), float(m_.group(3))) for m_ in
                re.finditer(r"b=(\w+)\s+rms [\d.e+-]+: error share in bands ([\d. ]+)\|.*?Rayleigh ratio ([\d.]+)", t)}
        cell = lambda b_: " ".join(rows[b_][0]) + f"; {rows[b_][1]:.1f}" if b_ in rows else "--"
        L.append(f"{k} & {pr:.1f} & {top[2]}/{top[3]}/{top[4]} & {bh:.4f} & {cell('closure')} & {cell('CC1')}\\\\")
    L += [r"\bottomrule\end{tabular}}\end{table}"]
open(OUT, "w").write("\n".join(L) + "\n"); print("wrote", OUT, "nets", nets)
