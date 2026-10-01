"""Residual cut, round 2 (applied on top of patch_views.py): memoize the slices that
`_predict_core` / `_dslices` / `_hub2` take of PERMANENT pool buffers (the interleaved leg
slabs legs4["AP0"/"AP1"], the factor slabs fap4 / fap24, bufs["lap4"], the hub result, and
the rows of the term buffer abbuf) and hand those memoized views to the Strassen kernel, so
its quadrant / combo memo (round 1) also hits for the top-level operands.

    python patch_views2.py bundles/v29r/estimator.py bundles/v29r2/estimator.py

Only objects that ARE pool buffers (returned by _Pool.get / get_pair, alive for the whole
process) are used as memo bases; per-layer arrays (W, Qc, QU) are sliced fresh as before.
Same ops, same order, same data: FLOPs and outputs must be bit-identical.
"""
import sys

REPL = []

# helper on _Strassen
REPL.append(('''    def _quad(self, X, a, b):''', '''    def ps(self, base, key):
        """Memoized view of a PERMANENT pool buffer: key ("s", a, b) -> base[a:b],
        ("nn",) -> base[None], ("i", i) -> base[i]. The view is (re)registered for the
        derived-view memo above."""
        e = self._ext.get(id(base))
        if e is None or e[0] is not base:
            e = (base, {})
            self._ext[id(base)] = e
        r = e[1].get(key)
        if r is None:
            if key[0] == "s":
                r = base[key[1]:key[2]]
            elif key[0] == "nn":
                r = base[None]
            else:
                r = base[key[1]]
            e[1][key] = r
        v = self._vm.get(id(r))
        if v is None or v[0] is not r:
            self._reg(r)
        return r

    def _quad(self, X, a, b):'''))
REPL.append(('''        self._vm = {}   # residual cut: id(view) -> (view, {tag: derived views}), pooled views only
''', '''        self._vm = {}   # residual cut: id(view) -> (view, {tag: derived views}), pooled views only
        self._ext = {}  # residual cut (round 2): id(permanent pool buffer) -> (buffer, {key: view})
'''))

# call sites in _predict_core
REPL.append(('''                    smm.mm(Qc[None, None], fap4[2 * fa_off:2 * (fa_off + m1)],
                           legs4["AP1"][2 * kb:2 * ka], smm.level(n, r_old, n, s_sb))''',
             '''                    smm.mm(Qc[None, None], smm.ps(fap4, ("s", 2 * fa_off, 2 * (fa_off + m1))),
                           smm.ps(legs4["AP1"], ("s", 2 * kb, 2 * ka)), smm.level(n, r_old, n, s_sb))'''))
REPL.append(('''                    smm.mm(QU[None, None], fap24[:2 * kb], legs4["AP1"][:2 * kb],''',
             '''                    smm.mm(QU[None, None], smm.ps(fap24, ("s", 0, 2 * kb)), smm.ps(legs4["AP1"], ("s", 0, 2 * kb)),'''))
REPL.append(('''                    smm.mm(W[None, None], legs4["AP0"][2 * ka:2 * k + extra],
                           legs4["AP1"][2 * ka:2 * k + extra], smm.level(n, n, n, s_lev))''',
             '''                    smm.mm(W[None, None], smm.ps(legs4["AP0"], ("s", 2 * ka, 2 * k + extra)),
                           smm.ps(legs4["AP1"], ("s", 2 * ka, 2 * k + extra)), smm.level(n, n, n, s_lev))'''))
REPL.append(('''                    smm.mm(W[None, None], legs4["AP0"][0:2], legs4["AP1"][0:2],''',
             '''                    smm.mm(W[None, None], smm.ps(legs4["AP0"], ("s", 0, 2)), smm.ps(legs4["AP1"], ("s", 0, 2)),'''))
REPL.append(('''                                        sb1=(fap4[2 * fa_off:2 * (fa_off + ka - kb)] if ka > kb else None),
                                        sb2=(fap24[:2 * kb] if kb > 0 else None), s_sb=s_sb)''',
             '''                                        sb1=(smm.ps(fap4, ("s", 2 * fa_off, 2 * (fa_off + ka - kb))) if ka > kb else None),
                                        sb2=(smm.ps(fap24, ("s", 0, 2 * kb)) if kb > 0 else None), s_sb=s_sb)'''))
REPL.append(('''                        fnp.copyto(abbuf[t_], b2[b_])
                    elif b_ in ("ones2", "ones1"):
                        fnp.copyto(abbuf[t_], b2[a_])
                    else:
                        fnp.multiply(b2[a_], b2[b_], out=abbuf[t_])''',
             '''                        fnp.copyto(smm.ps(abbuf, ("i", t_)), b2[b_])
                    elif b_ in ("ones2", "ones1"):
                        fnp.copyto(smm.ps(abbuf, ("i", t_)), b2[a_])
                    else:
                        fnp.multiply(b2[a_], b2[b_], out=smm.ps(abbuf, ("i", t_)))'''))
# _hub2 and _dslices
REPL.append(('''        self._smm.hub(bufs["lap4"][2 * k0:2 * k], apb4[2 * k0:2 * k], out[None],''',
             '''        smm = self._smm
        smm.hub(smm.ps(bufs["lap4"], ("s", 2 * k0, 2 * k)), smm.ps(apb4, ("s", 2 * k0, 2 * k)), smm.ps(out, ("nn",)),'''))
REPL.append(('''                smm.hub(bufs["lap4"][2 * kb:2 * ka], sb1, inner, smm.level(n, n, r1, s_sb))''',
             '''                smm.hub(smm.ps(bufs["lap4"], ("s", 2 * kb, 2 * ka)), sb1, inner, smm.level(n, n, r1, s_sb))'''))
REPL.append(('''                smm.hub(bufs["lap4"][:2 * kb], sb2, inner2, smm.level(n, n, r2_, s_sb))''',
             '''                smm.hub(smm.ps(bufs["lap4"], ("s", 0, 2 * kb)), sb2, inner2, smm.level(n, n, r2_, s_sb))'''))


def main():
    src, dst = sys.argv[1:3]
    t = open(src).read()
    for old, new in REPL:
        assert t.count(old) == 1, old[:160]
        t = t.replace(old, new)
    import os
    os.makedirs(os.path.dirname(dst), exist_ok=True)
    open(dst, "w").write(t)


if __name__ == "__main__":
    main()
