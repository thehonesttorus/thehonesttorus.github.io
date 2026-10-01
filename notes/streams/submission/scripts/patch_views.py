"""Residual cut for V29 (applied on top of patch_safe.py): memoize the VIEWS of the pooled
Strassen buffers instead of re-slicing them on every call.

    python patch_views.py bundles/v29/estimator.py bundles/v29r/estimator.py

Why: on the grader every flopscope call -- including a free basic slice -- is a client/server
round trip whose client-side bookkeeping (wrapper preamble, RemoteArray construction, its
weakref finaliser) lands partly in the residual clock. A steady-state V29 predict makes ~28k
calls, ~14.6k of them slices, and most of those re-slice the same persistent pooled buffers
(combo columns, quadrants, batch prefixes, swapaxes, [None], [0]).

What: `_Strassen` keeps (a) a slice cache for the pooled buffers it hands out and (b) a
registry of those persistent view objects with their derived views (quadrants, the 7 combo
columns, swapaxes, [None], [0]). Only registered objects (strong refs held, identity
checked) get memoized; any other operand is sliced fresh exactly as before. Any pool
growth clears both caches (old buffers are released, views rebuilt lazily).

Same ops on the same data in the same order: FLOPs and outputs must be bit-identical
(verified per MLP, see REPORT.md).
"""
import sys

REPL = []

# 1. caches + helpers
REPL.append(('''    def __init__(self, dtype, bmax):
        self.dtype = dtype
        self.bmax = int(bmax)
        self.pool = {}
''', '''    def __init__(self, dtype, bmax):
        self.dtype = dtype
        self.bmax = int(bmax)
        self.pool = {}
        self._sc = {}   # residual cut: cached slices of pooled buffers
        self._vm = {}   # residual cut: id(view) -> (view, {tag: derived views}), pooled views only

    def _reset_views(self):
        self._sc.clear()
        self._vm.clear()

    def _reg(self, o):
        self._vm[id(o)] = (o, {})
        return o

    def _derived(self, X, tag, make):
        e = self._vm.get(id(X))
        if e is None or e[0] is not X:
            return make()
        r = e[1].get(tag)
        if r is None:
            r = make()
            if isinstance(r, tuple):
                for o in r:
                    self._reg(o)
            else:
                self._reg(r)
            e[1][tag] = r
        return r

    def _quad(self, X, a, b):
        """(X[..,:a,:b], X[..,:a,b:], X[..,a:,:b], X[..,a:,b:]), memoized for pooled views."""
        return self._derived(X, ("q", a, b), lambda: (X[..., :a, :b], X[..., :a, b:],
                                                      X[..., a:, :b], X[..., a:, b:]))

    def _bsel(self, key, shape, sel):
        """self._buf(key, shape)[:sel] (sel int >= 1) or [0] (sel == 0 with shape[0] == 1)."""
        buf = self._buf(key, shape)
        k = (key, sel)
        s = self._sc.get(k)
        if s is None or s[0] is not buf:
            s = (buf, self._reg(buf[0] if sel == 0 else buf[:sel]))
            self._sc[k] = s
        return s[1]
'''))

REPL.append(('''        b = self.pool.get(key)
        if b is None or b.shape[0] < shape[0]:
            b = fnp.empty(shape, dtype=self.dtype)
            fnp.copyto(b, 0.0)   # metered page touch (first-MLP residual, F80)
            self.pool[key] = b
        return b
''', '''        b = self.pool.get(key)
        if b is None or b.shape[0] < shape[0]:
            b = fnp.empty(shape, dtype=self.dtype)
            fnp.copyto(b, 0.0)   # metered page touch (first-MLP residual, F80)
            self.pool[key] = b
            self._reset_views()
        return b
'''))

# 2. _buf2: cached, registered batch prefixes
REPL.append(('''        e = self.pool.get(key)
        b = int(shape5[0])
        if e is None or e[0].shape[0] < b:
            b5 = fnp.empty(shape5, dtype=self.dtype)
            fnp.copyto(b5, 0.0)   # metered page touch (first-MLP residual, F80)
            b4 = fnp.reshape(b5, (b, shape5[1] * shape5[2], shape5[3], shape5[4]))
            e = (b5, b4)
            self.pool[key] = e
        return e[0][:b], e[1][:b]
''', '''        e = self.pool.get(key)
        b = int(shape5[0])
        if e is None or e[0].shape[0] < b:
            b5 = fnp.empty(shape5, dtype=self.dtype)
            fnp.copyto(b5, 0.0)   # metered page touch (first-MLP residual, F80)
            b4 = fnp.reshape(b5, (b, shape5[1] * shape5[2], shape5[3], shape5[4]))
            e = (b5, b4)
            self.pool[key] = e
            self._reset_views()
        k = ("b2", key, b)
        s = self._sc.get(k)
        if s is None or s[0] is not e[0]:
            s = (e[0], self._reg(e[0][:b]), self._reg(e[1][:b]))
            self._sc[k] = s
        return s[1], s[2]
'''))

# 3. _combos: memoized combo columns
REPL.append(('''        buf, buf4 = self._buf2(key, (b, 7, P, h, q))
        X11, X12, X21, X22 = Q
        if kind == "L":
            fnp.add(X11, X22, out=buf[:, 0]); fnp.add(X21, X22, out=buf[:, 1])
            fnp.copyto(buf[:, 2], X11);       fnp.copyto(buf[:, 3], X22)
            fnp.add(X11, X12, out=buf[:, 4]); fnp.subtract(X21, X11, out=buf[:, 5])
            fnp.subtract(X12, X22, out=buf[:, 6])
        else:
            fnp.add(X11, X22, out=buf[:, 0]); fnp.copyto(buf[:, 1], X11)
            fnp.subtract(X12, X22, out=buf[:, 2]); fnp.subtract(X21, X11, out=buf[:, 3])
            fnp.copyto(buf[:, 4], X22);       fnp.add(X11, X12, out=buf[:, 5])
            fnp.add(X21, X22, out=buf[:, 6])
        return buf4
''', '''        buf, buf4 = self._buf2(key, (b, 7, P, h, q))
        X11, X12, X21, X22 = Q
        c = self._derived(buf, "cols", lambda: tuple(buf[:, i] for i in range(7)))
        if kind == "L":
            fnp.add(X11, X22, out=c[0]); fnp.add(X21, X22, out=c[1])
            fnp.copyto(c[2], X11);       fnp.copyto(c[3], X22)
            fnp.add(X11, X12, out=c[4]); fnp.subtract(X21, X11, out=c[5])
            fnp.subtract(X12, X22, out=c[6])
        else:
            fnp.add(X11, X22, out=c[0]); fnp.copyto(c[1], X11)
            fnp.subtract(X12, X22, out=c[2]); fnp.subtract(X21, X11, out=c[3])
            fnp.copyto(c[4], X22);       fnp.add(X11, X12, out=c[5])
            fnp.add(X21, X22, out=c[6])
        return buf4
'''))

# 4. _assemble: memoized out quadrants and product columns
REPL.append(('''    @staticmethod
    def _assemble(M, out):
        """M (b, 7, P, h, w) products -> out (b, P, m, w) quadrants (8 ops)."""
        h, w = M.shape[3], M.shape[4]
        C11, C12 = out[..., :h, :w], out[..., :h, w:]
        C21, C22 = out[..., h:, :w], out[..., h:, w:]
        fnp.add(M[:, 0], M[:, 3], out=C11); fnp.subtract(C11, M[:, 4], out=C11)
        fnp.add(C11, M[:, 6], out=C11)
        fnp.add(M[:, 2], M[:, 4], out=C12)
        fnp.add(M[:, 1], M[:, 3], out=C21)
        fnp.subtract(M[:, 0], M[:, 1], out=C22); fnp.add(C22, M[:, 2], out=C22)
        fnp.add(C22, M[:, 5], out=C22)
''', '''    def _assemble(self, M, out):
        """M (b, 7, P, h, w) products -> out (b, P, m, w) quadrants (8 ops)."""
        h, w = M.shape[3], M.shape[4]
        C11, C12, C21, C22 = self._quad(out, h, w)
        Mc = self._derived(M, "cols", lambda: tuple(M[:, i] for i in range(7)))
        fnp.add(Mc[0], Mc[3], out=C11); fnp.subtract(C11, Mc[4], out=C11)
        fnp.add(C11, Mc[6], out=C11)
        fnp.add(Mc[2], Mc[4], out=C12)
        fnp.add(Mc[1], Mc[3], out=C21)
        fnp.subtract(Mc[0], Mc[1], out=C22); fnp.add(C22, Mc[2], out=C22)
        fnp.add(C22, Mc[5], out=C22)
'''))

# 5. mm
REPL.append(('''        h, q, v = m // 2, kd // 2, w // 2
        XQ = (X[..., :h, :q], X[..., :h, q:], X[..., h:, :q], X[..., h:, q:])
        YQ = (Y[..., :q, :v], Y[..., :q, v:], Y[..., q:, :v], Y[..., q:, v:])
        if P >= STRASSEN_FUSE_P and not self._ok(h, q, v, lev - 1):
            # V28: fused leaf (deep levels only) -- one product at a time into the out quadrants
            self._leaf_mm(XQ, YQ, (out[..., :h, :v], out[..., :h, v:], out[..., h:, :v], out[..., h:, v:]),
                          bx, by, P, h, q, v)
            return
''', '''        h, q, v = m // 2, kd // 2, w // 2
        XQ = self._quad(X, h, q)
        YQ = self._quad(Y, q, v)
        if P >= STRASSEN_FUSE_P and not self._ok(h, q, v, lev - 1):
            # V28: fused leaf (deep levels only) -- one product at a time into the out quadrants
            self._leaf_mm(XQ, YQ, self._quad(out, h, v),
                          bx, by, P, h, q, v)
            return
'''))

REPL.append(('''        Lb = self._buf(("L1", P, h, q), (bx, P, h, q))[:bx]
        Rb = self._buf(("R1", P, q, v), (by, P, q, v))[:by]
        Mb = self._buf(("M1", P, h, v), (by, P, h, v))[:by]
''', '''        Lb = self._bsel(("L1", P, h, q), (bx, P, h, q), bx)
        Rb = self._bsel(("R1", P, q, v), (by, P, q, v), by)
        Mb = self._bsel(("M1", P, h, v), (by, P, h, v), by)
'''))

# 6. hub
REPL.append(('''            T = self._buf(("K", P, m, w), (k, P, m, w))[:k]
            fnp.matmul(X, fnp.swapaxes(Y, -1, -2), out=T)
            fnp.sum(T, axis=0, out=out)
            return
        h, q, v = m // 2, kd // 2, w // 2
        XQ = (X[..., :h, :q], X[..., :h, q:], X[..., h:, :q], X[..., h:, q:])
        # (Y^T) blocks expressed on the un-transposed Y (rows c, cols j)
        YQ = (Y[..., :v, :q], Y[..., v:, :q], Y[..., :v, q:], Y[..., v:, q:])
        if P >= STRASSEN_FUSE_P and not self._ok(h, q, v, lev - 1):
            self._leaf_hub(XQ, YQ, (out[..., :h, :v], out[..., :h, v:], out[..., h:, :v], out[..., h:, v:]),
                           k, P, h, q, v)
            return
        Xc = self._combos(X, XQ, "L", ("L", P, h, q))
        Yc = self._combos(Y, YQ, "R", ("R", P, v, q))
        Mb, Mb4 = self._buf2(("HM", P, h, v), (1, 7, P, h, v))
        self.hub(Xc, Yc, Mb4[0], lev - 1)
        self._assemble(Mb, out[None])
''', '''            T = self._bsel(("K", P, m, w), (k, P, m, w), k)
            fnp.matmul(X, self._derived(Y, "sw", lambda: fnp.swapaxes(Y, -1, -2)), out=T)
            fnp.sum(T, axis=0, out=out)
            return
        h, q, v = m // 2, kd // 2, w // 2
        XQ = self._quad(X, h, q)
        # (Y^T) blocks expressed on the un-transposed Y (rows c, cols j)
        g = self._quad(Y, v, q)
        YQ = (g[0], g[2], g[1], g[3])
        if P >= STRASSEN_FUSE_P and not self._ok(h, q, v, lev - 1):
            self._leaf_hub(XQ, YQ, self._quad(out, h, v),
                           k, P, h, q, v)
            return
        Xc = self._combos(X, XQ, "L", ("L", P, h, q))
        Yc = self._combos(Y, YQ, "R", ("R", P, v, q))
        Mb, Mb4 = self._buf2(("HM", P, h, v), (1, 7, P, h, v))
        self.hub(Xc, Yc, self._derived(Mb4, "i0", lambda: Mb4[0]), lev - 1)
        self._assemble(Mb, self._derived(out, "none", lambda: out[None]))
'''))

REPL.append(('''        Lb = self._buf(("L1", P, h, q), (k, P, h, q))[:k]
        Rb = self._buf(("R1", P, v, q), (k, P, v, q))[:k]
        T = self._buf(("K", P, h, v), (k, P, h, v))[:k]
        Mb = self._buf(("KS", P, h, v), (1, P, h, v))[0]
        def prod(L_, R_, dst):
            fnp.matmul(L_, fnp.swapaxes(R_, -1, -2), out=T)
            fnp.sum(T, axis=0, out=dst)
''', '''        Lb = self._bsel(("L1", P, h, q), (k, P, h, q), k)
        Rb = self._bsel(("R1", P, v, q), (k, P, v, q), k)
        T = self._bsel(("K", P, h, v), (k, P, h, v), k)
        Mb = self._bsel(("KS", P, h, v), (1, P, h, v), 0)
        def prod(L_, R_, dst):
            fnp.matmul(L_, self._derived(R_, "sw", lambda: fnp.swapaxes(R_, -1, -2)), out=T)
            fnp.sum(T, axis=0, out=dst)
'''))


def main():
    src, dst = sys.argv[1:3]
    t = open(src).read()
    for old, new in REPL:
        assert t.count(old) == 1, old[:120]
        t = t.replace(old, new)
    t = t.replace("# float64 covariance propagation. The suite-shape arithmetic is unchanged.\n",
                  "# float64 covariance propagation. The suite-shape arithmetic is unchanged.\n"
                  "# Second local modification (residual cut only, same ops in the same order): the\n"
                  "# _Strassen helper memoizes the views of its pooled buffers (see `_derived`).\n", 1)
    import os
    os.makedirs(os.path.dirname(dst), exist_ok=True)
    open(dst, "w").write(t)


if __name__ == "__main__":
    main()
