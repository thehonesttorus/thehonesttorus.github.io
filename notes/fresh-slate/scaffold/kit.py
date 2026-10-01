"""numpy -> flopscope porting kit: helpers for the primitives fresh-slate designs are likely to
need, each written in the form that is cheapest in units AND in flopscope calls. Measured costs at
n = 1024 are in KIT.md (probe_kit.py regenerates them). Copy what you need into your estimator.py:
the grader ships one file unless you point `whest package` at a folder, and it has no numpy.

Conventions: 1 unit = 2^31 FLOPs = one dense (1024 x 1024) @ (1024 x 1024) float32 product;
B = 1024 units. "calls" = flopscope ops (each ~0.02-0.04 ms of residual, see RESIDUAL.md).
"""
import flopscope as flops
import flopscope.numpy as fnp

F32 = fnp.float32


# ---- products -------------------------------------------------------------------------------
def mm(a, b, out=None):
    """Dense product, 2mkn - mn FLOPs, 1 call. Batched (k, m, n) @ (k, n, p) is one call at k
    times the price (no batching discount). Any float64 operand bills the whole op at 2x."""
    return fnp.matmul(a, b, out=out) if out is not None else fnp.matmul(a, b)


def gram_same(x):
    """x^T x as an aliased 2-operand einsum: HALF price (0.5 u at n = 1024), 1 call, output tagged
    symmetric (legitimately: it is a Gram of one object). x.T @ x (a view = another object) bills
    full price. Returns a fresh tagged array: copy it into an untagged pooled buffer before any
    out= writes of plain results into it."""
    return fnp.einsum("ji,jk->ik", x, x)


def gram_weighted_pos(x, d, tmp):
    """x^T diag(d) x for d >= 0 (elementwise): sqrt(d) * x into tmp, then the aliased Gram.
    0.5 u + O(n^2). For sign-indefinite d use two Grams over the rows with d > 0 and d <= 0."""
    fnp.multiply(x, fnp.sqrt(d)[:, None], out=tmp)
    return fnp.einsum("ji,jk->ik", tmp, tmp)


def sandwich(w, c, tmp, out):
    """w^T c w for any c (no symmetry tags involved): 2 dense products, 2.0 u, 2 calls.
    The tagged einsum('ij,ia,jb->ab', as_symmetric(c), w, w) bills 1.5 u but raises a float32
    SymmetryError on indefinite O(1) inputs (zeroed MLP): do not use it on anything but SPD
    inputs, and even then guard it."""
    fnp.matmul(c, w, out=tmp)
    fnp.matmul(w.T, tmp, out=out)
    return out


def diag_sandwich(w, c, tmp, out_vec):
    """diag(w^T c w) only: c @ w (1 u) then column sums of w * (c w): 1.0 u, 3 calls."""
    fnp.matmul(c, w, out=tmp)
    fnp.multiply(tmp, w, out=tmp)
    return fnp.sum(tmp, axis=0, out=out_vec)


# ---- Gaussian elementwise -------------------------------------------------------------------
def gauss_cdf_pdf(z):
    """Phi(z), phi(z) for a float32 vector/matrix. flopscope.stats.norm always bills float64
    (cdf 96, pdf 54 FLOPs per element as billed); cast back to float32. On (n,) vectors this is
    negligible (5e-5 u); on (n, n) matrices it is 0.047 + 0.026 u: prefer vector-level gates."""
    return (flops.stats.norm.cdf(z).astype(F32), flops.stats.norm.pdf(z).astype(F32))


def relu_moments(m, v):
    """E relu(Z), E relu(Z)^2, P(Z > 0) for Z ~ N(m, v) per element (vectors)."""
    s = fnp.sqrt(fnp.maximum(v, 1e-30))
    a = m / s
    Phi, phi = gauss_cdf_pdf(a)
    mean = m * Phi + s * phi
    second = (m * m + v) * Phi + m * s * phi
    return mean, second, Phi


# ---- linear algebra (see KIT.md for measured prices) ----------------------------------------
def chol(a):
    return fnp.linalg.cholesky(a)


def eigh(a):
    return fnp.linalg.eigh(a)


def slogdet(a):
    """Batched (k, m, m) works; n^3-class at n = 1024 (see KIT.md). det via LU internally."""
    return fnp.linalg.slogdet(a)


def solve(a, b):
    return fnp.linalg.solve(a, b)


def qr(a):
    return fnp.linalg.qr(a)


# ---- gathers --------------------------------------------------------------------------------
def gather_rows(x, idx):
    """x[idx] as fnp.take: billed 4 FLOPs per gathered element, 1 call. Basic slices are free."""
    return fnp.take(x, idx, axis=0)


# ---- pooled buffers -------------------------------------------------------------------------
class Pool:
    def __init__(self):
        self.bufs = {}

    def get(self, name, shape, dtype=None):
        dtype = dtype or F32
        key = (name, tuple(int(x) for x in shape), str(dtype))
        b = self.bufs.get(key)
        if b is None:
            b = fnp.empty(key[1], dtype=dtype)
            fnp.copyto(b, 0.0)
            self.bufs[key] = b
        return b


# ---- Strassen-Winograd families (lower-op engine from notes/streams/est-cost/sw2.py) ---------
# Usage:  S = Strassen(F32, bmax);  S.mm(X, Y, out, S.level(m, kd, w, L))
#   X (bx, P, m, kd) @ Y (by, P, kd, w) -> out (by, P, m, w), bx in {1, by}; for one product use
#   X[None, None] etc. hub: X (k, P, m, kd), Y (k, P, w, kd) -> out (P, m, w) = sum_k X_k Y_k^T.
# Price per (1024^3) product at L0..L5: 1.0 / 0.877 / 0.772 / 0.683 / 0.609 / 0.556 units; calls per
# family are independent of the batch: L5 79 (hub 87), L3 46; f32 relative error 8e-6 at L5.
# Pass cached grids (S.gridc(buf, buf4)) for pooled operands, else each call rebuilds them
# (+15 calls and a reshape billed per element).
STRASSEN_MIN = 32
STRASSEN_FUSE_P = 343


class _G:
    """Quadrant grid of a (b, P, m, k) operand: g = (b, P, 2, h, 2, q) and derived pair views,
    all with the batch axis first so that sub(b0, b1) slices them consistently."""
    __slots__ = ("g", "d", "a", "ac", "t", "b")

    def __init__(self, g, d=None, a=None, ac=None, t=None):
        self.g = g
        self.d = d     # [G00, G11]  (b, P, 2, h, q)
        self.a = a     # [G10, G01]
        self.ac = ac   # [G01, G10]
        self.t = t     # row-pair view (b, P, 2(c), h, 2(r), q)
        self.b = g.shape[0]

    @staticmethod
    def make(x, hubT=False):
        """x (b, P, m, k) -> grid. hubT: block-transposed grid (H[r][c] = x block (c, r))."""
        b, P, m, k = x.shape
        g = fnp.reshape(x, (b, P, 2, m // 2, 2, k // 2))
        if hubT:
            g = fnp.transpose(g, (0, 1, 4, 3, 2, 5))
        return _G(g)

    def full(self, need_t=False):
        g = self.g
        if self.d is None:
            self.d = fnp.moveaxis(fnp.diagonal(g, axis1=2, axis2=4), -1, 2)
            self.a = fnp.moveaxis(fnp.diagonal(g[:, :, ::-1], axis1=2, axis2=4), -1, 2)
            self.ac = fnp.moveaxis(fnp.diagonal(g[:, :, :, :, ::-1], axis1=2, axis2=4), -1, 2)
        if need_t and self.t is None:
            self.t = fnp.transpose(g, (0, 1, 4, 3, 2, 5))
        return self

    def sub(self, b0, b1):
        if b0 == 0 and b1 == self.b:
            return self
        s = slice(b0, b1)
        return _G(self.g[s], None if self.d is None else self.d[s], None if self.a is None else self.a[s],
                  None if self.ac is None else self.ac[s], None if self.t is None else self.t[s])


class Strassen:
    """plain: X (bx, P, m, kd) @ Y (by, P, kd, w) -> out (by, P, m, w), bx in {1, by}
       hub:   X (k, P, m, kd), Y (k, P, w, kd) -> out (P, m, w) = sum_k X_k Y_k^T
    Operands may be passed as arrays or as prebuilt grids (_G) of pooled buffers."""

    def __init__(self, dtype, bmax):
        self.dtype = dtype
        self.bmax = int(bmax)
        self.pool = {}
        self.grids = {}

    # ---- pools -------------------------------------------------------------------------
    def _alloc(self, shape):
        b = fnp.empty(shape, dtype=self.dtype)
        fnp.copyto(b, 0.0)   # metered page touch
        return b

    def _cb(self, key, b, P, h, q, hubT=False):
        """Combos buffer (b, P, 7, h, q) + its next-level operand view (b, 7P, h, q) and the
        grid of that view (built once per allocation)."""
        e = self.pool.get(key)
        if e is None or e[0].shape[0] < b:
            buf = self._alloc((b, P, 7, h, q))
            v4 = fnp.reshape(buf, (b, 7 * P, h, q))
            e = [buf, v4, {}]
            self.pool[key] = e
        return e

    def _grid_of(self, e, b, hubT=False, recurse=True, need_t=False):
        """Grid of entry e's next-level operand view, batch [:b] (cached per (hubT))."""
        gk = ("G", hubT)
        G = e[2].get(gk)
        if G is None:
            v4 = e[1]
            G = _G.make(v4, hubT=hubT)
            e[2][gk] = G
        if recurse:
            G.full(need_t)
        elif need_t:
            G.full(True)
        return G.sub(0, b)

    def _flat(self, key, shape):
        b = self.pool.get(key)
        if b is None or b.shape[0] < shape[0]:
            b = self._alloc(shape)
            self.pool[key] = b
        return b[:shape[0]]

    def gridc(self, base, view4, hubT=False):
        """Cached grid of a persistent pooled operand (identity key `base`, operand view `view4`
        of shape (b, P, m, k)); slice it with .sub(b0, b1). Pass as XG/YG/OG to mm/hub so the
        reshape/diagonal views (logged calls, reshape billed per element) are built once."""
        ent = self.grids.get((id(base), hubT))
        if ent is None or ent[0] is not base:
            ent = (base, _G.make(view4, hubT=hubT).full(True))
            self.grids[(id(base), hubT)] = ent
        return ent[1]

    def level(self, m, kd, w, lev):
        while lev > 0 and not self._ok(m, kd, w, lev):
            lev -= 1
        return lev

    @staticmethod
    def _ok(m, kd, w, lev):
        d = 2 ** lev
        return (lev > 0 and m % d == 0 and kd % d == 0 and w % d == 0
                and min(m, kd, w) // d >= STRASSEN_MIN)

    # ---- one level: combos --------------------------------------------------------------
    @staticmethod
    def _lcombos(G, C):
        """G grid of X, C (b, P, 7, h, q) <- [M1, M2, M5, M3, M4, M6, M7] left combos."""
        g = G.g
        fnp.add(g[:, :, :, :, 0], g[:, :, 1:2, :, 1], out=C[:, :, 0:2])   # [X11+X22, X21+X22]
        fnp.add(g[:, :, 0, :, 0], g[:, :, 0, :, 1], out=C[:, :, 2])        # X11+X12
        fnp.copyto(C[:, :, 3:5], G.d)                                      # [X11, X22]
        fnp.subtract(G.a, G.d, out=C[:, :, 5:7])                           # [X21-X11, X12-X22]

    @staticmethod
    def _rcombos(H, C):
        """H grid of Y (or block-transposed Y for the hub) -> right combos, same slot order."""
        h = H.g
        fnp.add(h[:, :, 0, :, 0], h[:, :, 1, :, 1], out=C[:, :, 0])        # Y11+Y22
        fnp.copyto(C[:, :, 1:3], H.d)                                      # [Y11, Y22]
        fnp.subtract(H.ac, H.d[:, :, ::-1], out=C[:, :, 3:5])              # [Y12-Y22, Y21-Y11]
        fnp.add(h[:, :, :, :, 0], h[:, :, :, :, 1], out=C[:, :, 5:7])      # [Y11+Y12, Y21+Y22]

    @staticmethod
    def _assemble(M, o):
        """M (b, P, 7, h, w) products -> out grid o (b, P, 2, h, 2, w)."""
        C11, C12, C21, C22 = o[:, :, 0, :, 0], o[:, :, 0, :, 1], o[:, :, 1, :, 0], o[:, :, 1, :, 1]
        fnp.add(M[:, :, 0:2], M[:, :, 4:5], out=o[:, :, :, :, 0])   # [C11, C21] = [M1+M4, M2+M4]
        fnp.add(M[:, :, 3], M[:, :, 2], out=C12)                    # C12 = M3+M5
        fnp.subtract(M[:, :, 0], M[:, :, 1], out=C22)               # C22 = M1-M2
        fnp.add(C22, M[:, :, 3], out=C22)                           #     + M3
        fnp.add(C22, M[:, :, 5], out=C22)                           #     + M6
        fnp.subtract(C11, M[:, :, 2], out=C11)                      # C11 -= M5
        fnp.add(C11, M[:, :, 6], out=C11)                           #     += M7

    # ---- plain product ------------------------------------------------------------------
    def mm(self, X, Y, out, lev, XG=None, YG=None, OG=None):
        bx, P, m, kd = X.shape
        by, w = Y.shape[0], Y.shape[3]
        if not self._ok(m, kd, w, lev):
            fnp.matmul(X, Y, out=out)
            return
        h, q, v = m // 2, kd // 2, w // 2
        XG = (XG or _G.make(X)).full()
        YG = (YG or _G.make(Y)).full()
        OG = OG or _G.make(out)
        if P >= STRASSEN_FUSE_P and not self._ok(h, q, v, lev - 1):
            self._leaf_mm(XG, YG, OG, bx, by, P, h, q, v)
            return
        # the products reuse this level's left-combo buffer when the callee recurses (its
        # combos consume X before its assemble writes M), exactly as V29
        alias = self._ok(h // 2, q // 2, v // 2, lev - 2) and q == v
        eL = self._cb(("L", P, h, q), max(bx, by) if alias else bx, P, h, q)
        eR = self._cb(("R", P, q, v), by, P, q, v)
        self._lcombos(XG, eL[0][:bx])
        self._rcombos(YG, eR[0][:by])
        Mkey_e = eL if alias else self._cb(("M", P, h, v), by, P, h, v)
        Xn = eL[1][:bx]
        Yn = eR[1][:by]
        Mn = Mkey_e[1][:by]
        rec = self._ok(h, q, v, lev - 1)
        self.mm(Xn, Yn, Mn, lev - 1,
                XG=self._grid_of(eL, bx, recurse=rec) if rec else None,
                YG=self._grid_of(eR, by, recurse=rec) if rec else None,
                OG=self._grid_of(Mkey_e, by, recurse=False) if rec else None)
        self._assemble(Mkey_e[0][:by], OG.g)

    def _pairbuf(self, key, b, P, h, q):
        return self._flat(key, (b, P, 2, h, q))

    def _leaf_mm(self, XG, YG, OG, bx, by, P, h, q, v):
        """Fused leaf: 7 products in 4 matmuls (pairs + M1), 19 ops."""
        g, hh, o = XG.g, YG.g, OG.g
        C11, C12, C21, C22 = o[:, :, 0, :, 0], o[:, :, 0, :, 1], o[:, :, 1, :, 0], o[:, :, 1, :, 1]
        if OG.t is None:
            OG.t = fnp.transpose(o, (0, 1, 4, 3, 2, 5))
        row1r = OG.t[:, :, ::-1, :, 1]                       # [C22, C21]
        Lb = self._pairbuf(("PL", P, h, q), bx, P, h, q)
        Rb = self._pairbuf(("PR", P, q, v), by, P, q, v)
        Mb = self._pairbuf(("PM", P, h, v), by, P, h, v)
        # (M3, M4) = [X11, X22] @ [Y12-Y22, Y21-Y11] -> [C22, C21]
        fnp.subtract(YG.ac, YG.d[:, :, ::-1], out=Rb)
        fnp.matmul(XG.d, Rb, out=row1r)
        fnp.copyto(C12, C22)                                 # C12 = M3
        # M1 = (X11+X22)(Y11+Y22)
        fnp.add(g[:, :, 0, :, 0], g[:, :, 1, :, 1], out=Lb[:, :, 0])
        fnp.add(hh[:, :, 0, :, 0], hh[:, :, 1, :, 1], out=Rb[:, :, 0])
        fnp.matmul(Lb[:, :, 0], Rb[:, :, 0], out=Mb[:, :, 0])
        fnp.add(Mb[:, :, 0], C21, out=C11)                   # C11 = M1 + M4
        fnp.add(C22, Mb[:, :, 0], out=C22)                   # C22 = M3 + M1
        # (M2, M5) = [X21+X22, X11+X12] @ [Y11, Y22]
        fnp.add(g[:, :, ::-1, :, 0], g[:, :, ::-1, :, 1], out=Lb)
        fnp.matmul(Lb, YG.d, out=Mb)
        fnp.add(C21, Mb[:, :, 0], out=C21)                   # C21 = M4 + M2
        fnp.subtract(C22, Mb[:, :, 0], out=C22)              # C22 -= M2
        fnp.add(C12, Mb[:, :, 1], out=C12)                   # C12 += M5
        fnp.subtract(C11, Mb[:, :, 1], out=C11)              # C11 -= M5
        # (M6, M7) = [X21-X11, X12-X22] @ [Y11+Y12, Y21+Y22]
        fnp.subtract(XG.a, XG.d, out=Lb)
        fnp.add(hh[:, :, :, :, 0], hh[:, :, :, :, 1], out=Rb)
        fnp.matmul(Lb, Rb, out=Mb)
        fnp.add(C22, Mb[:, :, 0], out=C22)                   # C22 += M6
        fnp.add(C11, Mb[:, :, 1], out=C11)                   # C11 += M7

    # ---- hub ----------------------------------------------------------------------------
    def hub(self, X, Y, out, lev, XG=None, YG=None, OG=None):
        """X (k, P, m, kd), Y (k, P, w, kd) -> out (P, m, w); out passed as (1, P, m, w)."""
        k, P, m, kd = X.shape
        w = Y.shape[2]
        if not self._ok(m, kd, w, lev):
            T = self._flat(("K", P, m, w), (k, P, m, w))
            fnp.matmul(X, fnp.swapaxes(Y, -1, -2), out=T)
            fnp.sum(T, axis=0, out=out[0])
            return
        h, q, v = m // 2, kd // 2, w // 2
        XG = (XG or _G.make(X)).full()
        YG = (YG or _G.make(Y, hubT=True)).full()
        OG = OG or _G.make(out)
        if P >= STRASSEN_FUSE_P and not self._ok(h, q, v, lev - 1):
            self._leaf_hub(XG, YG, OG, k, P, h, q, v)
            return
        eL = self._cb(("L", P, h, q), k, P, h, q)
        eR = self._cb(("R", P, v, q), k, P, v, q)
        self._lcombos(XG, eL[0][:k])
        self._rcombos(YG, eR[0][:k])
        eM = self._cb(("HM", P, h, v), 1, P, h, v)
        rec = self._ok(h, q, v, lev - 1)
        self.hub(eL[1][:k], eR[1][:k], eM[1][:1], lev - 1,
                 XG=self._grid_of(eL, k) if rec else None,
                 YG=self._grid_of(eR, k, hubT=True) if rec else None,
                 OG=self._grid_of(eM, 1, recurse=False) if rec else None)
        self._assemble(eM[0][:1], OG.g)

    def _leaf_hub(self, XG, YG, OG, k, P, h, q, v):
        """Fused hub leaf: like _leaf_mm, each pair product = batched GEMM over (k, P, 2)
        into T, then a k-sum into its destination."""
        g, hh, o = XG.g, YG.g, OG.g
        C11, C12, C21, C22 = o[0, :, 0, :, 0], o[0, :, 0, :, 1], o[0, :, 1, :, 0], o[0, :, 1, :, 1]
        if OG.t is None:
            OG.t = fnp.transpose(o, (0, 1, 4, 3, 2, 5))
        row1r = OG.t[0, :, ::-1, :, 1]                       # [C22, C21]  (P, 2, h, v)
        Lb = self._pairbuf(("PL", P, h, q), k, P, h, q)
        Rb = self._pairbuf(("PRh", P, v, q), k, P, v, q)
        T = self._pairbuf(("PT", P, h, v), k, P, h, v)
        Mb = self._pairbuf(("PMh", P, h, v), 1, P, h, v)[0]

        def prod(Lx, Ry, dst, pair=True):
            Tt = T if pair else T[:, :, 0]
            fnp.matmul(Lx, fnp.swapaxes(Ry, -1, -2), out=Tt)
            fnp.sum(Tt, axis=0, out=dst)

        fnp.subtract(YG.ac, YG.d[:, :, ::-1], out=Rb)
        prod(XG.d, Rb, row1r)                                 # [C22, C21] = (M3, M4)
        fnp.copyto(C12, C22)
        fnp.add(g[:, :, 0, :, 0], g[:, :, 1, :, 1], out=Lb[:, :, 0])
        fnp.add(hh[:, :, 0, :, 0], hh[:, :, 1, :, 1], out=Rb[:, :, 0])
        prod(Lb[:, :, 0], Rb[:, :, 0], Mb[:, 0], pair=False)  # M1
        fnp.add(Mb[:, 0], C21, out=C11)
        fnp.add(C22, Mb[:, 0], out=C22)
        fnp.add(g[:, :, ::-1, :, 0], g[:, :, ::-1, :, 1], out=Lb)
        prod(Lb, YG.d, Mb)                                    # (M2, M5)
        fnp.add(C21, Mb[:, 0], out=C21)
        fnp.subtract(C22, Mb[:, 0], out=C22)
        fnp.add(C12, Mb[:, 1], out=C12)
        fnp.subtract(C11, Mb[:, 1], out=C11)
        fnp.subtract(XG.a, XG.d, out=Lb)
        fnp.add(hh[:, :, :, :, 0], hh[:, :, :, :, 1], out=Rb)
        prod(Lb, Rb, Mb)                                      # (M6, M7)
        fnp.add(C22, Mb[:, 0], out=C22)
        fnp.add(C11, Mb[:, 1], out=C11)
