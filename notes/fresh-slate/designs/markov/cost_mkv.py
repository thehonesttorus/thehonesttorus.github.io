"""Projected flopscope cost of MKV(w) at n=1024, L=16 in units (1 unit = 2^31 FLOPs).
Product price: 0.56 u (V29 Strassen-Winograd L5, notes/streams/costmodel/cost.py); n^2 elementwise op 0.0005 u
(transcendental x16).  Counts follow mkv.py with the grouping noted in DESIGN.md §3/§9."""
import sys
MM, EW = 0.56, 2 ** 20 / 2 ** 31
def units(w, two_site=True, L=16):
    tot = 0.0
    for l in range(L - 1):                       # transport l -> l+1
        S = min(l, w)                            # active sources carrying content into layer l
        last = (l == L - 2)
        field = (0.0 if last else 2 * MM) + 70 * EW        # W^T Cov(a) W (diag-only at the end) + Mehler/ReLU
        newsrc = 1 * MM                          # K = C D W
        per = 2 + 3 + (0 if last else 2) + (2 if two_site else 0)   # P,K; PPU,PKU,KKU; M; two-site
        tot += field + newsrc + S * per * MM + S * 40 * EW + 1 * MM + 200 * EW
    return tot
if __name__ == "__main__":
    for w in [1, 2, 3, 4, 8, 16]:
        print(f"w={w:2d}: {units(w):6.1f} u   (no two-site kappa4: {units(w, False):6.1f} u)")
