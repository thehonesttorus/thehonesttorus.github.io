"""Chaos spectrum of the depth-L output of a wide bias-free He ReLU MLP, and the ceiling it puts on exact-mean controls.

For x ~ N(0, I_d) and two inputs at correlation rho, the per-neuron output kernel is q f^(L)(rho) with
f(rho) = (sqrt(1 - rho^2) + (pi - arccos rho) rho) / pi (degree-1 arc-cosine map, q = E h^2). Its power series
f^(L)(rho) = sum_k b_k rho^k has b_k >= 0, sum b_k = 1, and b_k is the share of E h^2 in Wiener chaos k (for the
sphere, degree-k harmonics; the two agree as d -> inf). Hence the mean is b_0 (|mu|^2 / n = q b_0), the per-neuron
variance is q (1 - b_0), and any control built from input polynomials of degree <= k (all of which have exactly known
means) removes at most sum_(1 <= j <= k) b_j / (1 - b_0) of the variance. f maps the unit disc into itself, so f^(L)
is analytic there and b_k follows from a Cauchy integral on |rho| = r < 1 (FFT).

  python notes/stage24/theory24/chaos_spectrum.py"""
import numpy as np


def f(r):
    return (np.sqrt(1 - r * r) + (np.pi - np.arccos(r)) * r) / np.pi


def spectrum(L=16, r=0.97, N=1 << 15, K=400):
    phi = 2 * np.pi * np.arange(N) / N; z = r * np.exp(1j * phi); F = z.copy()
    for _ in range(L): F = f(F)
    b = (np.fft.fft(F) / N).real[:K] / r ** np.arange(K)
    return b


if __name__ == "__main__":
    for L in (1, 2, 4, 8, 12, 16):
        b = spectrum(L); v = 1 - b[0]; c = np.cumsum(b[1:]) / v
        x = 0.0
        for _ in range(L): x = f(x)
        print(f"L={L:2d}: b0 = {b[0]:.6f} (direct f^L(0) = {x:.6f}), variance 1-b0 = {v:.5f}; sum_k<400 b_k = {b.sum():.6f};"
              f" variance share in chaos 1, 2, 3, 4: {b[1] / v:.4f} {b[2] / v:.4f} {b[3] / v:.4f} {b[4] / v:.4f};"
              f" cumulative <=1/2/4/8/16/32/64: " + " ".join(f"{c[k - 1]:.3f}" for k in (1, 2, 4, 8, 16, 32, 64)))
    b = spectrum(16); k = np.arange(len(b))
    sl = np.polyfit(np.log(k[50:300]), np.log(b[50:300]), 1)[0]
    print(f"L=16 tail: b_k ~ k^{sl:.2f} over 50 <= k < 300 (f(1 - e) = 1 - e + (sqrt2/pi) e^(3/2) + ..., and composition keeps the e^(3/2) kink: b_k ~ k^(-5/2) asymptotically);"
          f" max variance reduction factor from exact-mean controls of degree <= k: " +
          " ".join(f"k={kk}: {1 / (1 - np.sum(b[1:kk + 1]) / (1 - b[0])):.2f}x" for kk in (1, 2, 3, 4, 8)))
