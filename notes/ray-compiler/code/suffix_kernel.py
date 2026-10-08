# A priori suffix-risk weights of the He ReLU suffix (infinite width; finite width moves them by O(D/n), about 1.5% here).
#   python suffix_kernel.py
# K_D = rho^(oD), rho(c) = (sqrt(1 - c^2) + (pi - arccos c) c) / pi, is the normalised two-input kernel of D random He
# layers. For a perturbation of the law at post-activation layer k, the expected output energy over random suffix weights
# is the double integral of k_D against the perturbation; pairs of independent samples sit at the typical cosine c_k
# (c_0 = rho(0) for x ~ N(0, I), c_(k+1) = rho(c_k)), and the j-th Taylor coefficient K_D^(j)(c_k) / j! weighs the part of
# the perturbation of order j in (c - c_k). Mean errors read K', second-order (covariance) errors K'', and so on.
# Jets: K_D(c_k + e) as a truncated power series in e, composed layer by layer.
import math
import numpy as np

J = 5


def rho_jet(c):                      # Taylor coefficients of rho at c, orders 0..J-1
    s = math.sqrt(1 - c * c)
    d = [(s + (math.pi - math.acos(c)) * c) / math.pi, 0.5 + math.asin(c) / math.pi, 1 / (math.pi * s),
         c / (math.pi * s ** 3), (1 + 2 * c * c) / (math.pi * s ** 5)]
    return [d[j] / math.factorial(j) for j in range(J)]


def compose(outer, inner):           # outer(inner(e)) for power series with inner[0] the expansion point of outer
    res = np.zeros(J); res[0] = outer[0]
    dlt = np.array(inner, float); dlt[0] = 0.0
    pw = np.zeros(J); pw[0] = 1.0
    for j in range(1, J):
        pw = np.convolve(pw, dlt)[:J]
        res += outer[j] * pw
    return res


L = 16
c = [1 / math.pi]
for _ in range(L - 1):
    c.append(rho_jet(c[-1])[0])
rows = []
for k in range(L - 1):               # error at post-activation layer k, D = L - 1 - k layers to the output
    D = L - 1 - k
    jet = np.array([c[k], 1.0, 0, 0, 0])
    for i in range(D):
        jet = compose(rho_jet(jet[0]), jet)
    rows.append((k, D, c[k], jet))
out = ["post-activation layer k, layers to the output D, typical pair cosine c_k, and the Taylor coefficients of K_D at c_k",
       "  k  D   c_k     K_D     K'      K''/2    K'''/6   K''''/24"]
for k, D, ck, jet in rows:
    out.append(f" {k:2d} {D:2d}  {ck:.4f}  {jet[0]:.4f}  {jet[1]:.4f}  {jet[2]:8.4f} {jet[3]:8.3f} {jet[4]:9.2f}")
out.append("Reading: K' is the mean-channel weight (the angular contraction prod rho'(c_j) along the suffix); the higher\n"
           "coefficients grow toward c -> 1 because rho has the (1 - c)^(3/2) branch point there.")
print("\n".join(out))
open(__file__.replace("code/suffix_kernel.py", "outputs/suffix_kernel.txt"), "w").write("\n".join(out) + "\n")
