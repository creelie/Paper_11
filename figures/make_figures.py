"""Figures of the paper: three-dimensional schematics rendered as opaque PNG files.

Run from the repository root:  python3 figures/make_figures.py
Reads code/siegel_dyons_out.json for the orbifold indices and writes
figures/fig_indices.png, figures/fig_caustic.png, figures/fig_verdicts.png,
figures/fig_weakcoupling.png, figures/fig_valley.png, figures/fig_rigidity.png.
"""
import json
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401
from PIL import Image

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
BURNT = (204 / 255, 85 / 255, 0.0)
GREY = (0.55, 0.57, 0.62)
YELLOW = (1.0, 241 / 255, 168 / 255)
GREEN = (0.30, 0.55, 0.40)
BLUE = (0.36, 0.52, 0.72)
DPI = 200
W, H = 8.0, 6.0  # inches -> 1600 x 1200 px


def save(fig, name):
    path = os.path.join(HERE, name)
    fig.savefig(path, dpi=DPI, facecolor="white")
    plt.close(fig)
    im = Image.open(path).convert("RGB").quantize(colors=256, method=Image.Quantize.MEDIANCUT)
    im.save(path, optimize=True)
    print(name, Image.open(path).size, Image.open(path).mode)


def fig_indices():
    data = json.load(open(os.path.join(ROOT, "code", "siegel_dyons_out.json")))["tor2_R_no_1_over_N"]
    fig = plt.figure(figsize=(W, H))
    ax = fig.add_subplot(111, projection="3d")
    P2 = [0, 2, 4]
    for q2 in range(0, 13):
        vals = data[f"0,{q2}"]  # Q.P = 0
        for j, p2 in enumerate(P2):
            d = vals[j]
            h = np.sign(d) * np.log10(1 + abs(d))
            col = BURNT if d < 0 else GREY
            z0 = min(0, h)
            ax.bar3d(q2 - 0.35, j - 0.35, z0, 0.7, 0.7, abs(h), color=col, edgecolor="black", linewidth=0.3, shade=True)
    ax.set_xlabel(r"$Q^2$", labelpad=8)
    ax.set_ylabel(r"$P^2$", labelpad=8)
    ax.set_yticks([0, 1, 2])
    ax.set_yticklabels(["0", "2", "4"])
    ax.set_zlabel(r"$\mathrm{sign}(d)\,\log_{10}(1+|d|)$", labelpad=8)
    ax.set_xticks(range(0, 13, 2))
    ax.view_init(elev=24, azim=-58)
    ax.set_box_aspect((2.2, 1.0, 1.1))
    fig.text(0.5, 0.94, r"Full dyon index $d(Q^2,P^2,0)$ of the $\mathbb{Z}_2$ toroidal orbifold", ha="center", fontsize=12)
    fig.text(0.5, 0.05, "burnt orange: negative full index; horizon index is a separate quantity; grey: positive", ha="center", fontsize=10)
    fig.subplots_adjust(left=0.02, right=0.98, bottom=0.06, top=0.95)
    save(fig, "fig_indices.png")


def fig_caustic():
    a, k = 0.5, 1.0
    tc = 1.0 / (a * k * k)
    x0 = np.linspace(0, 4 * np.pi, 1600)
    T0 = a * np.cos(k * x0)
    v0 = a * k * np.sin(k * x0) / np.sqrt(1 + a*a*k*k*np.sin(k*x0)**2)
    v0p = a * k * k * np.cos(k * x0) / (1+a*a*k*k*np.sin(k*x0)**2)**1.5
    ts = np.linspace(0, 0.97 * tc, 60)
    xs = np.linspace(0, 4 * np.pi, 320)
    Z = np.zeros((len(ts), len(xs)))
    for i, t in enumerate(ts):
        X = x0 + t * v0
        dens = 1.0 / (1 + t * v0p)
        order = np.argsort(X)
        Z[i] = np.interp(xs, X[order], dens[order])
    Z = np.minimum(Z, 12)
    TT, XX = np.meshgrid(ts / tc, xs, indexing="ij")
    fig = plt.figure(figsize=(W, H))
    ax = fig.add_subplot(111, projection="3d")
    ax.plot_surface(XX, TT, Z, cmap="copper_r", linewidth=0, antialiased=False, rstride=1, cstride=1, shade=True)
    # characteristics drawn on the floor
    for xi in np.linspace(0, 4 * np.pi, 41):
        tt = np.linspace(0, 1.0, 50)
        xx = xi + tt * tc * (a * k * np.sin(k * xi) / np.sqrt(1+a*a*k*k*np.sin(k*xi)**2))
        ax.plot(xx, tt, np.zeros_like(tt), color="black", linewidth=0.4)
    for m in (0, 1):
        xm = (2 * m + 1) * np.pi / k
        ax.plot([xm, xm], [0.97, 1.0], [0, 12], color=BURNT, linewidth=1.5)
    ax.set_xlabel(r"$x$", labelpad=8)
    ax.set_ylabel(r"$t/t_c$", labelpad=8)
    ax.set_zlabel(r"$\varepsilon(t,x)/\varepsilon_0$", labelpad=8)
    ax.set_xticks([0, np.pi, 2 * np.pi, 3 * np.pi, 4 * np.pi])
    ax.set_xticklabels(["0", r"$\pi/k$", r"$2\pi/k$", r"$3\pi/k$", r"$4\pi/k$"])
    ax.view_init(elev=28, azim=-62)
    ax.set_box_aspect((2.0, 1.2, 1.0))
    fig.text(0.5, 0.94, r"Tachyon dust from the datum $T_0 = T_* + a\cos kx$: density along the characteristics", ha="center", fontsize=12)
    fig.text(0.5, 0.05, r"the characteristics (floor) converge on the minima of $T_0$ and the density diverges at $t_c = 1/(ak^2)$", ha="center", fontsize=10)
    fig.subplots_adjust(left=0.02, right=0.98, bottom=0.06, top=0.95)
    save(fig, "fig_caustic.png")


def fig_verdicts():
    verdict = {1: "c", 2: "p", 3: "c", 4: "p", 5: "p", 6: "p", 7: "p", 8: "p", 9: "i", 10: "c", 11: "p", 12: "p", 13: "p",
               14: "p", 15: "c", 16: "p", 17: "i", 18: "c", 19: "p", 20: "p", 21: "p", 22: "p", 23: "p", 24: "p", 25: "p",
               26: "p", 27: "c", 28: "p", 29: "c", 30: "p", 31: "p", 32: "c", 33: "p", 34: "c", 35: "p", 36: "c"}
    areas = [("Duality", [1, 2, 3, 4, 5]), ("Dyons", [6, 7, 8, 9, 10, 36]), ("Tachyons", list(range(11, 22))),
             ("F-theory", [22, 23]), ("Entropy", list(range(24, 30))), ("SFT, soft, networks", list(range(30, 36)))]
    colour = {"p": GREEN, "c": BLUE, "i": YELLOW, "f": BURNT}
    height = {"p": 1.0, "c": 1.0, "i": 1.0, "f": 1.0}
    fig = plt.figure(figsize=(W, H))
    ax = fig.add_subplot(111, projection="3d")
    labels = []
    for ia, (name, claims) in enumerate(areas):
        for jc, n in enumerate(claims):
            v = verdict[n]
            ax.bar3d(ia - 0.4, jc - 0.4, 0, 0.8, 0.8, height[v], color=colour[v], edgecolor="black", linewidth=0.3, shade=True)
            labels.append((ia, jc, height[v], n))
    ax.set_xticks(range(len(areas)))
    ax.set_xticklabels([a for a, _ in areas], fontsize=8, ha="right", va="center")
    ax.set_yticks(range(0, 11, 2))
    ax.set_ylabel("position in its section", labelpad=10)
    ax.set_zticks([])
    ax.set_zlim(0, 1.6)
    ax.view_init(elev=62, azim=-38)
    ax.set_box_aspect((1.3, 2.0, 0.45))
    from matplotlib.patches import Patch
    ax.legend(handles=[Patch(color=GREEN, label="proved (24)"), Patch(color=BLUE, label="confirmed, unproved (10)"),
                       Patch(color=YELLOW, label="scope limited (2)")],
              loc="upper left", fontsize=8, frameon=False)
    fig.text(0.5, 0.94, "Selected statements: established, supported, scope limited", ha="center", fontsize=12)
    fig.subplots_adjust(left=0.02, right=0.98, bottom=0.04, top=0.95)
    fig.canvas.draw()
    from mpl_toolkits.mplot3d import proj3d
    for ia, jc, h, n in labels:  # numbers drawn in the 2D overlay so that no bar hides them
        x2, y2, _ = proj3d.proj_transform(ia, jc, h, ax.get_proj())
        fx, fy = fig.transFigure.inverted().transform(ax.transData.transform((x2, y2)))
        fig.text(fx, fy, str(n), ha="center", va="center", fontsize=8)
    save(fig, "fig_verdicts.png")


def fig_weakcoupling():
    P = np.polynomial.Polynomial
    x = P([0, 1])
    h, eta, chi = x**4 - 3 * x**2 + 1, x**8 + 2, 3 * x**12 + 5
    xs = np.linspace(-2.3, 2.3, 461)
    les = np.linspace(-3.0, -0.5, 51)
    XX, LE = np.meshgrid(xs, les, indexing="ij")
    E = 10.0 ** LE

    def logj(xv, ev):
        f = -3 * h(xv)**2 + ev * eta(xv)
        g = -2 * h(xv)**3 + ev * h(xv) * eta(xv) - ev**2 * chi(xv) / 12
        return np.log10(1 + np.abs(1728 * 4 * f**3 / (4 * f**3 + 27 * g**2)))
    Z = np.clip(logj(XX, E), 3, 11)  # values below 10^3 are drawn on the floor
    fig = plt.figure(figsize=(W, H))
    ax = fig.add_subplot(111, projection="3d")
    ax.plot_surface(XX, LE, Z, cmap="Blues", linewidth=0, antialiased=False, rstride=1, cstride=2, shade=True, alpha=0.6, vmin=0, vmax=11)
    for r in h.roots():
        ax.plot([r.real] * len(les), les, np.clip(logj(r.real, 10.0 ** les), 3, 11) + 0.08, color=BURNT, linewidth=2.6, zorder=10)
    for r in (eta**2 - h * chi).roots():
        if abs(r.imag) < 1e-9 and abs(r.real) <= 2.3:
            ax.plot([r.real] * len(les), les, np.full(len(les), 11.0), color="black", linewidth=1.4, linestyle="--", zorder=10)
            ax.plot([r.real] * len(les), les, np.full(len(les), 3.0), color="black", linewidth=1.0, linestyle=":", zorder=10)
    ax.set_xlabel(r"base coordinate $x$", labelpad=8)
    ax.set_ylabel(r"$\log_{10}\varepsilon$", labelpad=8)
    ax.set_zlabel(r"$\log_{10}(1+|j|)$", labelpad=6)
    ax.set_zlim(3, 11.5)
    ax.view_init(elev=24, azim=-72)
    ax.set_box_aspect((2.2, 1.2, 1.0))
    fig.text(0.5, 0.94, r"Sen's weak-coupling limit over $\mathbb{P}^1$: $|j|$ grows like $\varepsilon^{-2}$ away from the O7-planes", ha="center", fontsize=12)
    fig.text(0.5, 0.05, r"orange: O7-planes $h=0$, where $j\to1728$; black dashed: D7-branes $\eta^2=h\chi$, where $j=\infty$", ha="center", fontsize=10)
    fig.subplots_adjust(left=0.02, right=0.98, bottom=0.06, top=0.95)
    save(fig, "fig_weakcoupling.png")


def lovelock_GH(D, lam):
    from math import factorial
    n = D - 2

    def G(y):
        return sum(2 * m * lam[m - 1] * factorial(n) / factorial(n + 2 - 2 * m) * y ** (1 - m)
                   for m in range(1, len(lam) + 1) if n + 2 - 2 * m >= 0)

    def H(y):
        return sum(lam[m - 1] * factorial(n) / factorial(n - 2 * m) * y ** (-m)
                   for m in range(1, len(lam) + 1) if n - 2 * m >= 0)
    return G, H


def prester_coefficients(M):
    from math import factorial
    return [4.0 ** (1 - m) / factorial(m) for m in range(1, M + 1)]


def fig_valley():
    G, H = lovelock_GH(6, prester_coefficients(3))
    v1 = np.linspace(0.2, 0.85, 120)
    v2 = np.linspace(1.0, 8.0, 160)
    V1, V2 = np.meshgrid(v1, v2, indexing="ij")
    R = -V1**2 * H(V2) + V1 * G(V2)
    fig = plt.figure(figsize=(W, H_FIG))
    ax = fig.add_subplot(111, projection="3d")
    ax.plot_surface(V1, V2, R, cmap="viridis", linewidth=0, antialiased=False, rstride=2, cstride=2, shade=True, alpha=0.9)
    ax.plot([0.5] * len(v2), v2, np.ones_like(v2), color=BURNT, linewidth=2.4)
    ys = 1.5 + np.sqrt(3.0)
    ax.scatter([0.5], [ys], [1.0], color="black", s=30, depthshade=False)
    ax.text(0.52, ys, 1.25, r"$v_2=\frac{3}{2}+\sqrt{3}$", fontsize=10)
    ax.set_xlabel(r"$v_1$ ($AdS_2$)", labelpad=8)
    ax.set_ylabel(r"$v_2$ ($S^4$)", labelpad=8)
    ax.set_zlabel(r"$R(v_1,v_2)$", labelpad=6)
    ax.view_init(elev=26, azim=-130)
    ax.set_box_aspect((1.3, 1.8, 1.0))
    fig.text(0.5, 0.94, r"Prester's Lovelock action in $D=6$: $R\equiv\alpha'$ along $v_1=\alpha'/2$", ha="center", fontsize=12)
    fig.text(0.5, 0.05, r"every near-horizon solution lies on the orange line, so $S_{\rm BH}=4\pi\sqrt{nw}$; the solution is the saddle point", ha="center", fontsize=10)
    fig.subplots_adjust(left=0.02, right=0.98, bottom=0.06, top=0.95)
    save(fig, "fig_valley.png")


def fig_rigidity():
    Ds = np.arange(4, 17)
    lb = np.linspace(np.log10(0.25), np.log10(4.0), 61)
    ygrid = 10.0 ** np.linspace(-3, 3.5, 4000)
    Z = np.zeros((len(Ds), len(lb)))
    for i, D in enumerate(Ds):
        m = D // 2
        base = prester_coefficients(m)
        for k, l in enumerate(lb):
            lam = base[:-1] + [base[-1] * 10.0 ** l]
            G, Hf = lovelock_GH(int(D), lam)
            Z[i, k] = np.sqrt(np.min(G(ygrid) ** 2 / Hf(ygrid)) / 4)
    DD, LB = np.meshgrid(Ds, lb, indexing="ij")
    fig = plt.figure(figsize=(W, H_FIG))
    ax = fig.add_subplot(111, projection="3d")
    ax.plot_surface(DD, LB, Z, cmap="coolwarm", linewidth=0.2, edgecolor=(0, 0, 0, 0.25), antialiased=False, shade=True, alpha=0.95)
    ax.plot(Ds, np.zeros_like(Ds, dtype=float), np.ones(len(Ds)), color=BURNT, linewidth=2.6)
    ax.set_xlabel(r"dimension $D$", labelpad=8)
    ax.set_ylabel(r"$\log_{10}\beta$", labelpad=8)
    ax.set_zlabel(r"$S_{\rm BH}/4\pi\sqrt{nw}$", labelpad=8)
    ax.view_init(elev=22, azim=-55)
    ax.set_box_aspect((1.8, 1.4, 1.0))
    fig.text(0.5, 0.94, r"Top Lovelock coefficient $\lambda_{\lfloor D/2\rfloor}=\beta\,4^{1-m}/m!$: entropy at the minimum of $\Phi_\beta$", ha="center", fontsize=12)
    fig.text(0.5, 0.05, r"the statistical value $1$ is reached only at $\beta=1$ (orange), in every dimension", ha="center", fontsize=10)
    fig.subplots_adjust(left=0.02, right=0.98, bottom=0.06, top=0.95)
    save(fig, "fig_rigidity.png")


H_FIG = H

if __name__ == "__main__":
    fig_indices()
    fig_caustic()
    fig_verdicts()
    fig_weakcoupling()
    fig_valley()
    fig_rigidity()
