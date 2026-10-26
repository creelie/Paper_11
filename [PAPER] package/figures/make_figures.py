"""Figures of the paper: three-dimensional schematics rendered as opaque PNG files.

Run from the repository root:  python3 figures/make_figures.py
Reads code/siegel_dyons_out.json for the orbifold indices and writes
figures/fig_indices.png, figures/fig_caustic.png, figures/fig_verdicts.png.
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
               14: "p", 15: "c", 16: "p", 17: "i", 18: "c", 19: "p", 20: "p", 21: "p", 22: "p", 23: "i", 24: "p", 25: "p",
               26: "i", 27: "c", 28: "p", 29: "c", 30: "p", 31: "p", 32: "c", 33: "p", 34: "c", 35: "p", 36: "c"}
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
    ax.legend(handles=[Patch(color=GREEN, label="proved (22)"), Patch(color=BLUE, label="confirmed, unproved (10)"),
                       Patch(color=YELLOW, label="scope limited (4)")],
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


if __name__ == "__main__":
    fig_indices()
    fig_caustic()
    fig_verdicts()
