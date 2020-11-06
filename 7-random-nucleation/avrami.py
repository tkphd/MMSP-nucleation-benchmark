# Plotting script for PFHub Nucleation Benchmark
# Questions/comments to trevor.keller@nist.gov (Trevor Keller)

# Y = 1 - exp(-K*t**n)
# y0 = t**n*exp(-K*t**n)
# y1 = K*t**n*exp(-K*t**n)*log(t)

from math import floor
from matplotlib import style
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.optimize import curve_fit
from scipy.stats import chisquare, describe
from string import ascii_letters as letters

style.use("seaborn")

title = "PFHub Benchmark 8.3 Pre-seed"
tlim = [0, 600]
p0 = (5.0e-8, 3.0) # initial guess for non-linear solver

frames = [
    #pd.read_csv("run-a/free_energy.csv"),
    #pd.read_csv("run-b/free_energy.csv"),
    #pd.read_csv("run-c/free_energy.csv"),
    #pd.read_csv("run-d/free_energy.csv"),
    pd.read_csv("run-e/free_energy.csv"),
    pd.read_csv("run-f/free_energy.csv"),
    #pd.read_csv("run-g/free_energy.csv"),
    #pd.read_csv("run-h/free_energy.csv"),
    #pd.read_csv("run-i/free_energy.csv"),
    #pd.read_csv("run-j/free_energy.csv"),
    #pd.read_csv("run-k/free_energy.csv"),
    #pd.read_csv("run-l/free_energy.csv"),
    #pd.read_csv("run-m/free_energy.csv"),
    #pd.read_csv("run-n/free_energy.csv"),
    #pd.read_csv("run-o/free_energy.csv"),
    #pd.read_csv("run-p/free_energy.csv"),
    #pd.read_csv("run-q/free_energy.csv"),
    #pd.read_csv("run-r/free_energy.csv"),
    #pd.read_csv("run-s/free_energy.csv"),
    #pd.read_csv("run-t/free_energy.csv"),
]

figsize = (10, 6)

# === Equations ===

def f_jmak(t, K, n):
    # JMAK growth law, Y(t) = 1 - exp(-Ktⁿ)
    # where $n$ is the spatial dimension
    return 1.0 - np.exp(-K * t**n)

def df_jmak(t, K, n):
    # Jacobian: df/dp for p=(K, n)
    return np.array([
        t**n * np.exp(-K * t**n),
        K * t**n * np.exp(-K * t**n) * np.log(t),
        # -K * n * t**n * np.exp(-K * t**n) / (t - t0),
    ]).T

def jmak_x(x):
    return np.log(x)

def jmak_y(y):
    return np.log(-np.log(1 - y))

def sigfig(x, n):
    # Round a float, x, to n significant figures.
    # Source: https://github.com/corriander/python-sigfig
    n = int(n)

    e = np.floor(np.log10(np.abs(x)) - n + 1)  # exponent, 10 ** e
    shifted_dp = x / (10 ** e)  # decimal place shifted n d.p.
    return np.around(shifted_dp) * (10 ** e)  # round and revert

# === Avrami/JMAK Plots ===

plt.figure(figsize=figsize)
plt.title(title)
plt.xlabel("$\\log(t)$")
plt.ylabel("$\\log(-\\log(1-Y))$")

for i, df in enumerate(frames):
    df = df[df["time"] > 0]
    df = df[df["fraction"] > 0]
    df = df[df["fraction"] < 1]

    t = np.array(df["time"])
    y = np.array(df["fraction"])

    run = "Run {}".format(letters[i + 26])
    plt.plot(jmak_x(t), jmak_y(y), label=None)

# === Levenburg-Marquardt Least-Squares Fit ===

fit_t = np.array([])
fit_y = np.array([])

K = []
n = []
t0 = []

for i, df in enumerate(frames):
    df = df[df["time"] > 0]
    df = df[df["fraction"] > 0]
    df = df[df["fraction"] < 1]

    t = np.array(df["time"])
    y = np.array(df["fraction"])

    # Fit this dataset & print coeffs
    p, pcov = curve_fit(f_jmak, t, y, p0=p0, sigma=None,
                        method="lm", jac=df_jmak, maxfev=2000)
    print("Run", letters[i + 26], "coeffs:", p)

    K.append(p[0])
    n.append(p[1])

    fit_t = np.append(fit_t, t)
    fit_y = np.append(fit_y, y)

print()
p_naive = np.array([np.average(K), np.average(n)])
p_nstd = np.array([np.std(K), np.std(n)])
print("Individual fit: K={0:5.3e} n={1:5.3f}".format(*p_naive))
print("         stdev:   {0:5.3e}   {1:5.3f}".format(*p_nstd))

p, pcov = curve_fit(f_jmak, fit_t, fit_y, p0=p0, sigma=None,
                    method="lm", jac=df_jmak, maxfev=2000)
perr = np.sqrt(np.diag(pcov))

print()
print("Collective fit: K={0:5.3e} n={1:5.3f}".format(*p))
print("         error:   {0:5.3e}   {1:5.3f}".format(*perr))

fit_max = np.amax(fit_t)
fit_min = np.exp(floor(np.log(fit_max) / 3))

t_hat = np.linspace(fit_max, fit_min, 201)
y_hat = f_jmak(t_hat, *p)

jx = jmak_x(t_hat)
jy = jmak_y(y_hat)
eqn = "$1-\\exp\{-(%.2f \\pm %.2f) \\times 10^{-9} \\times t^{%.2f \\pm %.2f}\}$" % (
    sigfig(p[0] * 1e9, 4), sigfig(perr[0] * 1e9, 4),
    sigfig(p[1], 4), sigfig(perr[1], 4)
)
plt.plot(jx, jy, "-.k", label=eqn)

upr_p = p + perr
lwr_p = p - perr

upper = jmak_y(f_jmak(t_hat, *upr_p))
lower = jmak_y(f_jmak(t_hat, *lwr_p))

it = np.argsort(jx)
plt.fill_between(
    jx[it], upper[it], jy[it], edgecolor=None, facecolor="silver", zorder=1, label=None
)
plt.fill_between(
    jx[it], lower[it], jy[it], edgecolor=None, facecolor="silver", zorder=1, label=None
)

y_naive = f_jmak(t_hat, *p_naive)

jy_naive = jmak_y(y_naive)
eqn_naive = "$1-\\exp\{-(%.1f \\pm %.1f) \\times 10^{-9} \\times t^{%.2f \\pm %.2f}\}$" % (
    sigfig(p_naive[0] * 1e9, 4), sigfig(p_nstd[0] * 1e9, 4),
    sigfig(p_naive[1], 4), sigfig(p_nstd[1], 4)
)
plt.plot(jx, jy_naive, "-.b", label=eqn_naive)

tmin, tmax = plt.xlim()
plt.xlim([0, tmax])
plt.ylim([-10, 2])
legend = plt.legend(loc="best", fontsize="small")
plt.savefig("avrami.png", dpi=400, bbox_inches="tight")
plt.close()

# === Linear Plots ===

plt.figure(figsize=figsize)
plt.title(title)
plt.xlabel("$t$")
plt.ylabel("$Y$")

for i, df in enumerate(frames):
    df = df[df["time"] > 0]
    df = df[df["fraction"] > 0]
    run = "Run {}".format(letters[i + 26])
    plt.plot(df["time"], df["fraction"], label=None)

plt.plot(t_hat, y_hat, "-.k", label=eqn)

upr_p = p + perr
lwr_p = p - perr

upper = f_jmak(t_hat, *upr_p)
lower = f_jmak(t_hat, *lwr_p)

it = np.argsort(t_hat)
plt.fill_between(
    t_hat[it], upper[it], y_hat[it], edgecolor=None, facecolor="silver", zorder=1, label=None
)
plt.fill_between(
    t_hat[it], lower[it], y_hat[it], edgecolor=None, facecolor="silver", zorder=1, label=None
)

plt.plot(t_hat, y_naive, "-.b", label=eqn_naive)

"""
upr_p = p_naive + p_nstd
lwr_p = p_naive - p_nstd

upper = f_jmak(t_hat, *upr_p)
lower = f_jmak(t_hat, *lwr_p)

it = np.argsort(t_hat)
plt.fill_between(
    t_hat[it], upper[it], y_naive[it], edgecolor=None, facecolor="silver", zorder=1, label=None, alpha=0.5
)
plt.fill_between(
    t_hat[it], lower[it], y_naive[it], edgecolor=None, facecolor="silver", zorder=1, label=None, alpha=0.5
)
"""

tmin, tmax = plt.xlim()
plt.xlim([0, tmax])
plt.ylim([0, 1])
legend = plt.legend(loc="best", fontsize="small")
plt.savefig("linear.png", dpi=400, bbox_inches="tight")
plt.close()
