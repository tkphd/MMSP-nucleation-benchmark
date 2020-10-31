# Plotting script for PFHub Nucleation Benchmark
# Questions/comments to trevor.keller@nist.gov (Trevor Keller)

# Y      = 1 - exp(-K*(t - t0)**n)
# dY/dK  = (t - t0)**n*exp(-K*(t - t0)**n)
# dY/dn  = K*(t - t0)**n*exp(-K*(t - t0)**n)*log(t - t0)
# dY/dt0 = -K*n*(t - t0)**n*exp(-K*(t - t0)**n)/(t - t0)

from math import floor
from matplotlib import style
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.optimize import curve_fit
from scipy.stats import chisquare, describe

style.use("seaborn")

title = "PFHub Benchmark 8.3"
tlim = [0, 600]
p0 = (5.0e-8, 3.0, 0.0) # initial guess for non-linear solver

labels = ["Run A",
          "Run B",
          "Run C",
          "Run D",
          "Run E",
          "Run F",
          "Run G",
          "Run H",
          "Run I",
          #"Run J",
]

frames = [pd.read_csv("run-a/free_energy.csv"),
          pd.read_csv("run-b/free_energy.csv"),
          pd.read_csv("run-c/free_energy.csv"),
          pd.read_csv("run-d/free_energy.csv"),
          pd.read_csv("run-e/free_energy.csv"),
          pd.read_csv("run-f/free_energy.csv"),
          pd.read_csv("run-g/free_energy.csv"),
          pd.read_csv("run-h/free_energy.csv"),
          pd.read_csv("run-i/free_energy.csv"),
          #pd.read_csv("run-j/free_energy.csv"),
]

figsize=(10,6)

# === Equations ===

def f_jmak(t, K, n, t0):
    # JMAK growth law, Y(t) = 1 - exp(-Ktⁿ)
    # where $n$ is the spatial dimensionality
    return 1.0 - np.exp(-K * (t - t0)**n)

def df_jmak(t, K, n, t0):
    # Jacobian: df/dp for p=(K, n)
    return np.array([
        (t - t0)**n * np.exp(-K * (t - t0)**n),
        K * (t - t0)** n * np.exp(-K * (t - t0)**n) * np.log(t - t0),
        -K * n * (t - t0)**n * np.exp(-K * (t - t0)**n) / (t - t0),
    ]).T

def jmak_x(x):
    return np.log(x)

def jmak_y(y):
    return np.log(-np.log(1 - y))

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

    plt.plot(jmak_x(t), jmak_y(y), label=labels[i])

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
    print(labels[i], " coeffs: ", p)

    K.append(p[0])
    n.append(p[1])
    t0.append(p[2])

    fit_t = np.append(fit_t, t)
    fit_y = np.append(fit_y, y)

print()
p_naive = np.array([np.average(K), np.average(n), np.average(t0)])
p_nstd = np.array([np.std(K), np.std(n), np.std(t0)])
print("Individual fit: K={0:5.3e} n={1:5.3f} t0={2:5.3f}".format(*p_naive))
print("         stdev:   {0:5.3e}   {1:5.3f}    {2:5.3f}".format(*p_nstd))

p, pcov = curve_fit(f_jmak, fit_t, fit_y, p0=p0, sigma=None,
                    method="lm", jac=df_jmak, maxfev=2000)
perr = np.sqrt(np.diag(pcov))

print()
print("Collective fit: K={0:5.3e} n={1:5.3f}  t0={2:5.3f}".format(*p))
print("         error:   {0:5.3e}   {1:5.3f}     {2:5.3f}".format(*perr))

fit_max = np.amax(fit_t)
fit_min = np.exp(floor(np.log(fit_max) / 3))

t_hat = np.linspace(fit_max, fit_min, 201)
y_hat = f_jmak(t_hat, *p)

jx = jmak_x(t_hat)
jy = jmak_y(y_hat)
eqn = "$1-\\exp(-%5.3f \\times 10^{-9} \\times (t-%5.3f)^{%5.3f})$" % (p[0] * 1e9, p[2], p[1])
plt.plot(jx, jy, "-.k", label=eqn)

t_naive = np.linspace(fit_max, fit_min, 201)
y_naive = f_jmak(t_naive, *p_naive)

jx_naive = jmak_x(t_naive)
jy_naive = jmak_y(y_naive)
eqn_naive = "$1-\\exp(-%5.3f \\times 10^{-9} \\times (t-%5.3f)^{%5.3f})$" % (p_naive[0] * 1e9, p_naive[2], p_naive[1])
plt.plot(jx_naive, jy_naive, "-.b", label=eqn_naive)

tmin, tmax = plt.xlim()
plt.xlim([0, tmax])
plt.ylim([-10, 2])
plt.legend(loc="best")
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
    plt.plot(df["time"], df["fraction"], label=labels[i])

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

plt.plot(t_naive, y_naive, "-.b", label=eqn_naive)

upr_p = p_naive + p_nstd
lwr_p = p_naive - p_nstd

upper = f_jmak(t_naive, *upr_p)
lower = f_jmak(t_naive, *lwr_p)

it = np.argsort(t_naive)
plt.fill_between(
    t_naive[it], upper[it], y_naive[it], edgecolor=None, facecolor="silver", zorder=1, label=None, alpha=0.5
)
plt.fill_between(
    t_naive[it], lower[it], y_naive[it], edgecolor=None, facecolor="silver", zorder=1, label=None, alpha=0.5
)

tmin, tmax = plt.xlim()
plt.xlim([0, tmax])
plt.ylim([0, 1])
plt.legend(loc="best")
plt.savefig("linear.png", dpi=400, bbox_inches="tight")
plt.close()
