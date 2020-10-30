# Plotting script for PFHub Nucleation Benchmark
# Questions/comments to trevor.keller@nist.gov (Trevor Keller)

from math import floor
from matplotlib import style
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.optimize import curve_fit
from scipy.stats import chisquare, describe

style.use("seaborn")

title = "MODIFIED Benchmark 8.3"
tlim = [0, 600]
p0 = (5.0e-8, 3.0) # initial guess for non-linear solver

labels = ["Run A",
          "Run B",
          "Run C",
          "Run D",
          #"Run E",
          #"Run F",
          #"Run G",
          #"Run H",
          #"Run I",
          #"Run J",
          #"Run K",
          #"Run L",
]

frames = [pd.read_csv("run-a/free_energy.csv"),
          pd.read_csv("run-b/free_energy.csv"),
          pd.read_csv("run-c/free_energy.csv"),
          pd.read_csv("run-d/free_energy.csv"),
          #pd.read_csv("run-e/free_energy.csv"),
          #pd.read_csv("run-f/free_energy.csv"),
          #pd.read_csv("run-g/free_energy.csv"),
          #pd.read_csv("run-h/free_energy.csv"),
          #pd.read_csv("run-i/free_energy.csv"),
          #pd.read_csv("run-j/free_energy.csv"),
          #pd.read_csv("run-k/free_energy.csv"),
          #pd.read_csv("run-l/free_energy.csv"),
]

figsize=(10,6)

# === Equations ===

def f_jmak(t, K, n):
    # JMAK growth law, Y(t) = 1 - exp(-Ktⁿ)
    # where $n$ is the spatial dimensionality
    return 1.0 - np.exp(-K * t**n)

def df_jmak(t, K, n):
    # Jacobian: df/dp for p=(K, n)
    return np.array([
        t**n * np.exp(-K * t**n),
        K * t**n * np.log(t) * np.exp(-K * t**n)
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

    t = np.array(df["time"])
    y = np.array(df["fraction"])

    plt.plot(jmak_x(t), jmak_y(y), label=labels[i])

# === Levenburg-Marquardt Least-Squares Fit ===

fit_t = np.array([])
fit_y = np.array([])

K = []
n = []

for i, df in enumerate(frames):
    df = df[df["time"] > 0]
    df = df[df["fraction"] > 0]

    t = np.array(df["time"])
    y = np.array(df["fraction"])

    # Fit this dataset & print coeffs
    p, pcov = curve_fit(f_jmak, t, y, p0=p0, sigma=None,
                        method="lm", jac=df_jmak, maxfev=2000)
    print(labels[i], " coeffs: ", p)

    K.append(p[0])
    n.append(p[1])

    fit_t = np.append(fit_t, t)
    fit_y = np.append(fit_y, y)

print()
p_naive = np.array([np.average(K), np.average(n)])
p_nstd = np.array([np.std(K), np.std(n)])
print("Individual fit: K={0:.3e} n={1:.3e}".format(*p_naive))
print("         stdev:   {0:.3e}   {1:.3e}".format(*p_nstd))

p, pcov = curve_fit(f_jmak, fit_t, fit_y, p0=p0, sigma=None,
                    method="lm", jac=df_jmak, maxfev=2000)
perr = np.sqrt(np.diag(pcov))

print()
print("Collective fit: K={0:.3e} n={1:.3e}".format(*p))
print("         error:   {0:.3e}   {1:.3e}".format(*perr))

fit_max = np.amax(fit_t)
fit_min = np.exp(floor(np.log(fit_max) / 3))

t_hat = np.linspace(fit_max, fit_min, 201)
y_hat = f_jmak(t_hat, *p)

jx = jmak_x(t_hat)
jy = jmak_y(y_hat)
eqn = "$1-\\exp(-%.4g \\times t^{%.4g})$" % (p[0], p[1])
plt.plot(jx, jy, "-.k", label=eqn)

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

tmin, tmax = plt.xlim()
plt.xlim([0, tmax])
ymin, ymax = plt.ylim()
plt.ylim([0, ymax])
plt.legend(loc="best")
plt.savefig("linear.png", dpi=400, bbox_inches="tight")
plt.close()
