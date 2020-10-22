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

title = "PFHub Benchmark 8.3"
tlim = [0, 600]

labels = ["Run A",
          "Run B",
          "Run C",
          "Run D"]

frames = [pd.read_csv("run-a/free_energy.csv"),
          pd.read_csv("run-b/free_energy.csv"),
          pd.read_csv("run-c/free_energy.csv"),
          pd.read_csv("run-d/free_energy.csv")]

colors = ["red", "green", "blue", "orange"]

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

# === Avrami ===

plt.figure(figsize=figsize)
plt.title(title)
plt.xlabel("$\\log(t)$")
plt.ylabel("$\\log(-\\log(1-Y))$")

fit_t = np.array([])
fit_y = np.array([])

for i, df in enumerate(frames):
    df = df[df["time"] > 0]
    df = df[df["fraction"] > 0]
    x = jmak_x(df["time"])
    y = jmak_y(df["fraction"])
    plt.plot(x, y, color=colors[i], label=labels[i])

    fit_t = np.append(fit_t, df["time"])
    fit_y = np.append(fit_y, df["fraction"])

# === Levenburg-Marquardt Least-Squares Fit ===

p, pcov = curve_fit(f_jmak, fit_t, fit_y, sigma=None, method="lm", jac=df_jmak, maxfev=1000)
print("coeff: ", p, "; covar:")
coef = describe(pcov)
print(coef)

fit_max = np.amax(fit_t)
fit_min = np.exp(floor(np.log(fit_max) / 2) - 1)

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
