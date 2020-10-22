# Plotting script for PFHub Nucleation Benchmark
# Questions/comments to trevor.keller@nist.gov (Trevor Keller)

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

title = "PFHub Benchmark 8.1"
tlim = [0, 100]

labels = ["$r_0=0.99 r^*$",
          "$r_0=1.00 r^*$",
          "$r_0=1.01 r^*$"]

frames = [pd.read_csv("run-a/free_energy.csv"),
          pd.read_csv("run-b/free_energy.csv"),
          pd.read_csv("run-c/free_energy.csv")]

colors = ["red", "black", "blue"]

figsize=(10,6)

# === Energy ===

plt.figure(figsize=figsize)
plt.title(title)
plt.xlabel("Time $\\tilde{t}$")
plt.ylabel("Free Energy $\\tilde{\\mathcal{F}}$")

for i, df in enumerate(frames):
    plt.plot(df["time"], df["energy"], color=colors[i], label=labels[i])

plt.xlim(tlim)
plt.ylim([0, 4])
plt.legend(loc="best")
plt.savefig("free-energy.png", dpi=400, bbox_inches="tight")
plt.close()

# === Solid Fraction ===

plt.figure(figsize=figsize)
plt.title(title)
plt.xlabel("Time $\\tilde{t}$")
plt.ylabel("Solid Fraction $Y$")

for i, df in enumerate(frames):
    plt.plot(df["time"], df["fraction"], color=colors[i], label=labels[i])

plt.xlim(tlim)
plt.legend(loc="best")
plt.savefig("solid-fraction.png", dpi=400, bbox_inches="tight")
plt.close()
