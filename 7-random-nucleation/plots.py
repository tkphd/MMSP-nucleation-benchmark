# Plotting script for PFHub Nucleation Benchmark
# Questions/comments to trevor.keller@nist.gov (Trevor Keller)

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

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

# === Energy ===

plt.figure(figsize=figsize)
plt.title(title)
plt.xlabel("Time $\\tilde{t}$")
plt.ylabel("Free Energy $\\tilde{\\mathcal{F}}$")

for i, df in enumerate(frames):
    plt.plot(df["time"], df["energy"], color=colors[i], label=labels[i])

plt.xlim(tlim)
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

# === Particles ===

plt.figure(figsize=figsize)
plt.title(title)
plt.xlabel("Time $\\tilde{t}$")
plt.ylabel("Particle Count")

for i, df in enumerate(frames):
    plt.plot(df["time"], df["particles"], color=colors[i], label=labels[i])

plt.xlim(tlim)
plt.legend(loc="best")
plt.savefig("particle-count.png", dpi=400, bbox_inches="tight")
plt.close()
