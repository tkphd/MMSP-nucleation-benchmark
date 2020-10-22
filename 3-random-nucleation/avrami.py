# Plotting script for PFHub Nucleation Benchmark
# Questions/comments to trevor.keller@nist.gov (Trevor Keller)

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib import style
style.use("seaborn")

title = "PFHub Benchmark 8.3"

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

# === Avrami ===

plt.figure(figsize=figsize)
plt.title(title)
plt.xlabel("$\\log(t)$")
plt.ylabel("$\\log(-\\log(1-Y))$")

for i, df in enumerate(frames):
    df = df[df["time"] > 0]
    df = df[df["fraction"] > 0]
    x = np.log(df["time"])
    y = np.log(-np.log(1 - df["fraction"]))
    plt.plot(x, y, color=colors[i], label=labels[i])

# plt.ylim([-5, 1])
plt.legend(loc="best")
plt.savefig("avrami.png", dpi=400, bbox_inches="tight")
plt.close()
