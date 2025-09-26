import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
import seaborn as sns
import re
import matplotlib.animation as animation

methods = ["FORCE", "LF", "Richtymer", "Godunov"]
source_dirs = ["FORCE_Method/", "LF_Method/", "Richtymer_Method/", "SLIC_Method"]


def extract_time(filename):
    match = re.search(r"EulerResults_([0-9.]+)", filename)
    return float(match.group(1)) if match else float("inf")


column_names = ["x", "density", "velocity", "pressure"]
fig, ax = plt.subplots(figsize=(10, 6))
plt.xlabel("x")
plt.ylabel("Density")


for i, dir in enumerate(source_dirs):
    files = [f for f in os.listdir(dir) if not f.endswith(".png")]
    sorted_files = sorted(files, key=extract_time)
    final_file = sorted_files[-1]

    data = pd.read_csv(
        os.path.join(dir, final_file),
        sep=r"\s+",
        header=None,
        skip_blank_lines=True,
        names=column_names,
    )
    xdata = data["x"]
    density = data["density"]
    plt.plot(xdata, density, label=methods[i])

plt.legend()
plt.tight_layout()
plt.savefig("Trial.png")
