import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
import seaborn as sns
import re
import matplotlib.animation as animation

source_dir = "BrioWu/"
save_dir = "PlotsBrioWu/"

files = [f for f in os.listdir(source_dir) if not f.endswith(".png")]


def extract_time(filename):
    match = re.search(r"t_([0-9.]+)", filename)
    return float(match.group(1)) if match else float("inf")


sorted_files = sorted(files, key=extract_time)

column_names = [
    "x",
    "y",
    "density",
    "velocity",
    "v_y",
    "v_z",
    "pressure",
    "B_x",
    "B_y",
    "B_z",
    "psi",
    "J_z",
]

data = pd.read_csv(
    os.path.join(source_dir, sorted_files[-1]),
    sep=r"\s+",
    header=None,
    skip_blank_lines=True,
    names=column_names,
)
xdata = np.unique(data["x"])
num_x = len(xdata)

ydata = np.unique(data["y"])
num_y = len(ydata)

for i, column in enumerate(column_names[2:]):

    fig, ax = plt.subplots(figsize=(8, 6))
    plt.plot(xdata[:-4], data[column][:: num_y + 1])
    plt.ylabel(column)
    plt.xlabel("x")
    plt.tight_layout()
    plt.savefig(os.path.join(save_dir, column + "_plot"))
