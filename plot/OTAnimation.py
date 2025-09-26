import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
import seaborn as sns
import re
import matplotlib.animation as animation

source_dir = "GEMData/"
save_name = "GEMAni"

files = [f for f in os.listdir(source_dir) if not f.endswith(".png")][::5]


def extract_time(filename):
    match = re.search(r"t_([0-9.]+)", filename)
    return float(match.group(1)) if match else float("inf")


sorted_files = sorted(files, key=extract_time)

column_names = [
    "x",
    "y",
    "density",
    "v_x",
    "v_y",
    "v_z",
    "pressure",
    "B_x",
    "B_y",
    "B_z",
    "psi",
    "J_z",
]


def convert_to_2D(row, num_xvalues, num_yvalues):
    return row.to_numpy().reshape(num_xvalues, num_yvalues).T


fig, ax = plt.subplots(figsize=(8, 6))

# initialize with the first frame
first_data = pd.read_csv(
    os.path.join(source_dir, sorted_files[0]),
    sep=r"\s+",
    header=None,
    skip_blank_lines=True,
    names=column_names,
)

xdata = np.unique(first_data["x"])
ydata = np.unique(first_data["y"])
num_x = len(xdata)
num_y = len(ydata)

density = convert_to_2D(first_data["J_z"], num_x, num_y)
hm = sns.heatmap(density, cmap="magma", cbar=False, ax=ax)
ax.set_xticks([])
ax.set_yticks([])


def update(frame_idx):
    ax.clear()
    file = sorted_files[frame_idx]
    raw_data = pd.read_csv(
        os.path.join(source_dir, file),
        sep=r"\s+",
        header=None,
        skip_blank_lines=True,
        names=column_names,
    )
    density = convert_to_2D(raw_data["J_z"], num_x, num_y)
    sns.heatmap(density, cmap="magma", cbar=False, ax=ax)
    ax.set_xticks([])
    ax.set_yticks([])
    print(f"Completed Frame {frame_idx}/{len(sorted_files)}")
    return []


ani = animation.FuncAnimation(
    fig, update, frames=len(sorted_files), blit=False, repeat=False
)

plt.tight_layout()
ani.save(f"{save_name}.mp4", writer="ffmpeg", fps=10)
