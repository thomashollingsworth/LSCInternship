import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
import seaborn as sns
import re
import matplotlib.animation as animation

"""Data is arranged as  x0 y0
                        x1 y0
                        x2 y0
                        x3 y0
                        etc.

                        x1 y0
                        ...
                        """

source_dir = (
    "/Users/tomhollingsworth/Documents/LSC Internship/ResistiveMHDProject/ResistiveOT"
)
save_name = (
    "/Users/tomhollingsworth/Documents/LSC Internship/ResistiveMHDProject/Res_B_Ani"
)

files = [f for f in os.listdir(source_dir) if not f.endswith(".png")]


def extract_time(filename):
    match = re.search(r"t_([0-9.]+)", filename)
    return float(match.group(1)) if match else float("inf")


# Sort files by extracted time
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
]

Bx_arrays = []
By_arrays = []


def convert_to_2D(row, num_xvalues, num_yvalues):
    array = row[:].to_numpy().reshape((num_yvalues, num_xvalues))
    return array


for file in sorted_files:
    raw_data = pd.read_csv(
        os.path.join(
            source_dir,
            file,
        ),
        sep=r"\s+",
        header=None,
        skip_blank_lines=True,
        names=column_names,
    )

    num_x = len(np.unique(raw_data["x"]))
    num_y = len(np.unique(raw_data["y"]))

    Bx_data = convert_to_2D(raw_data["B_x"], num_x, num_y)
    By_data = convert_to_2D(raw_data["B_y"], num_x, num_y)
    Bx_arrays.append(Bx_data)
    By_arrays.append(By_data)

num_x = len(np.unique(raw_data["x"]))
num_y = len(np.unique(raw_data["y"]))
x_data = np.arange(num_x)
y_data = np.arange(num_y)

X, Y = np.meshgrid(x_data, y_data)


fig, ax = plt.subplots(figsize=(8, 6))


def init():
    plt.figure(figsize=(12, 7))
    ax.streamplot(
        X,
        Y,
        Bx_arrays[0],
        By_arrays[0],
        color="k",  # black streamlines
        density=1.2,
        linewidth=0.8,
    )
    plt.title("B Field")
    ax.set_xticks([])
    ax.set_yticks([])
    return []


def update(frame):
    ax.clear()
    ax.streamplot(
        X,
        Y,
        Bx_arrays[frame],
        By_arrays[frame],
        color=np.sqrt(Bx_arrays[frame] ** 2 + By_arrays[frame] ** 2),  # color by |B|
        cmap="plasma",
        density=1.2,
        linewidth=0.8,
    )
    plt.title("B Field")
    ax.set_xticks([])
    ax.set_yticks([])
    return []


ani = animation.FuncAnimation(
    fig, update, frames=len(Bx_arrays), init_func=init, blit=False, repeat=False
)

plt.tight_layout()
ani.save(f"{save_name}.mp4", writer="ffmpeg", fps=5)
