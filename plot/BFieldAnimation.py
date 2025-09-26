import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
import seaborn as sns
import re
import matplotlib.animation as animation
import scipy.sparse
import scipy.linalg


"""Plotting contours of A_z i.e. magnetic field lines
For a 2D system:    
    B_x = dA_z/dy
    B_y = -dA_z/dx

Laplacian(A_z)=-J_z (in reduced units)
Use scipy sparse solver to calc A_z (periodic BCs) """

source_dir = "GEMData/"
save_name = "GEMFieldAni"

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
# Data is stored in column major order (kinda annoying but oh well)


# initialize with the first frame
first_data = pd.read_csv(
    os.path.join(source_dir, sorted_files[0]),
    sep=r"\s+",
    header=None,
    skip_blank_lines=True,
    names=column_names,
)

xdata = np.unique(first_data["x"].to_numpy())
ydata = np.unique(first_data["y"].to_numpy())
num_x = len(xdata)
num_y = len(ydata)

X, Y = np.meshgrid(xdata, ydata)

dx = np.mean(np.diff(xdata))
dy = np.mean(np.diff(ydata))

print(f"dx={dx}, dy ={dy} \n num_x = {num_x}, num_y = {num_y}")


# Creating a suitable sparse Laplacian matrix in scipy
# Perform 1d Laplacian for x and y directions with customisable BCs then combine to form a single 2D laplacian
def laplacian_1d(n, dx, bc="dirichlet"):
    """Create 1D Laplacian with nx points and given BC."""

    # Construct 1D Laplacian array by diagonals "lil" format is used purely to speedup BC assignment
    L = (
        scipy.sparse.diags(
            [-np.ones(n - 1), 2 * np.ones(n), -np.ones(n - 1)],
            offsets=[-1, 0, 1],
            format="lil",
        )
        / dx**2
    )

    if bc == "dirichlet":
        L[0, :] = 0
        L[0, 0] = 1
        L[-1, :] = 0
        L[-1, -1] = 1
    elif bc == "neumann":
        L[0, :] = 0
        L[0, 0] = 1
        L[0, 1] = -1
        L[-1, :] = 0
        L[-1, -1] = 1
        L[-1, -2] = -1
    elif bc == "periodic":
        L[0, -1] = -1 / dx**2
        L[-1, 0] = -1 / dx**2
    else:
        raise ValueError("BC must be 'dirichlet', 'neumann', or 'periodic'")

    return L.tocsr()  # Return in csr format for fast spsolve


# Forming 2D laplacian and solving for A_z
Jz_corder = first_data["J_z"].to_numpy()


def solve_A_z(
    J_zcorder, nx=num_x, ny=num_y, dx=dx, dy=dy, bx="dirichlet", by="dirichlet"
):
    Jz_rorder = np.ravel(
        Jz_corder.reshape((nx, ny), order="F"), order="C"
    )  # Convert to C style row major order

    Lx = laplacian_1d(nx, dx, bx)
    Ly = laplacian_1d(ny, dy, by)

    Ix = scipy.sparse.identity(nx, format="csr")
    Iy = scipy.sparse.identity(ny, format="csr")

    L2 = scipy.sparse.kron(Iy, Lx) + scipy.sparse.kron(Ly, Ix)

    # Solve sparse system
    A_flat = scipy.sparse.linalg.spsolve(L2, Jz_rorder)
    A = A_flat.reshape((nx, ny))
    return A


A_data = []
for i, file in enumerate(sorted_files):
    data = pd.read_csv(
        os.path.join(source_dir, file),
        sep=r"\s+",
        header=None,
        skip_blank_lines=True,
        names=column_names,
    )
    Jz_corder = data["J_z"].to_numpy()
    A = solve_A_z(Jz_corder, bx="periodic", by="neumann")
    print(f"Solved A for step {i}/{len(sorted_files)}")
    A_data.append(A)


fig, ax = plt.subplots(figsize=(8, 6))
contours = ax.contour(X, Y, A_data[0].T, levels=20, colors="k", linestyles="solid")
ax.set_aspect("equal")
ax.axis("off")


def update(frame_idx):
    ax.clear()
    # Draw new contours
    plt.contour(X, Y, A_data[frame_idx].T, levels=20, colors="k", linestyles="solid")
    ax.set_aspect("equal")
    ax.axis("off")
    print(f"Completed Frame {frame_idx}/{len(sorted_files)}")
    return []


ani = animation.FuncAnimation(
    fig, update, frames=len(sorted_files), blit=False, repeat=False
)

plt.tight_layout()
ani.save(f"{save_name}.mp4", writer="ffmpeg", fps=10)
