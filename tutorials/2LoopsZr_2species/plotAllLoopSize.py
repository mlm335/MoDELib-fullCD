import os, re, sys, pathlib
import numpy as np
import matplotlib.pyplot as plt
sys.path.insert(0, '../../lib')

from readNodesMod2 import *
from readFile import *

plt.rc('text', usetex=True)
plt.rc('font', family='Times New Roman')
plt.rcParams['text.latex.preamble'] = r'\usepackage{amsmath}'


def sorted_evl_ids(evl_folder):
    ids = []
    for fname in os.listdir(evl_folder):
        m = re.match(r"evl_(\d+)\.txt$", fname)
        if m:
            ids.append(int(m.group(1)))
    return sorted(ids)


def get_loop_radii(folder):
    """
    Reads all evl_*.txt files and returns:
        time_SI: simulation time in seconds
        radii_nm: array with shape [n_steps, n_loops]
    """

    folder = pathlib.Path(folder)
    evl_folder = folder / "evl"
    F_file = folder / "F" / "F_0.txt"
    mat_file = folder / "inputFiles" / "Zr.txt"

    F = np.loadtxt(F_file)

    b_SI = get_scalar(str(mat_file), "b_SI")
    mu0_SI = get_scalar(str(mat_file), "mu0_SI")
    rho_SI = get_scalar(str(mat_file), "rho_SI")

    v_dd2SI = np.sqrt(mu0_SI / rho_SI)
    t_dd2SI = b_SI / v_dd2SI

    run_ids = F[:, 0].astype(int)
    time_SI = F[:, 1] * t_dd2SI

    all_radii = []

    for runID in run_ids:
        evl_file_no_ext = evl_folder / f"evl_{runID}"
        nd = readNODEStxt(str(evl_file_no_ext))

        loop_areas = np.array(nd.loopsArea)

        # MoDELib area is in b^2 units, so radius = sqrt(A*b^2/pi)
        radii_m = np.sqrt(loop_areas * b_SI**2 / np.pi)
        radii_nm = radii_m * 1e9

        all_radii.append(radii_nm)

    # Pad in case loops disappear/change number
    max_loops = max(len(r) for r in all_radii)
    radii_nm = np.full((len(all_radii), max_loops), np.nan)

    for i, r in enumerate(all_radii):
        radii_nm[i, :len(r)] = r

    return time_SI, radii_nm


# ============================================================
# USER INPUT
# ============================================================

folder = pathlib.Path.cwd()   # run from tutorial folder
# folder = "/Users/matthewmaron/Documents/MoDELib2-05262026/tutorials/OswaldRippeningSemiLumped"

time_SI, radii_nm = get_loop_radii(folder)

# ============================================================
# PLOT
# ============================================================

fig, ax = plt.subplots(figsize=(10, 7))

for loop_id in range(radii_nm.shape[1]):
    ax.plot(time_SI, radii_nm[:, loop_id], marker='o', markersize=3, label=fr"Loop {loop_id}")

ax.set_xlabel(r"Time [s]", fontsize=24)
ax.set_ylabel(r"Loop radius [nm]", fontsize=24)
ax.tick_params(axis="both", which="major", direction="in", top=True, right=True, labelsize=18)
ax.grid(True, linestyle="dotted")
ax.legend(fontsize=14, ncol=2)
ax.set_ylim(0,50)
plt.tight_layout()
plt.savefig("all_loop_radii_vs_time.png", dpi=300)
plt.savefig("all_loop_radii_vs_time.pdf", dpi=300)
plt.show()