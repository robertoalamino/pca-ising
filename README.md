# PCA and the ferromagnetic Ising model

Code companion to Roberto C. Alamino, *When Does PCA Recover the
Magnetisation? The Principal Component of the Ferromagnetic Ising
Model*, accepted for publication in the *Journal of Statistical
Mechanics: Theory and Experiment*.

The paper studies when the first principal component (PC1) of Ising
spin configurations aligns with the magnetisation direction. The code
covers the periodic square lattice and the binary Cayley tree.

## Setup

The scripts were run with Python 3.15.2. From the repository root,
create and activate a virtual environment. On Windows with Git Bash:

```bash
py -3.15 -m venv .venv
source .venv/Scripts/activate
python -m pip install -r requirements.txt
```

On Linux or macOS, use `python3.15 -m venv .venv` and
`source .venv/bin/activate`, then run the same installation command.

## Reproduce the analyses

Run the commands from the repository root with the virtual environment
active. Figure numbers refer to the paper. Where a workflow lists
multiple commands, run them in the order shown.

| Paper figure | Analysis | Commands |
| --- | --- | --- |
| Figure 1 | PC1 overlap with the magnetisation direction as a function of periodic square-lattice size. | `python scripts/square_lattice_size_sweep.py --full`<br>`python scripts/plot_square_lattice_size_sweep.py` |
| Figure 2 | Site-by-site PC1 on a periodic 20 × 20 square lattice. | `python scripts/square_lattice_pc1.py` |
| Figure 3 | PC1 overlap as a function of Cayley tree size, with an exact covariance baseline. | `python scripts/cayley_tree_size_sweep.py --full` |
| Figures 4 and 5 | Radial map and generation profile of PC1 on a 255-site Cayley tree. Both views use the same sampled PC1 vector. | `python scripts/cayley_tree_pc1.py` |
| Figure 6 | Comparison of simulated and exact PC1 ratios near the root as the buffer grows. The script reads the archived simulation profiles in `data/roi/df_profiles.csv` and recalculates the exact curves. | `python scripts/roi_mix.py` |
| Figure 7 | Extended calculation of the exact PC1 ratios for buffer sizes 1–50. | `python scripts/roi_sim.py` |
| Figure 8 | Temperature dependence of PC1 overlap, relative spectral gap, and mean absolute magnetisation on both geometries. | `python scripts/temperature_sweep.py --full-lattice`<br>`python scripts/temperature_sweep.py --full-tree`<br>`python scripts/plot_temperature_sweep.py` |

The stochastic workflows generate fresh samples using fixed seeds.
For Figure 3, the overlap decreases with tree size, but the largest
trees give higher sampled overlaps than those shown in the paper.
For Figure 6, the script reconstructs the simulation points from
the included profiles; it does not rerun the simulations that
generated those profiles.