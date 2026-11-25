# PCA and the ferromagnetic Ising model

Code companion to: <br><br>
**When Does PCA Recover the Magnetisation? The Principal Component of the Ferromagnetic Ising
Model**<br>
*Roberto C. Alamino<br>
J. Stat. Mech (2026) 113405*<br>
[DOI 10.1088/1742-5468/aec121](https://doi.org/10.1088/1742-5468/aec121)


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

See [Simulation and analysis algorithms](docs/algorithms.md) for the
samplers, covariance estimators, and plotted observables.

The document [Mathematical clarifications](docs/mathematical-clarifications.md) presents refinements of the work concerning the scope of the finite-tree result, a lower bound for the infinite-size argument, and notes on covariance and temperature limits which were prepared during the proof and not published together.

## Reproduce the Analyses

Run the commands from the repository root with the virtual environment
active. Figure numbers refer to the paper. Where a workflow lists
multiple commands, run them in the order shown.

| Paper figure | Analysis | Commands |
| --- | --- | --- |
| Figure 1 | PC1 overlap with the magnetisation direction as a function of periodic square-lattice size. The companion and the archived paper routine aggregate the log residual differently; see the [Figure 1 note](docs/algorithms.md#figure-1-square-lattice-size-sweep). | `python scripts/square_lattice_size_sweep.py --full`<br>`python scripts/plot_square_lattice_size_sweep.py` |
| Figure 2 | Site-by-site PC1 on a periodic 20 × 20 square lattice. | `python scripts/square_lattice_pc1.py` |
| Figure 3 | PC1 overlap as a function of Cayley tree size, with an exact covariance baseline. | `python scripts/cayley_tree_size_sweep.py --full` |
| Figures 4 and 5 | Radial map and generation profile of PC1 on a 255-site Cayley tree. Both views use the same sampled PC1 vector. | `python scripts/cayley_tree_pc1.py` |
| Figure 6 | Comparison of simulated and exact PC1 ratios near the root as the buffer grows. The script reads the archived simulation profiles in `data/roi/df_profiles.csv` and recalculates the exact curves. | `python scripts/roi_mix.py` |
| Figure 7 | Extended calculation of the exact PC1 ratios for buffer sizes 1–50. | `python scripts/roi_sim.py` |
| Figure 8 | Temperature dependence of PC1 overlap, relative spectral gap, and mean absolute magnetisation on both geometries. | `python scripts/temperature_sweep.py --full-lattice`<br>`python scripts/temperature_sweep.py --full-tree`<br>`python scripts/plot_temperature_sweep.py` |

Workflows that perform simulations generate fresh samples using fixed seeds.

**Figure 3 — Cayley tree size sweep.** The companion reproduces the
declining trend in PC1 overlap, with higher overlaps at the largest
sizes shown. See the [Figure 3 comparison](docs/algorithms.md#figure-3-cayley-tree-size-sweep)
for the numerical differences and the limits of the sampling explanation.

**Figure 6 — Mixed region-of-interest comparison.** The script
reconstructs the simulation points from the included profiles and
recalculates the exact curves. It does not rerun the simulations
that generated those profiles.

## Included files and generated outputs

- `scripts/` contains the simulation, analysis, and plotting code
  referenced in the figure table above.
- `requirements.txt` records the Python dependencies used by the scripts.
- `data/roi/df_profiles.csv` contains archived PC1 summaries by tree
  generation, buffer size, and simulation repetition. It is included
  because the Figure 6 workflow uses these summaries to reconstruct
  the simulation markers and error bars. The script divides each
  generation mean by the root value from the same run and recalculates
  the exact curves independently.
- `.gitignore` keeps local virtual environments, Python caches, and
  generated outputs out of Git.

The figure commands generate plots and numerical data under `outputs/`.
In particular, the workflow for Figures 4 and 5 exports
`outputs/cayley_tree_pc1/pc1_by_node.csv` for the radial map and
`outputs/cayley_tree_pc1/pc1_by_generation.csv` for the generation
profile. The Figure 6 workflow exports the plotted statistics to
`outputs/roi_mix/roi_mix_summary.csv`.