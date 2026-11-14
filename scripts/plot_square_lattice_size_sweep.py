import csv
from collections import defaultdict

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from square_lattice_size_sweep import ROOT


def main():
    output = ROOT / "outputs" / "square_lattice_size_sweep"
    groups = defaultdict(dict)

    with (output / "full_runs.csv").open(encoding="utf-8", newline="") as stream:
        for row in csv.DictReader(stream):
            size = int(row["L"])
            repetition = int(row["repetition"])
            residual = float(row["overlap_residual"])
            if repetition in groups[size]:
                raise ValueError(f"Duplicate run: L={size}, repetition={repetition}")
            if int(row["N"]) != size**2 or not np.isfinite(residual) or residual < 0:
                raise ValueError(f"Invalid row: {row}")
            groups[size][repetition] = residual

    sizes = sorted(groups)
    if sizes != list(range(3, 40)) or any(
        set(groups[size]) != set(range(20)) for size in sizes
    ):
        raise ValueError("Expected repetitions 0–19 at every L from 3 to 39")

    fig, ax = plt.subplots(figsize=(7, 4.5))
    missing_lower = []
    zero_means = []

    for size in sizes:
        values = np.array(list(groups[size].values()))
        mean = values.mean()
        sd = values.std(ddof=1)
        if mean == 0:
            zero_means.append(size)
            continue

        centre = np.log(mean)
        upper = np.log(mean + sd) - centre
        if mean > sd:
            lower = centre - np.log(mean - sd)
        else:
            lower = 0.0
            missing_lower.append(size)

        ax.errorbar(
            size**2, centre, yerr=[[lower], [upper]],
            fmt="o", markersize=4, capsize=3, color="C0",
        )

    ax.set(xlabel="Number of sites N", ylabel=r"$\ln(\overline{1-\alpha})$")
    ax.grid(alpha=0.25)
    fig.tight_layout()
    path = output / "overlap_vs_size_log.png"
    fig.savefig(path, dpi=200)
    plt.close(fig)

    print(f"Wrote {path}")
    print(f"Sizes with a nonpositive lower bound: {missing_lower}")
    print(f"Sizes with zero mean residual: {zero_means}")


if __name__ == "__main__":
    main()