import csv
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


OUTPUT = Path(__file__).resolve().parents[1] / "outputs" / "temperature_sweep"
METRICS = ("overlap", "relative_gap", "mean_abs_m")


def summarise_runs(path, beta_max):
    betas = np.linspace(0.001, beta_max, 100)
    runs = np.full((100, 20, len(METRICS)), np.nan)
    seen = set()

    with path.open(encoding="utf-8", newline="") as stream:
        for row in csv.DictReader(stream):
            i = int(row["beta_index"])
            rep = int(row["repetition"])
            if not (0 <= i < 100 and 0 <= rep < 20):
                raise ValueError(f"{path.name}: unexpected run index {(i, rep)}")
            if (i, rep) in seen:
                raise ValueError(f"{path.name}: duplicate run {(i, rep)}")
            if not np.isclose(float(row["beta"]), betas[i], rtol=0, atol=1e-10):
                raise ValueError(f"{path.name}: unexpected beta at index {i}")

            runs[i, rep] = [float(row[key]) for key in METRICS]
            seen.add((i, rep))

    if len(seen) != 2000 or not np.isfinite(runs).all():
        raise ValueError(f"{path.name}: expected 2,000 complete, finite runs")

    return betas, runs.mean(axis=1), runs.std(axis=1, ddof=1)


def main():
    panels = (
        ("Square lattice (15 × 15)", "square_lattice_runs.csv", 1.0,
         np.log(1 + np.sqrt(2)) / 2),
        ("Cayley tree (255 sites)", "cayley_tree_runs.csv", 5.0,
         np.arctanh(0.5)),
    )
    styles = (
        ("Overlap", "tab:blue"),
        ("Relative gap", "tab:red"),
        (r"Mean $|m|$", "tab:green"),
    )

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2), sharey=True)
    for ax, (title, filename, beta_max, beta_c) in zip(axes, panels):
        beta, means, sds = summarise_runs(OUTPUT / filename, beta_max)
        reduced = (beta - beta_c) / beta_c

        for j, (label, colour) in enumerate(styles):
            mean, sd = means[:, j], sds[:, j]
            ax.fill_between(reduced, mean - sd, mean + sd, color=colour, alpha=0.18)
            ax.plot(reduced, mean, "o-", color=colour, markersize=3,
                    markerfacecolor="white", label=label)

        ax.axvline(0, color="0.6", linestyle=":", linewidth=1)
        ax.set(title=title, xlabel=r"$(\beta-\beta_c)/\beta_c$")
        ax.grid(alpha=0.25)
        ax.legend(loc="upper left" if ax is axes[0] else "lower right")

    fig.tight_layout()
    path = OUTPUT / "ogm_reproduction.png"
    fig.savefig(path, dpi=200)
    plt.close(fig)
    print(f"Wrote {path} from 2,000 runs per panel")


if __name__ == "__main__":
    main()