import csv
from collections import defaultdict
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
BUFFERS = range(1, 14)
ROI_DEPTH = 4


def load_ratios(path):
    profiles = defaultdict(dict)
    with path.open(encoding="utf-8-sig", newline="") as stream:
        for row in csv.DictReader(stream):
            key = int(row["buffer"]), int(row["repetition"])
            level = int(row["level"])
            if level in profiles[key]:
                raise ValueError(f"Duplicate entry for {key}, level {level}")
            profiles[key][level] = float(row["pc1_mean"])

    expected = {(b, r) for b in BUFFERS for r in range(20)}
    if profiles.keys() != expected:
        raise ValueError("Expected 20 repetitions for each buffer size 1–13")

    ratios = defaultdict(list)
    for (buffer, repetition), profile in profiles.items():
        if not set(range(ROI_DEPTH + 1)) <= profile.keys():
            raise ValueError(f"Missing levels for buffer {buffer}, run {repetition}")
        root = profile[0]
        if not np.isfinite(root) or np.isclose(root, 0):
            raise ValueError(f"Invalid root PC1 for buffer {buffer}, run {repetition}")
        for level in range(1, ROI_DEPTH + 1):
            ratios[buffer, level].append(profile[level] / root)
    return ratios


def exact_ratios(depth):
    t = np.tanh(1.0)  # beta = J = 1
    sizes = 2.0 ** np.arange(depth + 1)
    weights = np.sqrt(sizes)
    reduced = np.zeros((depth + 1, depth + 1))

    # Build the covariance in the orthonormal basis of generation indicators.
    # At generation k there are 2**k independent spin innovations,
    # each contributing to its descendants with a factor t per edge.
    for k in range(depth + 1):
        decay = t ** np.arange(depth + 1 - k)
        v = weights[k:] * decay
        variance = 1.0 if k == 0 else 1.0 - t * t
        reduced[k:, k:] += (variance / sizes[k]) * np.outer(v, v)

    _, vectors = np.linalg.eigh(reduced)
    per_site = vectors[:, -1] / weights
    return per_site[1:ROI_DEPTH + 1] / per_site[0]


def main():
    ratios = load_ratios(ROOT / "data" / "roi" / "df_profiles.csv")
    exact = np.array([exact_ratios(ROI_DEPTH + b) for b in BUFFERS])
    output_dir = ROOT / "outputs" / "roi_mix"
    output_dir.mkdir(parents=True, exist_ok=True)

    table_path = output_dir / "roi_mix_summary.csv"
    with table_path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow((
            "buffer", "level", "simulation_mean_ratio", "simulation_sd_ratio",
            "repetitions", "exact_ratio",
        ))
        for index, buffer in enumerate(BUFFERS):
            for level in range(1, ROI_DEPTH + 1):
                values = np.asarray(ratios[buffer, level])
                writer.writerow((
                    buffer, level, values.mean(), values.std(ddof=1),
                    values.size, exact[index, level - 1],
                ))
    print(f"Wrote {table_path}")

    fig, ax = plt.subplots(figsize=(8, 4.5))
    for level in range(1, ROI_DEPTH + 1):
        samples = np.array([ratios[b, level] for b in BUFFERS])
        colour = f"C{level - 1}"
        ax.errorbar(
            BUFFERS, samples.mean(axis=1), yerr=samples.std(axis=1, ddof=1),
            fmt="o", capsize=4, color=colour, label=f"Simulation ℓ={level}",
        )
        ax.plot(
            BUFFERS, exact[:, level - 1], "--",
            color=colour, label=f"Exact ℓ={level}",
        )

    ax.axhline(1.0, color="0.5", linestyle=":")
    ax.set(xlabel=r"Buffer size $m_b$", ylabel=r"$R_\ell=v_\ell/v_0$")
    ax.set_xticks(list(BUFFERS))
    ax.grid(alpha=0.25)
    ax.legend(loc="center left", bbox_to_anchor=(1, 0.5))
    fig.tight_layout()
    path = output_dir / "roi_mix.png"
    fig.savefig(path, dpi=200)
    plt.close(fig)
    print(f"Plotted 13 buffer sizes × 4 levels, using 20 repetitions per point.")
    print(f"Wrote {path}")


if __name__ == "__main__":
    main()