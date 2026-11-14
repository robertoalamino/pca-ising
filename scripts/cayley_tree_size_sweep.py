import argparse
import csv

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from roi_mix import ROOT


T = np.tanh(1.0)  # beta = J = 1


def exact_overlap(depth):
    """Compute the population overlap in the generation basis."""
    sizes = 2.0 ** np.arange(depth + 1)
    weights = np.sqrt(sizes)
    covariance = np.zeros((depth + 1, depth + 1))

    for k in range(depth + 1):
        v = weights[k:] * T ** np.arange(depth + 1 - k)
        variance = 1.0 if k == 0 else 1.0 - T**2
        covariance[k:, k:] += (variance / sizes[k]) * np.outer(v, v)

    _, vectors = np.linalg.eigh(covariance)
    return abs(vectors[:, -1] @ weights) / np.linalg.norm(weights)


def sampled_overlap(depth, n_samples, rng):
    n_sites = 2 ** (depth + 1) - 1
    x = np.empty((n_samples, n_sites), dtype=np.float64)
    x[:, 0] = rng.choice([-1, 1], size=n_samples)

    for level in range(1, depth + 1):
        first = 2**level - 1
        stop = 2 ** (level + 1) - 1
        parents = (np.arange(first, stop) - 1) // 2
        signs = np.where(
            rng.random((n_samples, stop - first)) < (1 + T) / 2,
            1.0, -1.0,
        )
        x[:, first:stop] = x[:, parents] * signs

    # Centred PCA. The leading right singular vector is PC1.
    x -= x.mean(axis=0)
    _, _, vt = np.linalg.svd(x, full_matrices=False)
    return abs(vt[0].sum()) / np.sqrt(n_sites)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--full", action="store_true")
    args = parser.parse_args()

    mode = "full" if args.full else "smoke"
    depths = range(1, 11) if args.full else (2, 4, 7)
    repetitions, n_samples = (20, 500) if args.full else (3, 100)
    output = ROOT / "outputs" / "cayley_tree_size_sweep"
    output.mkdir(parents=True, exist_ok=True)
    summary = []

    with (output / f"{mode}_runs.csv").open(
        "w", encoding="utf-8", newline=""
    ) as stream:
        writer = csv.writer(stream)
        writer.writerow(("depth", "N", "repetition", "sampled_overlap", "exact_overlap"))

        for depth in depths:
            n_sites = 2 ** (depth + 1) - 1
            exact = exact_overlap(depth)
            values = []

            for repetition in range(repetitions):
                rng = np.random.default_rng([2026, depth, repetition])
                overlap = sampled_overlap(depth, n_samples, rng)
                values.append(overlap)
                writer.writerow((depth, n_sites, repetition, overlap, exact))

            mean = np.mean(values)
            sd = np.std(values, ddof=1)
            summary.append((n_sites, mean, sd, exact))
            stream.flush()
            print(
                f"N={n_sites:4d} | sampled={mean:.6f} ± {sd:.6f} "
                f"| exact={exact:.6f}",
                flush=True,
            )

    sites, means, sds, exacts = zip(*summary)
    fig, ax = plt.subplots(figsize=(7, 4.5))
    ax.errorbar(sites, means, yerr=sds, fmt="o-", capsize=4, label="Sampled PC1")
    ax.plot(sites, exacts, "--", label="Exact covariance")
    ax.set(xlabel="Number of sites N", ylabel=r"Overlap $\alpha$")
    ax.grid(alpha=0.25)
    ax.legend()
    fig.tight_layout()
    fig.savefig(output / f"{mode}_overlap_vs_size.png", dpi=200)
    plt.close(fig)


if __name__ == "__main__":
    main()