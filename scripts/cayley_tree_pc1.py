import csv
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

DEPTH = 7
BETA = 1.0
J = 1.0
N_SAMPLES = 20_000
SEED = 2026


def main():
    rng = np.random.default_rng(SEED)
    n_sites = 2 ** (DEPTH + 1) - 1
    parent = (np.arange(1, n_sites) - 1) // 2
    generation = np.zeros(n_sites, dtype=int)
    distance = np.zeros((n_sites, n_sites), dtype=int)

    correlation = np.tanh(BETA * J)
    p_same = (1 + correlation) / 2
    spins = np.empty((N_SAMPLES, n_sites), dtype=np.int8)
    spins[:, 0] = rng.choice([-1, 1], size=N_SAMPLES)

    for node, p in enumerate(parent, start=1):
        generation[node] = generation[p] + 1
        spins[:, node] = spins[:, p] * np.where(
            rng.random(N_SAMPLES) < p_same, 1, -1
        )

        # Every path from this node to an earlier node passes through its parent.
        distance[node, :node] = distance[p, :node] + 1
        distance[:node, node] = distance[node, :node]

    sampled_pc1 = leading_direction(np.cov(spins, rowvar=False))
    exact_pc1 = leading_direction(correlation**distance)

    print(f"Sites: {n_sites} | Samples: {N_SAMPLES} | Seed: {SEED}")
    print(f"Sampled overlap: {sampled_pc1.sum() / np.sqrt(n_sites):.6f}")
    print(f"Exact overlap:   {exact_pc1.sum() / np.sqrt(n_sites):.6f}")
    rows = []
    print("\nGeneration   Sites   Sampled mean   Within-layer SD   Exact mean")
    for level in range(DEPTH + 1):
        mask = generation == level
        values = sampled_pc1[mask]
        count = int(mask.sum())
        mean = float(values.mean())
        sd = float(values.std())
        exact_mean = float(exact_pc1[mask].mean())
        rows.append((level, count, mean, sd, exact_mean))

        print(
            f"{level:10d} {count:7d} {mean:14.6f} "
            f"{sd:17.6f} {exact_mean:12.6f}"
        )

    output_dir = Path("outputs") / "cayley_tree_pc1"
    output_dir.mkdir(parents=True, exist_ok=True)

    csv_path = output_dir / "pc1_by_generation.csv"
    with csv_path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(
            ("generation", "sites", "sampled_mean", "within_layer_sd", "exact_mean")
        )
        writer.writerows(rows)

    levels, _, means, sds, exact_means = zip(*rows)
    fig, ax = plt.subplots(figsize=(7, 4))
    ax.errorbar(
        levels, means, yerr=sds, fmt="o-", capsize=4,
        markerfacecolor="white", label="Sampled PC1",
    )
    ax.plot(levels, exact_means, "--", color="0.35", label="Exact covariance")
    ax.set(xlabel="Generation", ylabel="Mean PC1", xticks=levels)
    ax.grid(alpha=0.3)
    ax.legend()
    fig.tight_layout()

    plot_path = output_dir / "pc1_by_generation.png"
    fig.savefig(plot_path, dpi=200)
    plt.close(fig)
    print(f"\nWrote {csv_path} and {plot_path}")

    radial_path = output_dir / "pc1_radial_tree.png"
    plot_radial_tree(parent, generation, sampled_pc1, radial_path)
    print(f"Wrote {radial_path}")

def plot_radial_tree(parent, generation, pc1, output_path):
    # Put each generation on its own ring. Align each parent
    # with the centre of the angular sector occupied by its descendants.
    x = np.zeros(len(pc1))
    y = np.zeros(len(pc1))
    for level in range(1, int(generation.max()) + 1):
        nodes = np.flatnonzero(generation == level)
        angles = 2 * np.pi * (np.arange(len(nodes)) + 0.5) / len(nodes)
        x[nodes] = level * np.cos(angles)
        y[nodes] = level * np.sin(angles)

    fig, ax = plt.subplots(figsize=(7, 7))
    for level in range(1, int(generation.max()) + 1):
        ax.add_patch(plt.Circle(
            (0, 0), level, fill=False, linestyle="--",
            linewidth=0.6, color="0.8", zorder=0,
        ))

    for node, p in enumerate(parent, start=1):
        ax.plot([x[p], x[node]], [y[p], y[node]],
                color="0.7", linewidth=0.5, zorder=1)

    colours = dict(cmap="plasma_r", vmin=pc1.min(), vmax=pc1.max())
    points = ax.scatter(x[1:], y[1:], c=pc1[1:], s=14, zorder=2, **colours)
    ax.scatter(x[:1], y[:1], c=pc1[:1], s=55, marker="D",
               edgecolors="black", linewidths=0.5, zorder=3, **colours)
    fig.colorbar(points, ax=ax, shrink=0.75, label="PC1 component")
    ax.set_aspect("equal")
    ax.set_axis_off()
    fig.savefig(output_path, dpi=200, bbox_inches="tight")
    plt.close(fig)

def leading_direction(covariance):
    _, vectors = np.linalg.eigh(covariance)
    pc1 = vectors[:, -1]
    return pc1 if pc1.sum() >= 0 else -pc1


if __name__ == "__main__":
    main()