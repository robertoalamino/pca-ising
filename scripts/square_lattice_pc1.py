from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from square_lattice_size_sweep import colour_masks

ROOT = Path(__file__).resolve().parents[1]
L = 20
BETA = 1.0
J = 1.0
BURN_IN = 100_000
N_SAMPLES = 5_000
GAP = 10
SEED = 2026


def metropolis_samples(size, beta, burn_in, n_samples, gap, seed, *, lazy=False):
    """Sample a periodic square lattice at zero field."""
    if size < 3:
        raise ValueError("Periodic lattice size must be at least 3")

    rng = np.random.default_rng(seed)
    spins = rng.choice(np.array([-1, 1], dtype=np.int8), size=(size, size))

    colours = colour_masks(size)
    # Possible ΔE / J values are -8, -4, 0, 4, and 8.
    acceptance = np.exp(-beta * J * np.array([0, 0, 0, 4, 8]))
    samples = np.empty((n_samples, size * size), dtype=np.int8)

    def sweep():
        for mask in colours:
            neighbours = (
                np.roll(spins, 1, axis=0) + np.roll(spins, -1, axis=0)
                + np.roll(spins, 1, axis=1) + np.roll(spins, -1, axis=1)
            )
            delta = 2 * spins[mask] * neighbours[mask]
            probability = acceptance[(delta + 8) // 4]
            if lazy:
                probability = 0.5 * probability
            spins[mask] *= np.where(rng.random(delta.size) < probability, -1, 1)

    for _ in range(burn_in):
        sweep()
    for i in range(n_samples):
        for _ in range(gap):
            sweep()
        samples[i] = spins.ravel()

    return samples


def main():
    samples = metropolis_samples(L, BETA, BURN_IN, N_SAMPLES, GAP, SEED)
    magnetisation = samples.mean(axis=1)

    # At zero field, each configuration and its global flip have equal
    # Gibbs weight. The paired sample has zero mean and second moment X.T @ X / n.
    x = samples.astype(np.float64)
    covariance = x.T @ x / N_SAMPLES
    _, vectors = np.linalg.eigh(covariance)
    pc1 = vectors[:, -1]
    if pc1.sum() < 0:
        pc1 = -pc1

    print(f"L={L} | beta={BETA} | samples={N_SAMPLES} | seed={SEED}")
    print(f"Raw chain mean |m|: {np.abs(magnetisation).mean():.6f}")
    print(f"Raw chain m range: [{magnetisation.min():.6f}, {magnetisation.max():.6f}]")
    print(f"PC1 overlap with uniform direction: {pc1.mean() * L:.9f}")
    print(f"Variance across PC1 components: {pc1.var():.3e}")

    output = ROOT / "outputs" / "square_lattice_pc1"
    output.mkdir(parents=True, exist_ok=True)
    grid = pc1.reshape(L, L)
    np.savetxt(output / "pc1_grid.csv", grid, delimiter=",")

    fig, ax = plt.subplots(figsize=(5, 4.5))
    image = ax.imshow(grid, cmap="viridis_r", interpolation="nearest")
    fig.colorbar(image, ax=ax, label="PC1 component")
    ax.set(xlabel="Column", ylabel="Row", title=r"Square lattice PC1 ($\beta=1$)")
    fig.tight_layout()
    fig.savefig(output / "pc1_square_lattice.png", dpi=200)
    plt.close(fig)


if __name__ == "__main__":
    main()