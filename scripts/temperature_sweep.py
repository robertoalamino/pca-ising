import numpy as np
import argparse
import csv

from pathlib import Path

from square_lattice_pc1 import metropolis_samples


def observables(samples):
    """Return PC1 overlap, relative spectral gap, and mean |m|."""
    x = samples.astype(np.float64)
    # Spin-flip symmetrisation gives zero mean and this covariance.
    covariance = x.T @ x / len(x)
    eigenvalues, eigenvectors = np.linalg.eigh(covariance)

    pc1 = eigenvectors[:, -1]
    overlap = abs(pc1.sum()) / np.sqrt(x.shape[1])
    relative_gap = (eigenvalues[-1] - eigenvalues[-2]) / eigenvalues[-1]
    mean_abs_m = np.abs(x.mean(axis=1)).mean()
    return np.array([overlap, relative_gap, mean_abs_m])

def tree_samples(beta, n_samples, seed, depth=7):
    """Draw independent zero-field Ising configurations on a binary tree."""
    rng = np.random.default_rng(seed)
    n_sites = 2 ** (depth + 1) - 1
    samples = np.empty((n_samples, n_sites), dtype=np.int8)
    samples[:, 0] = rng.choice([-1, 1], size=n_samples)
    p_same = (1 + np.tanh(beta)) / 2  # J = 1

    for level in range(1, depth + 1):
        first = 2**level - 1
        stop = 2 ** (level + 1) - 1
        parents = (np.arange(first, stop) - 1) // 2
        signs = np.where(
            rng.random((n_samples, stop - first)) < p_same, 1, -1
        )
        samples[:, first:stop] = samples[:, parents] * signs

    return samples

def run_full_lattice():
    root = Path(__file__).resolve().parents[1]
    output = root / "outputs" / "temperature_sweep"
    output.mkdir(parents=True, exist_ok=True)
    path = output / "square_lattice_runs.csv"

    betas = np.linspace(0.001, 1.0, 100)
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow((
            "beta_index", "beta", "repetition",
            "overlap", "relative_gap", "mean_abs_m",
        ))

        for index, beta in enumerate(betas):
            runs = []
            for repetition in range(20):
                samples = metropolis_samples(
                    size=15,
                    beta=beta,
                    burn_in=5_000,
                    n_samples=500,
                    gap=5,
                    seed=[2026, 0, index, repetition],
                    lazy=True,
                )
                values = observables(samples)
                runs.append(values)
                writer.writerow((index, beta, repetition, *values))

            stream.flush()
            overlap, gap, magnetisation = np.mean(runs, axis=0)
            print(
                f"{index + 1:3d}/100 | beta={beta:.4f} | "
                f"overlap={overlap:.4f} | gap={gap:.4f} | "
                f"mean |m|={magnetisation:.4f}",
                flush=True,
            )

    print(f"Wrote {path}")

def run_full_tree():
    root = Path(__file__).resolve().parents[1]
    output = root / "outputs" / "temperature_sweep"
    output.mkdir(parents=True, exist_ok=True)
    path = output / "cayley_tree_runs.csv"

    betas = np.linspace(0.001, 5.0, 100)
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow((
            "beta_index", "beta", "repetition",
            "overlap", "relative_gap", "mean_abs_m",
        ))

        for index, beta in enumerate(betas):
            runs = []
            for repetition in range(20):
                samples = tree_samples(
                    beta, 500, seed=[2026, 1, index, repetition], depth=7,
                )
                values = observables(samples)
                runs.append(values)
                writer.writerow((index, beta, repetition, *values))

            stream.flush()
            overlap, gap, magnetisation = np.mean(runs, axis=0)
            print(
                f"{index + 1:3d}/{len(betas)} | beta={beta:.4f} | "
                f"overlap={overlap:.4f} | gap={gap:.4f} | "
                f"mean |m|={magnetisation:.4f}",
                flush=True,
            )

    print(f"Wrote {path}")

def main():
    parser = argparse.ArgumentParser()
    modes = parser.add_mutually_exclusive_group()
    modes.add_argument("--full-lattice", action="store_true")
    modes.add_argument("--full-tree", action="store_true")
    args = parser.parse_args()

    if args.full_lattice:
        run_full_lattice()
        return
    if args.full_tree:
        run_full_tree()
        return
    print(" beta    overlap    relative gap    mean |m|")
    for index, beta in enumerate((0.001, 0.2, 0.44, 0.6, 1.0)):
        runs = [
            observables(metropolis_samples(
                size=15,
                beta=beta,
                burn_in=5_000,
                n_samples=500,
                gap=5,
                seed=[2026, index, repetition],
                lazy=True,
            ))
            for repetition in range(3)
        ]
        overlap, gap, magnetisation = np.mean(runs, axis=0)
        print(
            f"{beta:5.3f}    {overlap:7.4f}         "
            f"{gap:7.4f}     {magnetisation:7.4f}",
            flush=True,
        )

    print("\n255-site binary tree")
    print(" beta    overlap    relative gap    mean |m|")
    for index, beta in enumerate((0.001, 0.2, 0.44, 0.6, 1.0)):
        runs = [
            observables(tree_samples(beta, 500, [2026, 1, index, repetition]))
            for repetition in range(3)
        ]
        overlap, gap, magnetisation = np.mean(runs, axis=0)
        print(
            f"{beta:5.3f}    {overlap:7.4f}         "
            f"{gap:7.4f}     {magnetisation:7.4f}",
            flush=True,
        )

if __name__ == "__main__":
    main()