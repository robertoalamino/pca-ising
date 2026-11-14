import argparse
import csv
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]


def colour_masks(size):
    # Colour each periodic axis, then combine the two axis colours.
    axis = np.arange(size) % 2
    n_colours = 2
    if size % 2:
        axis[-1] = 2
        n_colours = 3
    grid = (axis[:, None] + axis[None, :]) % n_colours
    return [grid == colour for colour in range(n_colours)]


def overlap_residual(size, rng, burn_in, n_samples, gap=5):
    spins = rng.choice(np.array([-1, 1], dtype=np.int8), (size, size))
    masks = colour_masks(size)
    acceptance = np.exp(-np.array([0, 0, 0, 4, 8]))  # beta = J = 1

    def sweep():
        for mask in masks:
            neighbours = (
                np.roll(spins, 1, axis=0) + np.roll(spins, -1, axis=0)
                + np.roll(spins, 1, axis=1) + np.roll(spins, -1, axis=1)
            )
            delta = 2 * spins[mask] * neighbours[mask]
            probability = acceptance[(delta + 8) // 4]
            spins[mask] *= np.where(rng.random(delta.size) < probability, -1, 1)

    for _ in range(burn_in):
        sweep()

    samples = np.empty((n_samples, size * size), dtype=float)
    for i in range(n_samples):
        for _ in range(gap):
            sweep()
        samples[i] = spins.ravel()

    # Spin-flip symmetrisation gives a zero-mean covariance proportional
    # to X.T @ X. Use power iteration without forming that matrix.
    uniform = np.full(size * size, 1.0 / size)
    pc1 = uniform.copy()
    for _ in range(50):
        updated = samples.T @ (samples @ pc1)
        updated /= np.linalg.norm(updated)
        change = np.linalg.norm(updated - pc1)
        pc1 = updated
        if change < 1e-10:
            break
    if pc1 @ uniform < 0:
        pc1 = -pc1

    # For unit vectors this equals 1 - overlap, without subtracting
    # two numbers close to one.
    return 0.5 * np.sum((pc1 - uniform) ** 2)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--full", action="store_true")
    args = parser.parse_args()

    sizes = range(3, 40) if args.full else (3, 4)
    repetitions, burn_in, n_samples = (20, 5000, 500) if args.full else (2, 100, 100)
    mode = "full" if args.full else "smoke"
    output = ROOT / "outputs" / "square_lattice_size_sweep"
    output.mkdir(parents=True, exist_ok=True)
    path = output / f"{mode}_runs.csv"

    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(("L", "N", "repetition", "overlap_residual"))
        for size in sizes:
            residuals = []
            for repetition in range(repetitions):
                rng = np.random.default_rng([123, size, repetition])
                residual = overlap_residual(size, rng, burn_in, n_samples)
                residuals.append(residual)
                writer.writerow((size, size**2, repetition, residual))
            stream.flush()
            print(f"L={size:2d} | mean residual={np.mean(residuals):.6e}", flush=True)

    print(f"Wrote {path}")


if __name__ == "__main__":
    main()