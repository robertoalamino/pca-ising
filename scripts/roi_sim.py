import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from roi_mix import ROOT, ROI_DEPTH, exact_ratios


def main():
    buffers = np.arange(1, 51)
    ratios = np.array([exact_ratios(ROI_DEPTH + b) for b in buffers])

    output_dir = ROOT / "outputs" / "roi_sim"
    output_dir.mkdir(parents=True, exist_ok=True)

    csv_path = output_dir / "roi_exact_ratios.csv"
    np.savetxt(
        csv_path,
        np.column_stack((buffers, ratios)),
        delimiter=",",
        header="buffer,level_1,level_2,level_3,level_4",
        comments="",
        fmt=["%d"] + ["%.12g"] * ROI_DEPTH,
    )

    fig, ax = plt.subplots(figsize=(8, 4.5))
    for level in range(1, ROI_DEPTH + 1):
        ax.plot(
            buffers, ratios[:, level - 1],
            marker="o", markersize=3, label=fr"$\ell={level}$",
        )

    ax.axhline(1.0, color="0.5", linestyle="--", linewidth=1)
    ax.set(xlabel=r"Buffer size $m_b$", ylabel=r"$R_\ell=v_\ell/v_0$")
    ax.grid(alpha=0.25)
    ax.legend(loc="center left", bbox_to_anchor=(1, 0.5))
    fig.tight_layout()

    plot_path = output_dir / "roi_sim.png"
    fig.savefig(plot_path, dpi=200)
    plt.close(fig)

    print(f"Wrote {csv_path} and {plot_path}")
    print("Ratios at buffer 50:", np.round(ratios[-1], 6))


if __name__ == "__main__":
    main()