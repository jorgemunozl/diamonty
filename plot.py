import matplotlib.pyplot as plt
import numpy as np


class ConvergencePlot:
    """Single-panel scatter plot of |ΔE| vs convergence parameter."""

    def __init__(self, x, xlabel, criterion, task_label):
        self.x = x
        self.xlabel = xlabel
        self.criterion = criterion
        self.task_label = task_label

        # ΔE between consecutive steps (meV)
        self.de_meV = np.abs(np.diff(self.x)) * 1000
        self.x_de = self.x[1:]

    def set_energies(self, energies):
        self.de_meV = np.abs(np.diff(energies)) * 1000

    def plot(self, filename):
        fig, ax = plt.subplots(figsize=(8, 5))
        # fig.suptitle(self.task_label, fontweight="bold")

        ax.plot(self.x_de, self.de_meV, "o-", color="steelblue", lw=1.8, ms=5, zorder=3)
        ax.axhline(
            self.criterion,
            color="#e74c3c",
            ls="--",
            lw=1.5,
            label=f"Criterion: {self.criterion:.0f} meV",
        )

        ax.set_xlabel(self.xlabel)
        ax.set_ylabel("$|\\Delta E|$  (meV)")
        ax.grid(True, ls="--", alpha=0.4)
        ax.legend(fontsize=9)

        plt.tight_layout()
        plt.savefig(filename, dpi=150, bbox_inches="tight")
        print(f"Saved: {filename}")
        plt.close(fig)


class PlotRunner:
    """Loads data and runs all convergence plots."""

    def __init__(self, data_path="analysis_data.npz"):
        self.data = np.load(data_path)

    # ── Task 1 ────────────────────────────────────────────────
    def task1_encut(self):
        cp = ConvergencePlot(
            x=self.data["cutoff_encut"],
            xlabel="ENCUT  (eV)",
            criterion=4.0,
            task_label="Task 1 – ENCUT Convergence (PBE)",
        )
        cp.set_energies(self.data["cutoff_e0"])
        cp.plot("encut_convergence.png")

    # ── Task 2 ────────────────────────────────────────────────
    def task2_kpoints(self):
        cp = ConvergencePlot(
            x=self.data["kpoint_mesh"],
            xlabel=r"k-mesh  (1/$\mathring{{\mathrm{A}}}^3)$",
            criterion=1.0,
            task_label=None,
        )
        cp.set_energies(self.data["kpoint_e0"])
        cp.plot("kpoint_convergence.png")

    # ── Run all ───────────────────────────────────────────────
    def run_all(self):
        self.task1_encut()
        self.task2_kpoints()


if __name__ == "__main__":
    runner = PlotRunner(data_path="analysis_data.npz")
    runner.run_all()
