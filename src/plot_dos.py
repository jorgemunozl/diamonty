"""Generate DOS plots for slides — PBE and HSE06."""

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from constants import HSE06_DOS_DIR, PBE_DOS_DIR
from utils import read_doscar

configs = [
    ("PBE", PBE_DOS_DIR, "dos_pbe.png"),
    ("HSE06", HSE06_DOS_DIR, "dos_hse06.png"),
]

for functional, dos_dir, filename in configs:
    print(f"\n--- Density of States ({functional}) ---")
    energy, dos_total, e_fermi = read_doscar(dos_dir / "DOSCAR")

    energy_shifted = energy - e_fermi

    print(f"  E_Fermi = {e_fermi:.4f} eV")
    print(f"  Energy range: [{energy_shifted[0]:.2f}, {energy_shifted[-1]:.2f}] eV")

    # Band gap from DOS
    mask = (energy > e_fermi - 15) & (energy < e_fermi + 15)
    e_win = energy[mask]
    dos_win = dos_total[mask]

    below = e_win < e_fermi
    above = e_win > e_fermi

    nb = np.where(dos_win[below] > 1e-6)[0]
    vbm = float(e_win[below][nb[-1]] - e_fermi) if len(nb) > 0 else 0.0

    na = np.where(dos_win[above] > 1e-6)[0]
    cbm = float(e_win[above][na[0]] - e_fermi) if len(na) > 0 else 0.0

    band_gap = cbm - vbm
    print(f"  Band gap from DOS: {band_gap:.3f} eV")
    print(f"  VBM: {vbm:.3f} eV, CBM: {cbm:.3f} eV (rel. to E_F)")

    # ── Plot ──
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.plot(energy_shifted, dos_total, color="steelblue", linewidth=1.2)
    ax.axvline(x=0, color="gray", linestyle="--", linewidth=0.8, label="$E_F = 0$")
    ax.axvspan(
        vbm, cbm, color="lightgreen", alpha=0.25, label=f"Gap = {band_gap:.2f} eV"
    )

    ax.set_xlabel("$E - E_F$ (eV)", fontsize=11)
    ax.set_ylabel("DOS (states/eV)", fontsize=11)
    ax.set_title(f"DOS — Diamond ({functional})", fontsize=12)
    ax.legend(fontsize=9)
    ax.set_xlim(energy_shifted[0], energy_shifted[-1])

    fig.tight_layout()
    path = f"slides/img/{filename}"
    fig.savefig(path, dpi=150)
    plt.close(fig)
    print(f"  Saved {path}")
