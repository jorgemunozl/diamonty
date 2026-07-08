from dataclasses import dataclass, field
from pathlib import Path
from typing import Literal, Optional

import matplotlib.pyplot as plt
import numpy as np

from constants import (
    HSE06_DOS_DIR,
    HSE06_RELAX,
    PBE_DOS_DIR,
    PBE_RELAX,
    cutoff_dirs_pbe,
    kdensity_dirs_pbe,
)
from utils import read_doscar, read_oszicar_energy


@dataclass
class Config:
    functional: Literal["PBE", "HSE06"] = field(
        default="PBE",
        metadata={"description": "Functional to use for the calculation"},
    )
    cutoff_dirs: list[Path] = field(
        default_factory=lambda: cutoff_dirs_pbe,
        metadata={"description": "Directories containing cutoff energy data"},
    )
    kdensity_dirs: list[Path] = field(
        default_factory=lambda: kdensity_dirs_pbe,
        metadata={"description": "Directories containing k-point density data"},
    )
    pseudo_relax_dir: Path = field(
        default=PBE_RELAX,
        metadata={"description": "Directory containing PBE relaxation results"},
    )
    dos_dir: Path = field(
        default=PBE_DOS_DIR,
        metadata={"description": "Directory containing DOSCAR data"},
    )

    def __post_init__(self):
        if self.functional == "HSE06":
            self.pseudo_relax_dir = HSE06_RELAX
            self.dos_dir = HSE06_DOS_DIR


class Printer:
    def __init__(self, config: Config):
        self.config = config

    def encut(self, cutoff_data):
        print(f"\n{'ENCUT (eV)':<15} {'E0 (eV)':<20} {'ΔE (meV)':<15} {'Converged?'}")
        print("-" * 65)
        prev_e = None
        converged_at = None
        for i, (encut, e0) in enumerate(cutoff_data):
            if prev_e is not None:
                de = abs(e0 - prev_e) * 1000
                conv = "YES" if de < 4.0 else "NO"
                if de < 4.0 and converged_at is None:
                    converged_at = encut
                print(f"  {encut:<13} {e0:<20.8f} {de:<15.3f} {conv}")
            else:
                print(f"  {encut:<13} {e0:<20.8f} {'--':<15} --")
            prev_e = e0

        rec_encut = converged_at + 50 if converged_at else 500
        print(f"\nConverged ENCUT >= {converged_at} eV (ΔE < 4 meV criterion)")
        print(f"Recommended ENCUT = {rec_encut} eV (with 50 eV safety margin)")

    def k_points(self, kpoint_data):
        print(
            f"\n{'Index':<8} {'k-mesh':<16} {'E0 (eV)':<20} {'ΔE (meV)':<15} {'Converged?'}"
        )
        print("-" * 72)

        prev_e = None
        converged_at_k = None
        for i, (idx, nk, e0) in enumerate(kpoint_data):
            if prev_e is not None:
                de = abs(e0 - prev_e) * 1000
                conv = "YES" if de < 1.0 else "NO"
                if de < 1.0 and converged_at_k is None:
                    converged_at_k = f"{nk}x{nk}x{nk}"
                print(f"  {idx:<6} {nk}x{nk}x{nk:<9} {e0:<20.8f} {de:<15.3f} {conv}")
            else:
                print(f"  {idx:<6} {nk}x{nk}x{nk:<9} {e0:<20.8f} {'--':<15} --")
            prev_e = e0

        print(f"\nConverged k-mesh >= {converged_at_k} (ΔE < 1 meV criterion)")

    def dos(self):
        """ """


class Diamonty:
    def __init__(self, config: Config):
        self.config = config

    def encut(self):
        """
        Returns the cutoff energy for the calculation.
        """
        cutoff_data = []
        for d in self.config.cutoff_dirs:
            encut = int(d.name)
            e0 = read_oszicar_energy(d / "OSZICAR")
            if e0 is not None:
                cutoff_data.append((encut, e0))
        return cutoff_data

    def kpoint(self):
        """
        Returns the k-point mesh for the calculation.
        """
        kpoint_data = []
        for d in self.config.kdensity_dirs:
            idx = int(d.name)
            with open(d / "KPOINTS") as f:
                kp_lines = f.readlines()
            mesh = [int(x) for x in kp_lines[3].split()]
            nk = mesh[0]

            e0 = read_oszicar_energy(d / "OSZICAR")
            if e0 is not None:
                kpoint_data.append((idx, nk, e0))
        return kpoint_data

    def dos(self):
        """
        Read DOSCAR, compute band gap, and plot DOS vs E - E_Fermi.
        """
        print(f"\n--- Density of States ({self.config.functional}) ---")
        energy, dos_total, e_fermi = read_doscar(self.config.dos_dir / "DOSCAR")

        # Shift energy so E_Fermi = 0
        energy_shifted = energy - e_fermi

        print(f"  E_Fermi = {e_fermi:.4f} eV")
        print(f"  Energy range: [{energy_shifted[0]:.2f}, {energy_shifted[-1]:.2f}] eV")
        print(f"  NEDOS = {len(energy)}")

        # Band gap from DOS (look within ±15 eV of E_F)
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

    def plot_dos(self, energy_shifted, dos_total, vbm, cbm, band_gap):
        fig, ax = plt.subplots(figsize=(8, 5))
        ax.plot(energy_shifted, dos_total, color="steelblue", linewidth=1.5)
        ax.axvline(x=0, color="gray", linestyle="--", linewidth=0.8, label="$E_F$ = 0")

        # Shade the band gap
        ax.axvspan(
            vbm, cbm, color="lightgreen", alpha=0.2, label=f"Gap = {band_gap:.2f} eV"
        )

        ax.set_xlabel("$E - E_F$ (eV)", fontsize=12)
        ax.set_ylabel("DOS (states/eV)", fontsize=12)
        ax.set_title(
            f"Density of States — Diamond ({self.config.functional})", fontsize=13
        )
        ax.legend(fontsize=10)
        ax.set_xlim(energy_shifted[0], energy_shifted[-1])

        fig.tight_layout()
        plt.show()

        return energy_shifted, dos_total, e_fermi, vbm, cbm, band_gap

    def extract_lattice(self):
        """
        Extracts the lattice parameters from the OUTCAR file.
        """
        print(f"  Pseudo relaxation directory: {self.config.pseudo_relax_dir}")
        with open(self.config.pseudo_relax_dir) as f:
            contcar = f.readlines()

        scale = float(contcar[1].strip())
        a1 = scale * np.array([float(x) for x in contcar[2].split()])
        a2 = scale * np.array([float(x) for x in contcar[3].split()])
        a3 = scale * np.array([float(x) for x in contcar[4].split()])

        lat_const_pbe = np.linalg.norm(a1) * np.sqrt(2)
        vol_pbe = np.abs(np.dot(a1, np.cross(a2, a3)))

        print(f"\n--- Lattice Parameters ---")
        print(f"Pseudopotential : {self.config.functional}")
        print(f"  Primitive cell vectors (Angstrom):")
        print(f"    a1 = {a1}")
        print(f"    a2 = {a2}")
        print(f"    a3 = {a3}")
        print(f"  Primitive cell volume: {vol_pbe:.4f} Å³")
        print(f"  Cubic lattice constant a = {lat_const_pbe:.4f} Å")
        print(f"  Lattice constant: {lat_const_pbe:.6f} Å")
        print(f"  Volume: {vol_pbe:.6f} Å³")
        return lat_const_pbe, vol_pbe


if __name__ == "__main__":
    config = Config(functional="PBE")
    # config = Config(functional="HSE06")
    diamonty = Diamonty(config)
    # diamonty.kpoint()
    diamonty.extract_lattice()
    # diamonty.dos()
