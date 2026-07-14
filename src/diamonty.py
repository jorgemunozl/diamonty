from dataclasses import dataclass, field
from pathlib import Path
from typing import Literal

import matplotlib.pyplot as plt
import numpy as np
from pymatgen.electronic_structure.plotter import BSPlotter
from pymatgen.io.vasp.outputs import Vasprun

from constants import (
    HSE06_BAND,
    HSE06_DOS_DIR,
    HSE06_RELAX,
    MU_C,
    MU_N,
    PBE_BAND,
    PBE_DOS_DIR,
    PBE_RELAX,
    PLOTS_SLIDE,
    SUPERCELL_CONV,
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
    band_path: Path = field(
        default=PBE_BAND,
        metadata={"description": "Path to band structure vasprun.xml"},
    )

    def __post_init__(self):
        if self.functional == "HSE06":
            self.pseudo_relax_dir = HSE06_RELAX
            self.dos_dir = HSE06_DOS_DIR
            self.band_path = HSE06_BAND


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
        """Returns the cutoff energy convergence data."""
        cutoff_data = []
        for d in self.config.cutoff_dirs:
            encut = int(d.name)
            e0 = read_oszicar_energy(d / "OSZICAR")
            if e0 is not None:
                cutoff_data.append((encut, e0))
        return cutoff_data

    def kpoint(self):
        """Returns the k-point mesh convergence data."""
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
        """Read DOSCAR, compute band gap from DOS."""
        print(f"\n--- Density of States ({self.config.functional}) ---")
        energy, dos_total, e_fermi = read_doscar(self.config.dos_dir / "DOSCAR")
        energy_shifted = energy - e_fermi
        print(f"  E_Fermi = {e_fermi:.4f} eV")
        print(f"  Energy range: [{energy_shifted[0]:.2f}, {energy_shifted[-1]:.2f}] eV")
        print(f"  NEDOS = {len(energy)}")

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
        """Plot DOS with gap shading."""
        fig, ax = plt.subplots(figsize=(8, 5))
        ax.plot(energy_shifted, dos_total, color="steelblue", linewidth=1.5)
        ax.axvline(x=0, color="gray", linestyle="--", linewidth=0.8, label="$E_F$ = 0")
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
        """Extract lattice parameters from CONTCAR."""
        print(f"  Pseudo relaxation directory: {self.config.pseudo_relax_dir}")
        with open(self.config.pseudo_relax_dir) as f:
            contcar = f.readlines()
        scale = float(contcar[1].strip())
        a1 = scale * np.array([float(x) for x in contcar[2].split()])
        a2 = scale * np.array([float(x) for x in contcar[3].split()])
        a3 = scale * np.array([float(x) for x in contcar[4].split()])
        lat_const = np.linalg.norm(a1) * np.sqrt(2)
        vol = np.abs(np.dot(a1, np.cross(a2, a3)))

        print(f"\n--- Lattice Parameters ---")
        print(f"Pseudopotential : {self.config.functional}")
        print(f"  Primitive cell vectors (Angstrom):")
        print(f"    a1 = {a1}")
        print(f"    a2 = {a2}")
        print(f"    a3 = {a3}")
        print(f"  Primitive cell volume: {vol:.4f} Å³")
        print(f"  Cubic lattice constant a = {lat_const:.4f} Å")
        print(f"  Lattice constant: {lat_const:.6f} Å")
        print(f"  Volume: {vol:.6f} Å³")
        return lat_const, vol

    def band_structure(self):
        """
        Extract band structure from vasprun.xml using pymatgen.

        Returns
        -------
        bs : BandStructureSymmLine
        gap : dict  (keys: 'direct', 'transition', 'energy')
        """
        print(f"\n--- Band Structure ({self.config.functional}) ---")
        v = Vasprun(str(self.config.band_path), parse_potcar_file=False)
        bs = v.get_band_structure()
        gap = bs.get_band_gap()
        print(
            f"  Gap: {gap['energy']:.3f} eV ({'direct' if gap['direct'] else 'indirect'})"
        )
        print(f"  Transition: {gap['transition']}")
        print(f"  VBM: {bs.get_vbm()['energy']:.3f} eV")
        print(f"  CBM: {bs.get_cbm()['energy']:.3f} eV")
        return bs, gap

    def plot_band_structure(self, bs, ylim=(-15, 20)):
        """
        Plot band structure along high-symmetry k-path.

        Uses pymatgen's BSPlotter for k-path labels, then
        overrides font sizes and shades the band gap.
        """
        plotter = BSPlotter(bs)
        plotter.get_plot(ylim=ylim)
        fig = plt.gcf()
        fig.set_size_inches(6, 4.5)
        ax = plt.gca()
        ax.set_title("")
        ax.set_xlabel("k-path", fontsize=10)
        ax.set_ylabel(r"$E - E_F$ (eV)", fontsize=10)
        ax.tick_params(labelsize=8)
        for label in ax.get_xticklabels():
            label.set_fontsize(8)
        ax.legend(fontsize=8, loc="upper right")

        fig.tight_layout(pad=0.5)
        fig.savefig(PLOTS_SLIDE / f"band_{self.config.functional}.pdf")
        return fig

    def formation_energy(self):
        """
        Compute NV center formation energy vs supercell size.

        E_form = E_def - E_perf - μ_C + μ_N
        """
        base = SUPERCELL_CONV
        print(f"\n--- NV Center Formation Energy ---")
        print(f"  μ_C = {MU_C:.6f} eV")
        print(f"  μ_N = {MU_N:.6f} eV")
        print(
            f"\n{'Supercell':<12} {'Atoms':<8} {'E_def (eV)':<18} "
            f"{'E_perf (eV)':<18} {'E_form (eV)':<14} {'ΔE (meV)':<12}"
        )
        print("-" * 78)

        results = []
        prev_form = None
        for n in range(1, 8):
            v_def = Vasprun(
                str(base / "defect" / str(n) / "vasprun.xml"), parse_potcar_file=False
            )
            v_perf = Vasprun(
                str(base / "pristine" / str(n) / "vasprun.xml"), parse_potcar_file=False
            )
            n_atoms = v_perf.final_structure.num_sites
            e_def = v_def.final_energy
            e_perf = v_perf.final_energy
            e_form = e_def - e_perf - MU_C + MU_N
            de_str = (
                f"{(e_form - prev_form) * 1000:.1f}" if prev_form is not None else "—"
            )
            print(
                f"  {n}x{n}x{n:<8} {n_atoms:<8} {e_def:<18.6f} {e_perf:<18.6f} "
                f"{e_form:<14.6f} {de_str:<12}"
            )
            results.append((f"{n}x{n}x{n}", n_atoms, e_def, e_perf, e_form))
            prev_form = e_form
        return results

    def plot_formation_energy(self, results):
        """Plot formation energy vs supercell size."""
        n_atoms = [r[1] for r in results]
        e_form = [r[4] for r in results]
        fig, ax = plt.subplots(figsize=(7, 4.5))
        ax.plot(n_atoms, e_form, "o-", color="steelblue", linewidth=1.5, markersize=6)
        ax.axhline(
            y=e_form[-1],
            color="gray",
            linestyle="--",
            linewidth=0.8,
            label=f"Converged: {e_form[-1]:.3f} eV",
        )
        ax.set_xlabel("Number of atoms", fontsize=11)
        ax.set_ylabel("Formation energy (eV)", fontsize=11)
        ax.set_title("NV Center Formation Energy — Supercell Convergence", fontsize=12)
        ax.legend(fontsize=15)
        fig.tight_layout()
        # plt.show()
        plt.savefig(PLOTS_SLIDE / "formation_energy.pdf", format="pdf")
        return fig


if __name__ == "__main__":
    config = Config(functional="HSE06")
    # config = Config(functional="PBE")
    diamonty = Diamonty(config)

    # Formation energy
    # results = diamonty.formation_energy()
    # diamonty.plot_formation_energy(results)

    # Band structure
    bs, gap = diamonty.band_structure()
    diamonty.plot_band_structure(bs)
