from dataclasses import dataclass, field
from pathlib import Path
from typing import Literal

import matplotlib.pyplot as plt
import numpy as np
from pymatgen.electronic_structure.core import Spin
from pymatgen.electronic_structure.plotter import BSPlotter, DosPlotter
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
    POINT_DEFECT,
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
        """Read DOS from vasprun.xml using pymatgen."""
        print(f"\n--- Density of States ({self.config.functional}) ---")
        v = Vasprun(str(self.config.dos_dir / "vasprun.xml"), parse_potcar_file=False)
        cdos = v.complete_dos
        gap = cdos.get_gap()
        print(f"  E_Fermi = {cdos.efermi:.4f} eV")
        print(f"  Gap: {gap:.3f} eV (from DOS)")
        cbm, vbm = cdos.get_cbm_vbm()
        print(f"  VBM: {vbm:.4f}, CBM: {cbm:.4f} eV")
        return cdos

    def plot_dos(self, cdos=None, xlim=(-20, 20)):
        """Plot DOS with pymatgen DosPlotter, saves to slides/img/."""
        if cdos is None:
            cdos = self.dos()

        plotter = DosPlotter(zero_at_efermi=True)
        plotter.add_dos_dict(cdos.get_element_dos())
        plotter.get_plot(xlim=xlim)
        fig = plt.gcf()
        fig.set_size_inches(6, 4)
        ax = plt.gca()
        ax.set_title(f"DOS — Diamond ({self.config.functional})", fontsize=6)
        ax.tick_params(labelsize=8)
        ax.legend(fontsize=8, loc="upper right")
        fig.tight_layout(pad=0.5)
        fig.savefig(f"slides/img/dos_{self.config.functional.lower()}.pdf")
        return fig

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

        # Shade the band gap
        efermi = bs.efermi
        vbm = bs.get_vbm()["energy"] - efermi
        cbm = bs.get_cbm()["energy"] - efermi
        ax.axhspan(
            vbm,
            cbm,
            color="lightgreen",
            alpha=0.2,
            zorder=0,
            label=f"$E_g = {cbm - vbm:.2f}$ eV",
        )

        ax.set_xlabel("k-path", fontsize=10)
        ax.set_ylabel(r"$E - E_F$ (eV)", fontsize=10)
        ax.tick_params(labelsize=8)
        for label in ax.get_xticklabels():
            label.set_fontsize(8)
        ax.legend(fontsize=8, loc="upper right")

        fig.tight_layout(pad=0.5)
        fig.savefig(PLOTS_SLIDE / f"band_{self.config.functional}.pdf")

    def defect_levels(self):
        """
        Extract single-particle defect levels for all charge states.

        Returns
        -------
        vbm, cbm : float  (perfect-cell band edges in eV)
        levels : dict  {charge: {Spin.up: [(E_rel, occ), ...], Spin.down: [...]}}
        """
        base = POINT_DEFECT
        v_perf = Vasprun(str(base / "perfect" / "vasprun.xml"), parse_potcar_file=False)
        all_e = np.array(list(v_perf.eigenvalues.values()))
        vbm = float(all_e[all_e < v_perf.efermi].max())
        cbm = float(all_e[all_e >= v_perf.efermi].min())

        charge_map = {
            "N_C-V_C_-3": -3,
            "N_C-V_C_-2": -2,
            "N_C-V_C_-1": -1,
            "N_C-V_C_0": 0,
            "N_C-V_C_1": 1,
            "N_C-V_C_2": 2,
        }

        print(f"\n--- Defect Levels (NV Center) ---")
        print(
            f"  Perfect cell: VBM = {vbm:.3f} eV, CBM = {cbm:.3f} eV, gap = {cbm - vbm:.3f} eV"
        )

        levels = {}
        for folder, q in charge_map.items():
            v = Vasprun(str(base / folder / "vasprun.xml"), parse_potcar_file=False)
            efermi = v.efermi
            levels[q] = {}
            for spin, arr in v.eigenvalues.items():
                flat = arr.flatten()
                mask = (flat > vbm) & (flat < cbm)
                gap_e_rel = flat[mask] - vbm
                occ = flat[mask] < efermi
                levels[q][spin] = list(zip(gap_e_rel, occ))
                n_occ = sum(occ)
                print(
                    f"  q={q:+d}  {spin.name}: {len(gap_e_rel)} gap states ({n_occ} occupied)"
                )

        return vbm, cbm, levels

    def plot_defect_levels(self, vbm, cbm, levels):
        """
        Plot single-particle defect level diagram.
        """
        from matplotlib.lines import Line2D

        fig, ax = plt.subplots(figsize=(8, 5))
        gap = cbm - vbm

        # Band edges
        ax.axhline(y=0, color="black", linewidth=1.2)
        ax.axhline(y=gap, color="black", linewidth=1.2)
        ax.text(-3.4, -0.15, "VBM", fontsize=9, va="top")
        ax.text(-3.4, gap + 0.1, "CBM", fontsize=9, va="bottom")

        color_up = "steelblue"
        color_dn = "crimson"

        for q in sorted(levels.keys()):
            for spin, states in levels[q].items():
                is_up = spin == Spin.up
                color = color_up if is_up else color_dn
                x = q + (0.15 if is_up else -0.15)
                for e_rel, occ in states:
                    ax.plot(
                        x,
                        e_rel,
                        "o" if occ else "^",
                        color=color,
                        markerfacecolor=color if occ else "none",
                        markersize=8,
                        markeredgewidth=1.2,
                    )

        ax.set_xticks(sorted(levels.keys()))
        ax.set_xticklabels([f"{q:+d}" for q in sorted(levels.keys())], fontsize=10)
        ax.set_xlabel("Charge state", fontsize=11)
        ax.set_ylabel("E − VBM (eV)", fontsize=11)
        ax.set_title("NV Center — Single-Particle Defect Levels", fontsize=12)
        ax.set_ylim(-1, gap + 1)

        legend = [
            Line2D(
                [0],
                [0],
                marker="o",
                color="w",
                markerfacecolor=color_up,
                markersize=8,
                label="↑ occupied",
            ),
            Line2D(
                [0],
                [0],
                marker="^",
                color="w",
                markerfacecolor="none",
                markeredgecolor=color_up,
                markersize=8,
                markeredgewidth=1.2,
                label="↑ empty",
            ),
            Line2D(
                [0],
                [0],
                marker="o",
                color="w",
                markerfacecolor=color_dn,
                markersize=8,
                label="↓ occupied",
            ),
            Line2D(
                [0],
                [0],
                marker="^",
                color="w",
                markerfacecolor="none",
                markeredgecolor=color_dn,
                markersize=8,
                markeredgewidth=1.2,
                label="↓ empty",
            ),
        ]
        ax.legend(handles=legend, fontsize=7, loc="upper right", ncol=2)

        fig.tight_layout()
        fig.savefig("slides/img/defect_levels.png", dpi=150)
        return fig
    def plot_defect_eigenvalues(self, charge_state="N_C-V_C_0", ylim=None):
        """
        LSPD-style eigenvalue plot for a single charge state at Gamma.

        Spin-up on left, spin-down on right. Colors:
        blue = occupied (> 0.9), red = empty (< 0.1), green = partial.
        """
        from pymatgen.electronic_structure.core import Spin

        base = POINT_DEFECT
        v = Vasprun(str(base / charge_state / "vasprun.xml"), parse_potcar_file=False)
        eigen = v.eigenvalues
        efermi = v.efermi

        # Get VBM/CBM from perfect
        v_perf = Vasprun(str(base / "perfect" / "vasprun.xml"), parse_potcar_file=False)
        all_e = np.array(list(v_perf.eigenvalues.values()))
        perf_vbm = all_e[all_e < v_perf.efermi].max()
        perf_cbm = all_e[all_e >= v_perf.efermi].min()

        if ylim is None:
            ylim = (perf_vbm - 2, perf_cbm + 2)

        fig, (ax_up, ax_dn) = plt.subplots(1, 2, figsize=(5, 6), sharey=True)

        for ax, spin in [(ax_up, Spin.up), (ax_dn, Spin.down)]:
            arr = eigen[spin]  # shape: (nkpts, nbands, 2) -> [energy, occ]
            energies = arr[:, :, 0].flatten()
            occupations = arr[:, :, 1].flatten()

            for e, occ in zip(energies, occupations):
                if occ > 0.9:
                    c, f = "steelblue", "steelblue"
                elif occ < 0.1:
                    c, f = "crimson", "none"
                else:
                    c, f = "green", "green"
                ax.plot(0, e, "o", color=c, markerfacecolor=f,
                        markersize=10, markeredgewidth=1.2, zorder=2)

            # VBM/CBM lines
            ax.axhline(y=perf_vbm, color="black", linewidth=0.8, linestyle="--")
            ax.axhline(y=perf_cbm, color="black", linewidth=0.8, linestyle="--")
            ax.set_title(f"spin {spin.name}", fontsize=11)
            ax.set_xticks([])
            ax.set_ylim(ylim)
        fig.suptitle(f"NV Center — KS Eigenvalues ({charge_state})", fontsize=12)
        fig.tight_layout()
        fig.savefig(f"slides/img/eigen_{charge_state}.pdf")
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
    # bs, gap = diamonty.band_structure()
    # diamonty.plot_band_structure(bs)

    # Defect levels
    # vbm, cbm, levels = diamonty.defect_levels()
    # diamonty.plot_defect_levels(levels)

    # DOS
    dos = diamonty.dos()
    diamonty.plot_dos(dos)

    # LDOS
    # ldos = diamonty.ldos()
    # diamonty.plot_ldos(ldos)
