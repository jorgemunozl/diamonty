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
    ZPL_EXCITED_DIR,
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

    def plot_ldos(self, cdos=None, xlim=(-20, 20)):
        """Plot orbital-projected DOS (s, p, d) — saves to slides/img/."""
        if cdos is None:
            cdos = self.dos()

        spd = cdos.get_element_spd_dos("C")
        plotter = DosPlotter(zero_at_efermi=True)
        plotter.add_dos_dict(spd)
        plotter.get_plot(xlim=xlim)
        fig = plt.gcf()
        fig.set_size_inches(6, 4)
        ax = plt.gca()
        ax.set_title(f"LDOS — Diamond ({self.config.functional})", fontsize=11)
        ax.tick_params(labelsize=8)
        ax.legend(fontsize=8, loc="upper right")
        fig.tight_layout(pad=0.5)
        fig.savefig(f"slides/img/ldos_{self.config.functional.lower()}.pdf")
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

    def formation_energy_diagram(self):
        """
        Compute and plot the thermodynamic formation energy diagram
        for the NV center as a function of Fermi level.

        E^f(q, E_F) = E_def(q) - E_perf - Σ n_i μ_i + q·(E_VBM + E_F) + E_corr(q)

        For NV center: n_C = -1 (remove 1 C), n_N = +1 (add 1 N)
        So: -n_C·μ_C - n_N·μ_N = +μ_C - μ_N
        Actually: E_form = E_def - E_perf - (-1)*μ_C - (+1)*μ_N
                   = E_def - E_perf + μ_C - μ_N
        """
        base = POINT_DEFECT
        band_base = HSE06_BAND.parent.parent / "band"

        # ── 1. VBM from HSE06 band structure ──
        v_band = Vasprun(str(HSE06_BAND), parse_potcar_file=False)
        all_e = np.array(list(v_band.eigenvalues.values()))
        e_fermi_band = v_band.efermi
        vbm = float(all_e[all_e < e_fermi_band].max())
        cbm = float(all_e[all_e >= e_fermi_band].min())
        band_gap = cbm - vbm
        print(f"\n--- Formation Energy Diagram ---")
        print(f"  VBM = {vbm:.4f} eV")
        print(f"  CBM = {cbm:.4f} eV")
        print(f"  Band gap = {band_gap:.4f} eV")

        # ── 2. Chemical potentials from competing phases ──
        # μ_C from diamond: use HSE06 relax from 4.Properties_of_Diamond
        diamond_vasprun = (
            BASE.parent
            / "diamonty"
            / "4.Properties_of_Diamond"
            / "vasp"
            / "HSE06"
            / "relax"
            / "vasprun.xml"
        )
        v_diamond = Vasprun(str(diamond_vasprun), parse_potcar_file=False)
        e_diamond = v_diamond.final_energy
        n_c = sum(1 for site in v_diamond.final_structure if site.species_string == "C")
        mu_C = e_diamond / n_c
        print(f"  μ_C (from diamond, {n_c} C atoms) = {mu_C:.6f} eV/atom")

        # μ_N from N₂ molecule (cpd folder in 6.Formation_Energy_Diagram)
        n2_vasprun = (
            BASE.parent
            / "diamonty"
            / "6.Formation_Energy_Diagram"
            / "cpd"
            / "mol_N2"
            / "vasprun.xml"
        )
        v_n2 = Vasprun(str(n2_vasprun), parse_potcar_file=False)
        e_n2 = v_n2.final_energy
        n_n = sum(1 for site in v_n2.final_structure if site.species_string == "N")
        mu_N = e_n2 / n_n
        print(f"  μ_N (from N₂ molecule, {n_n} N atoms) = {mu_N:.6f} eV/atom")

        # These are the C-rich (diamond) and N₂-rich chemical potentials.
        # For a full stability diagram, you would also need the C-poor limit
        # (μ_C from graphite or another carbon phase).

        # ── 3. Defect energies ──
        v_perf = Vasprun(str(base / "perfect" / "vasprun.xml"), parse_potcar_file=False)
        e_perf = v_perf.final_energy
        n_atoms = v_perf.final_structure.num_sites
        print(f"  E_perf = {e_perf:.6f} eV  ({n_atoms} atoms)")

        charge_map = {
            -3: "N_C-V_C_-3",
            -2: "N_C-V_C_-2",
            -1: "N_C-V_C_-1",
            0: "N_C-V_C_0",
            1: "N_C-V_C_1",
            2: "N_C-V_C_2",
        }

        # Energy corrections (meV -> eV)
        corrections = {
            -3: 2.906,
            -2: 1.362,
            -1: 0.388,
            0: 0.0,
            1: 0.173,
            2: 0.799,
        }

        # For NV center: remove 1 C, add 1 N
        # E_form = E_def - E_perf + μ_C - μ_N + q·(E_VBM + E_F) + E_corr
        # The +μ_C - μ_N accounts for: -(-1)*μ_C - (+1)*μ_N = +μ_C - μ_N
        # Wait, let me re-derive:
        # E_form = E_def - E_perf - Σ n_i μ_i + q·(E_VBM + E_F) + E_corr
        # n_C = -1 (removed), n_N = +1 (added)
        # - Σ n_i μ_i = -[(-1)*μ_C + (+1)*μ_N] = μ_C - μ_N

        e_form_const = 2 * mu_C - mu_N  # constant part from chemical potentials

        print(
            f"\n{'Charge':<8} {'E_def (eV)':<16} {'E_corr (eV)':<14} {'E_form(q,0)':<16}"
        )
        print("-" * 60)

        e_defs = {}
        e_form_q0 = {}
        for q in sorted(charge_map.keys()):
            folder = charge_map[q]
            v = Vasprun(str(base / folder / "vasprun.xml"), parse_potcar_file=False)
            e_def = v.final_energy
            e_defs[q] = e_def

            # Formation energy at E_F = 0 (i.e., at VBM)
            e_form_0 = e_def - e_perf + e_form_const + q * vbm + corrections[q]
            e_form_q0[q] = e_form_0

            print(
                f"  {q:+d}      {e_def:<16.6f} {corrections[q]:<14.3f} {e_form_0:<16.4f}"
            )

        # ── 4. Plot ──
        fig, ax = plt.subplots(figsize=(8, 5.5))

        e_f_grid = np.linspace(0, band_gap, 200)
        colors = ["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b"]

        for i, q in enumerate(sorted(charge_map.keys())):
            e_form = e_form_q0[q] + q * e_f_grid
            label = f"q = {q:+d}"
            lw = 2.5 if q == -1 else 1.5
            ls = "--" if q == 0 else "-"
            ax.plot(
                e_f_grid,
                e_form,
                color=colors[i],
                linewidth=lw,
                linestyle=ls,
                label=label,
            )

        # Find charge transition levels (where lines cross)
        print(f"\n{'Transition':<15} {'E_F (eV)':<12} {'E_form (eV)'}")
        print("-" * 45)
        for i, q1 in enumerate(sorted(charge_map.keys())):
            for q2 in sorted(charge_map.keys()):
                if q2 <= q1:
                    continue
                # E_form(q1) = E_form(q2) at transition
                # e_form_q0[q1] + q1*E_F = e_form_q0[q2] + q2*E_F
                # E_F = (e_form_q0[q2] - e_form_q0[q1]) / (q1 - q2)
                e_f_trans = (e_form_q0[q2] - e_form_q0[q1]) / (q1 - q2)
                e_form_trans = e_form_q0[q1] + q1 * e_f_trans
                if 0 <= e_f_trans <= band_gap:
                    print(
                        f"  {q1:+d}/{q2:+d}         {e_f_trans:<12.3f} {e_form_trans:<.3f}"
                    )
                    ax.axvline(x=e_f_trans, color="gray", linewidth=0.5, linestyle=":")
                    ax.plot(e_f_trans, e_form_trans, "ko", markersize=4)

        ax.set_xlabel("Fermi level E$_F$ (eV)", fontsize=11)
        ax.set_ylabel("Formation energy (eV)", fontsize=11)
        ax.set_title("NV Center — Formation Energy Diagram", fontsize=12)
        ax.legend(fontsize=9, loc="upper left", ncol=2)
        ax.set_xlim(0, band_gap)
        ax.set_ylim(bottom=min(e_form_q0.values()) - 1)

        fig.tight_layout()
        fig.savefig(PLOTS_SLIDE / "formation_energy_diagram.pdf")
        print(f"\n  Saved {PLOTS_SLIDE / 'formation_energy_diagram.pdf'}")
        return fig

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

    def spin_states(self):
        """
        Extract total magnetization (magnetic moment) for each charge state
        from OUTCAR files. This gives the spin state of the NV center.

        Spin state = magnetization / 2  (since each unpaired electron contributes 1 μ_B)
        """
        base = POINT_DEFECT

        charge_order = [
            ("perfect", "perfect"),
            ("N_C-V_C_2", "+2"),
            ("N_C-V_C_1", "+1"),
            ("N_C-V_C_0", "0"),
            ("N_C-V_C_-1", "-1"),
            ("N_C-V_C_-2", "-2"),
            ("N_C-V_C_-3", "-3"),
        ]

        print(
            f"\n{'System':<12} {'Charge':<8} {'Magnet. (μB)':<14} {'Spin (S)':<10} {'Unpaired e⁻'}"
        )
        print("-" * 60)

        results = {}
        for folder, q_label in charge_order:
            outcar = base / folder / "OUTCAR"
            mag = None
            with open(outcar) as f:
                for line in f:
                    if "number of electron" in line and "magnetization" in line:
                        mag = float(line.split()[-1])
            spin_s = round(mag / 2, 1) if mag is not None else "?"
            unpaired = round(mag) if mag is not None else "?"
            results[q_label] = mag
            print(f"  {folder:<12} {q_label:<8} {mag:<14} {spin_s:<10} {unpaired}")

        return results

    def plot_spin_states(self, spin_data):
        """Plot magnetization vs charge state."""
        charges = [k for k in spin_data if k != "perfect"]
        mags = [spin_data[c] for c in charges]

        fig, ax = plt.subplots(figsize=(6, 4))
        colors = ["steelblue" if m > 0 else "crimson" for m in mags]
        bars = ax.bar(charges, mags, color=colors, edgecolor="black", linewidth=0.8)

        for bar, mag in zip(bars, mags):
            ax.text(
                bar.get_x() + bar.get_width() / 2,
                bar.get_height() + 0.05,
                f"S = {mag / 2:.1f}",
                ha="center",
                va="bottom",
                fontsize=10,
            )

        ax.set_xlabel("Charge state", fontsize=11)
        ax.set_ylabel("Total magnetization (μ_B)", fontsize=11)
        ax.set_title("NV Center — Spin States", fontsize=12)
        ax.set_ylim(0, max(mags) + 0.5)

        fig.tight_layout()
        fig.savefig(PLOTS_SLIDE / "spin_states.pdf")
        print(f"  Saved {PLOTS_SLIDE / 'spin_states.pdf'}")
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
                ax.plot(
                    0,
                    e,
                    "o",
                    color=c,
                    markerfacecolor=f,
                    markersize=10,
                    markeredgewidth=1.2,
                    zorder=2,
                )

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

    def zpl(self):
        """
        Compute the Zero-Phonon Line (ZPL) for NV⁻ → NV⁻* optical transition.

        ZPL = E_excited - E_ground

        Uses constrained-occupation DFT (ΔSCF method):
        - Ground state: NV⁻ from Point_defect
        - Excited state: NV⁻* from 7.Optical_Transitions/ZPL/

        The excited state promotes one electron from an occupied
        defect level to a higher empty defect level (spin-down channel).
        """
        # ── Ground state: NV⁻ (q = -1) ──
        v_ground = Vasprun(
            str(POINT_DEFECT / "N_C-V_C_-1" / "vasprun.xml"),
            parse_potcar_file=False,
        )
        e_ground = v_ground.final_energy
        mag_ground = v_ground.final_structure.site_properties.get("magmom", [0])

        # ── Excited state: NV⁻* ──
        v_excited = Vasprun(
            str(ZPL_EXCITED_DIR / "vasprun.xml"),
            parse_potcar_file=False,
        )
        e_excited = v_excited.final_energy

        zpl = e_excited - e_ground

        # Also extract the KS eigenvalues to see which electron was promoted
        perf = Vasprun(
            str(POINT_DEFECT / "perfect" / "vasprun.xml"),
            parse_potcar_file=False,
        )
        all_e = np.array(list(perf.eigenvalues.values()))
        vbm = float(all_e[all_e < perf.efermi].max())

        print(f"\n--- Zero-Phonon Line (ZPL) — NV⁻ Center ---")
        print(f"  E_ground  = {e_ground:.4f} eV")
        print(f"  E_excited = {e_excited:.4f} eV")
        print(f"  ZPL       = {zpl:.4f} eV  ({1239.84 / zpl:.1f} nm)")
        print(f"  Experiment = 1.945 eV (637 nm)")
        print(f"  Deviation  = {(zpl - 1.945) * 1000:.0f} meV")

        # Show gap states for both to identify the transition
        print(f"\n{'State':<20} {'Spin':<10} {'E (eV)':<12} {'E - VBM':<12} {'Occ'}")
        print("-" * 70)
        for label, v in [("Ground (NV⁻)", v_ground), ("Excited (NV⁻*)", v_excited)]:
            for spin, arr in v.eigenvalues.items():
                flat_e = arr[:, :, 0].flatten()
                flat_occ = arr[:, :, 1].flatten()
                mask = (flat_e > vbm - 0.5) & (flat_e < vbm + 4.0)
                for e, occ in zip(flat_e[mask], flat_occ[mask]):
                    if 0.01 < occ < 0.99 or (vbm < e < vbm + 3.5):
                        print(
                            f"  {label:<20} {spin.name:<10} {e:<12.4f} {e - vbm:<12.4f} {occ:.2f}"
                        )

        return e_ground, e_excited, zpl

    def plot_zpl(self, e_ground, e_excited, zpl):
        """Plot energy level diagram for the optical transition."""
        fig, ax = plt.subplots(figsize=(5, 4))

        y_ground = 0
        y_excited = zpl

        # Ground state level
        ax.barh(y_ground, 1, height=0.3, color="steelblue", label="NV⁻ (ground)")
        # Excited state level
        ax.barh(y_excited, 1, height=0.3, color="crimson", label="NV⁻* (excited)")

        # Arrow
        ax.annotate(
            "",
            xy=(0.5, y_excited - 0.15),
            xytext=(0.5, y_ground + 0.15),
            arrowprops=dict(arrowstyle="->", color="darkgreen", lw=2),
        )
        ax.text(
            0.55,
            (y_ground + y_excited) / 2,
            f"ZPL = {zpl:.2f} eV\n({1239.84 / zpl:.0f} nm)",
            fontsize=11,
            color="darkgreen",
            va="center",
        )

        ax.set_ylabel("Energy (eV)")
        ax.set_title("NV⁻ Optical Transition — Zero-Phonon Line")
        ax.set_yticks([y_ground, y_excited])
        ax.set_yticklabels([f"{e_ground:.2f}", f"{e_excited:.2f}"])
        ax.set_xticks([])
        ax.legend(fontsize=9, loc="upper right")
        ax.set_xlim(0, 1.0)

        fig.tight_layout()
        fig.savefig(PLOTS_SLIDE / "zpl.pdf")
        print(f"  Saved {PLOTS_SLIDE / 'zpl.pdf'}")
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
            e_form = e_def - e_perf - 2 * MU_C + MU_N
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
