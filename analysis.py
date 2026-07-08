#!/usr/bin/env python3
"""
Analysis of Diamond Bulk Properties from VASP calculations.

Tasks:
1. Cutoff energy convergence (ENCUT) - criterion: 4 meV
2. K-point mesh convergence - criterion: 1 meV
3. PBE: lattice parameter, DOS, LDOS, band gap, band structure
4. HSE06: lattice parameter, DOS, LDOS, band gap, band structure
"""

import json
import os
from pathlib import Path

import numpy as np

BASE = Path(__file__).parent / "vasp"
SLIDES_DATA = Path(__file__).parent.parent / "diamonty" / "slides" / "data"

# ── Experimental references ─────────────────────────────────
# Diamond, cubic (Fd-3m), 2 atoms/primitive cell.
# Lattice constant: 3.56683 Å at 300 K, rounded to 3.567 Å.
#   TODO: cite source
# Indirect band gap: 5.47 eV.
#   TODO: cite source
EXP_A = 3.567  # cubic lattice constant (Å)
EXP_GAP = 5.47  # indirect band gap (eV)


def read_oszicar_energy(path):
    """Extract E0 from the final line of OSZICAR."""
    with open(path) as f:
        lines = f.readlines()
    for line in reversed(lines):
        if "E0=" in line:
            return float(line.split("E0=")[1].split()[0])
    return None


def read_doscar(path):
    """Read DOSCAR file, return energy array, DOS array, and E_fermi."""
    with open(path) as f:
        lines = f.readlines()
    parts = lines[5].split()
    e_max = float(parts[0])
    e_min = float(parts[1])
    nedos = int(parts[2])
    e_fermi = float(parts[3])

    energy = []
    dos_total = []
    for line in lines[6 : 6 + nedos]:
        vals = [float(x) for x in line.split()]
        energy.append(vals[0])
        dos_total.append(vals[1] + vals[2])

    return np.array(energy), np.array(dos_total), e_fermi


def read_eigenval(path):
    """Read EIGENVAL file, return band energies array and occupancies."""
    with open(path) as f:
        lines = f.readlines()

    header_parts = [x for x in lines[5].split() if x]
    n_electrons = int(header_parts[0])
    n_kpts = int(header_parts[1])
    n_bands = int(header_parts[2])

    kpt_energies = []
    line_idx = 7
    for ik in range(n_kpts):
        line_idx += 1  # skip k-point header line
        band_energies = []
        for ib in range(n_bands):
            parts = lines[line_idx].split()
            band_energies.append(float(parts[1]))
            line_idx += 1
        kpt_energies.append(band_energies)
        line_idx += 1  # skip blank line

    kpt_energies = np.array(kpt_energies)
    n_occ = n_electrons // 2  # number of fully occupied bands

    vbm = np.max(kpt_energies[:, :n_occ])
    cbm = np.min(kpt_energies[:, n_occ:])
    band_gap = cbm - vbm

    return n_electrons, n_kpts, n_bands, n_occ, vbm, cbm, band_gap, kpt_energies


# ============================================================
# TASK 1: Cutoff Convergence (ENCUT)
# ============================================================
print("=" * 70)
print("TASK 1: CUTOFF ENERGY CONVERGENCE (ENCUT)")
print("Criterion: Energy difference < 4 meV")
print("=" * 70)

cutoff_dirs = sorted(BASE.glob("convergence/cutoff/*"), key=lambda p: int(p.name))

cutoff_data = []
for d in cutoff_dirs:
    encut = int(d.name)
    e0 = read_oszicar_energy(d / "OSZICAR")
    if e0 is not None:
        cutoff_data.append((encut, e0))

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

# ============================================================
# TASK 2: K-point Convergence
# ============================================================
print("\n" + "=" * 70)
print("TASK 2: K-POINT MESH CONVERGENCE")
print("Criterion: Energy difference < 1 meV")
print("=" * 70)

kdensity_dirs = sorted(BASE.glob("convergence/kdensity/*"), key=lambda p: int(p.name))

kpoint_data = []
for d in kdensity_dirs:
    idx = int(d.name)
    with open(d / "KPOINTS") as f:
        kp_lines = f.readlines()
    mesh = [int(x) for x in kp_lines[3].split()]
    nk = mesh[0]

    e0 = read_oszicar_energy(d / "OSZICAR")
    if e0 is not None:
        kpoint_data.append((idx, nk, e0))

print(f"\n{'Index':<8} {'k-mesh':<16} {'E0 (eV)':<20} {'ΔE (meV)':<15} {'Converged?'}")
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


# ============================================================
# TASK 3: PBE Analysis
# ============================================================
print("\n" + "=" * 70)
print("TASK 3: PBE FUNCTIONAL ANALYSIS")
print("=" * 70)

# --- Lattice parameter from CONTCAR ---
pbe_relax_dir = BASE / "PBE" / "relax"
with open(pbe_relax_dir / "CONTCAR") as f:
    contcar = f.readlines()

scale = float(contcar[1].strip())
a1 = scale * np.array([float(x) for x in contcar[2].split()])
a2 = scale * np.array([float(x) for x in contcar[3].split()])
a3 = scale * np.array([float(x) for x in contcar[4].split()])

lat_const_pbe = np.linalg.norm(a1) * np.sqrt(2)
vol_pbe = np.abs(np.dot(a1, np.cross(a2, a3)))

print(f"\n--- Lattice Parameters (PBE) ---")
print(f"  Primitive cell vectors (Angstrom):")
print(f"    a1 = {a1}")
print(f"    a2 = {a2}")
print(f"    a3 = {a3}")
print(f"  Primitive cell volume: {vol_pbe:.4f} Å³")
print(f"  Cubic lattice constant a = {lat_const_pbe:.4f} Å")
print(f"  Experimental a = {EXP_A} Å")
print(f"  Deviation: {100 * (lat_const_pbe - EXP_A) / EXP_A:+.2f}%")

# --- DOS Analysis ---
pbe_dos_dir = BASE / "PBE" / "dos"
print(f"\n--- Density of States (PBE) ---")
energy, dos_total, e_fermi = read_doscar(pbe_dos_dir / "DOSCAR")

print(f"  E_Fermi = {e_fermi:.4f} eV")
print(f"  Energy range: [{energy[0]:.2f}, {energy[-1]:.2f}] eV")
print(f"  NEDOS = {len(energy)}")

# Band gap from DOS
mask = (energy > e_fermi - 15) & (energy < e_fermi + 15)
e_win = energy[mask]
dos_win = dos_total[mask]

below = e_win < e_fermi
above = e_win > e_fermi

nb = np.where((dos_win[below] > 1e-6))[0]
vbm = (e_win[below][nb[-1]] - e_fermi) if len(nb) > 0 else 0.0

na = np.where((dos_win[above] > 1e-6))[0]
cbm = (e_win[above][na[0]] - e_fermi) if len(na) > 0 else 0.0

print(f"  Band gap from DOS: {cbm - vbm:.3f} eV")
print(f"  VBM: {vbm:.3f} eV, CBM: {cbm:.3f} eV (rel. to E_F)")

# --- Band Structure ---
pbe_band_dir = BASE / "PBE" / "band"
print(f"\n--- Band Structure (PBE) ---")
ne_pbe, nk_pbe, nb_pbe, no_pbe, vbm_pbe_b, cbm_pbe_b, bg_pbe, bands_pbe = read_eigenval(
    pbe_band_dir / "EIGENVAL"
)

print(f"  N_electrons = {ne_pbe}, N_kpts = {nk_pbe}, N_bands = {nb_pbe}")
print(f"  Occupied bands: {no_pbe}")
print(f"  VBM = {vbm_pbe_b:.4f} eV")
print(f"  CBM = {cbm_pbe_b:.4f} eV")
print(f"  Band gap = {bg_pbe:.4f} eV")
print(f"  VBM rel. E_F: {vbm_pbe_b - e_fermi:.4f} eV")
print(f"  CBM rel. E_F: {cbm_pbe_b - e_fermi:.4f} eV")

# --- LDOS ---
print(f"\n--- Local DOS (PBE) ---")
procar_pbe = pbe_dos_dir / "PROCAR"
if procar_pbe.exists():
    print(f"  PROCAR: {os.path.getsize(procar_pbe)} bytes")
    with open(procar_pbe) as f:
        first = f.readline().strip()
    print(f"  Header: {first}")
    # Extract atom-projected DOS info (first PROCAR block)
    with open(procar_pbe) as f:
        all_procar = f.readlines()
    # Find "tot" line for total projection
    tot_lines = [l for l in all_procar[:500] if "tot" in l.lower()]
    if tot_lines:
        print(f"  Total projections (first k-point): {tot_lines[0].strip()}")
else:
    print("  PROCAR not found")


# ============================================================
# TASK 4: HSE06 Analysis
# ============================================================
print("\n" + "=" * 70)
print("TASK 4: HSE06 FUNCTIONAL ANALYSIS")
print("=" * 70)

# --- Lattice parameter ---
hse06_relax_dir = BASE / "HSE06" / "relax"
with open(hse06_relax_dir / "CONTCAR") as f:
    contcar_hse = f.readlines()

scale_h = float(contcar_hse[1].strip())
a1_h = scale_h * np.array([float(x) for x in contcar_hse[2].split()])
a2_h = scale_h * np.array([float(x) for x in contcar_hse[3].split()])
a3_h = scale_h * np.array([float(x) for x in contcar_hse[4].split()])

lat_const_hse = np.linalg.norm(a1_h) * np.sqrt(2)
vol_hse = np.abs(np.dot(a1_h, np.cross(a2_h, a3_h)))

print(f"\n--- Lattice Parameters (HSE06) ---")
print(f"  Primitive cell vectors (Angstrom):")
print(f"    a1 = {a1_h}")
print(f"    a2 = {a2_h}")
print(f"    a3 = {a3_h}")
print(f"  Primitive cell volume: {vol_hse:.4f} Å³")
print(f"  Cubic lattice constant a = {lat_const_hse:.4f} Å")
print(f"  Experimental a = {EXP_A} Å")
print(f"  Deviation: {100 * (lat_const_hse - EXP_A) / EXP_A:+.2f}%")

# --- DOS ---
hse06_dos_dir = BASE / "HSE06" / "dos"
print(f"\n--- Density of States (HSE06) ---")
energy_h, dos_total_h, e_fermi_h = read_doscar(hse06_dos_dir / "DOSCAR")

print(f"  E_Fermi = {e_fermi_h:.4f} eV")
print(f"  Energy range: [{energy_h[0]:.2f}, {energy_h[-1]:.2f}] eV")

mask_h = (energy_h > e_fermi_h - 15) & (energy_h < e_fermi_h + 15)
e_win_h = energy_h[mask_h]
dos_win_h = dos_total_h[mask_h]

below_h = e_win_h < e_fermi_h
above_h = e_win_h > e_fermi_h

nb_h = np.where((dos_win_h[below_h] > 1e-6))[0]
vbm_h = (e_win_h[below_h][nb_h[-1]] - e_fermi_h) if len(nb_h) > 0 else 0.0

na_h = np.where((dos_win_h[above_h] > 1e-6))[0]
cbm_h = (e_win_h[above_h][na_h[0]] - e_fermi_h) if len(na_h) > 0 else 0.0

print(f"  Band gap from DOS: {cbm_h - vbm_h:.3f} eV")
print(f"  VBM: {vbm_h:.3f} eV, CBM: {cbm_h:.3f} eV (rel. to E_F)")

# --- Band Structure ---
hse_band_dir = BASE / "HSE06" / "band"
print(f"\n--- Band Structure (HSE06) ---")
ne_h, nk_h, nb_h, no_h, vbm_h_b, cbm_h_b, bg_h, bands_hse = read_eigenval(
    hse_band_dir / "EIGENVAL"
)

print(f"  N_electrons = {ne_h}, N_kpts = {nk_h}, N_bands = {nb_h}")
print(f"  Occupied bands: {no_h}")
print(f"  VBM = {vbm_h_b:.4f} eV")
print(f"  CBM = {cbm_h_b:.4f} eV")
print(f"  Band gap = {bg_h:.4f} eV")
print(f"  VBM rel. E_F: {vbm_h_b - e_fermi_h:.4f} eV")
print(f"  CBM rel. E_F: {cbm_h_b - e_fermi_h:.4f} eV")

# --- LDOS ---
print(f"\n--- Local DOS (HSE06) ---")
procar_hse = hse06_dos_dir / "PROCAR"
if procar_hse.exists():
    print(f"  PROCAR: {os.path.getsize(procar_hse)} bytes")
    with open(procar_hse) as f:
        first = f.readline().strip()
    print(f"  Header: {first}")
    with open(procar_hse) as f:
        all_p = f.readlines()
    tot_lines = [l for l in all_p[:500] if "tot" in l.lower()]
    if tot_lines:
        print(f"  Total projections (first k-point): {tot_lines[0].strip()}")
else:
    print("  PROCAR not found")


# ============================================================
# SUMMARY
# ============================================================
print("\n" + "=" * 70)
print("FINAL SUMMARY")
print("=" * 70)

print(f"""
+------------------------------+-------------------+-------------------+
| Property                     | PBE               | HSE06             |
+------------------------------+-------------------+-------------------+
| Lattice constant a (Å)       | {lat_const_pbe:.3f}            | {lat_const_hse:.3f}            |
| Primitive volume (Å³)        | {vol_pbe:.3f}            | {vol_hse:.3f}            |
| a deviation from exp.        | {100 * (lat_const_pbe - EXP_A) / EXP_A:+.2f}%            | {100 * (lat_const_hse - EXP_A) / EXP_A:+.2f}%            |
| Band gap - EIGENVAL (eV)     | {bg_pbe:.3f}            | {bg_h:.3f}            |
| VBM (eV)                     | {vbm_pbe_b:.3f}            | {vbm_h_b:.3f}            |
| CBM (eV)                     | {cbm_pbe_b:.3f}            | {cbm_h_b:.3f}            |
| E_Fermi - DOS (eV)           | {e_fermi:.3f}            | {e_fermi_h:.3f}            |
+------------------------------+-------------------+-------------------+

EXPERIMENTAL REFERENCES:
  Lattice constant: a = {EXP_A} Å  [TODO: cite source]
  Band gap (indirect): ~{EXP_GAP} eV  [TODO: cite source]

CONVERGENCE PARAMETERS:
  Cutoff energy (ENCUT) >= {converged_at} eV (ΔE < 4 meV)
  Recommended ENCUT = {rec_encut} eV
  K-point mesh >= {converged_at_k} (ΔE < 1 meV)
""")

# Save data for plotting
np.savez(
    "analysis_data.npz",
    cutoff_encut=[d[0] for d in cutoff_data],
    cutoff_e0=[d[1] for d in cutoff_data],
    kpoint_mesh=[d[1] for d in kpoint_data],
    kpoint_e0=[d[2] for d in kpoint_data],
    pbe_energy=energy,
    pbe_dos=dos_total,
    pbe_fermi=e_fermi,
    hse_energy=energy_h,
    hse_dos=dos_total_h,
    hse_fermi=e_fermi_h,
    pbe_bands=bands_pbe,
    pbe_e_fermi=e_fermi,
    hse_bands=bands_hse,
    hse_e_fermi=e_fermi_h,
    pbe_lat_const=lat_const_pbe,
    hse_lat_const=lat_const_hse,
    pbe_band_gap=bg_pbe,
    hse_band_gap=bg_h,
)

print("Data saved to analysis_data.npz for plotting.")

# ── JSON export for Typst tables ────────────────────────────
json_data = {
    "exp": {
        "a": EXP_A,
        "gap": EXP_GAP,
    },
    "pbe": {
        "a": round(lat_const_pbe, 4),
        "gap": round(bg_pbe, 4),
        "vbm": round(vbm_pbe_b, 4),
        "cbm": round(cbm_pbe_b, 4),
        "fermi": round(e_fermi, 4),
    },
    "hse06": {
        "a": round(lat_const_hse, 4),
        "gap": round(bg_h, 4),
        "vbm": round(vbm_h_b, 4),
        "cbm": round(cbm_h_b, 4),
        "fermi": round(e_fermi_h, 4),
    },
    "convergence": {
        "encut_converged": converged_at,
        "encut_recommended": rec_encut,
        "kpoint_converged": converged_at_k,
    },
    "encut_table": [
        {"encut": d[0], "e0": round(d[1], 6), "de": None}
        if i == 0
        else {
            "encut": d[0],
            "e0": round(d[1], 6),
            "de": round(abs(d[1] - cutoff_data[i - 1][1]) * 1000, 2),
        }
        for i, d in enumerate(cutoff_data)
    ],
    "kpoint_table": [
        {"mesh": f"{d[1]}×{d[1]}×{d[1]}", "e0": round(d[2], 6), "de": None}
        if i == 0
        else {
            "mesh": f"{d[1]}×{d[1]}×{d[1]}",
            "e0": round(d[2], 6),
            "de": round(abs(d[2] - kpoint_data[i - 1][2]) * 1000, 2),
        }
        for i, d in enumerate(kpoint_data)
    ],
}

SLIDES_DATA.mkdir(parents=True, exist_ok=True)
with open(SLIDES_DATA / "analysis_data.json", "w") as f:
    json.dump(json_data, f, indent=2)

print(f"Data saved to {SLIDES_DATA / 'analysis_data.json'} for Typst tables.")
