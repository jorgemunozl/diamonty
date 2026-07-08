import numpy as np


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
