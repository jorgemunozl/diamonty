"""Extract orbital-projected LDOS from PROCAR for slides."""

from pathlib import Path

from constants import HSE06_DOS_DIR, PBE_DOS_DIR


def parse_procar(path: Path):
    """Parse PROCAR lm-decomposed. Return {ion: {orbital: weight}}"""
    with open(path) as f:
        lines = f.readlines()

    # Read header line 2
    header = lines[1].split()
    n_kpts = int(header[3])
    n_bands = int(header[7])
    n_ions = int(header[11])

    orbitals = ["s", "py", "pz", "px", "dxy", "dyz", "dz2", "dxz", "dx2-y2"]
    ion_orbitals = {i: {orb: 0.0 for orb in orbitals} for i in range(1, n_ions + 1)}

    line_idx = 2
    for _ in range(n_kpts):
        # Skip to next k-point line
        while line_idx < len(lines) and "k-point" not in lines[line_idx]:
            line_idx += 1
        if line_idx >= len(lines):
            break

        kpt_line = lines[line_idx]
        for sep in ("weight =", "weight="):
            if sep in kpt_line:
                kpt_weight = float(kpt_line.split(sep)[1].strip())
                break
        else:
            kpt_weight = 1.0 / n_kpts
        line_idx += 1

        for _ in range(n_bands):
            # Skip to band line
            while line_idx < len(lines) and "band" not in lines[line_idx]:
                line_idx += 1
            if line_idx >= len(lines):
                break

            band_line = lines[line_idx]
            occ = float(band_line.split("occ.")[1].strip())
            line_idx += 1

            # Skip blank line and "ion ..." header, land on first digit line
            while line_idx < len(lines):
                stripped = lines[line_idx].strip()
                if stripped and stripped[0].isdigit():
                    break
                line_idx += 1

            if line_idx >= len(lines):
                break

            # Read ion rows
            for ion in range(1, n_ions + 1):
                parts = lines[line_idx].split()
                for j, orb in enumerate(orbitals):
                    ion_orbitals[ion][orb] += float(parts[j + 1]) * kpt_weight * occ
                line_idx += 1

            # Skip 'tot' line
            line_idx += 1

    return ion_orbitals, orbitals


def summarize(ion_orbitals, orbitals):
    """Total s, p, d per ion."""
    result = {}
    for ion, data in ion_orbitals.items():
        s = data["s"]
        p = data["py"] + data["pz"] + data["px"]
        d = sum(data[orb] for orb in orbitals if orb.startswith("d"))
        total = s + p + d
        result[f"C{ion}"] = {
            "s": round(s, 2),
            "p": round(p, 2),
            "total": round(total, 2),
            "s/p": round(s / p, 3) if p > 0 else 0,
        }
    return result


for functional, dos_dir in [("PBE", PBE_DOS_DIR), ("HSE06", HSE06_DOS_DIR)]:
    proc_path = dos_dir / "PROCAR"
    ion_data, orbitals = parse_procar(proc_path)
    summary = summarize(ion_data, orbitals)

    print(f"\n=== {functional} ===")
    for label, vals in summary.items():
        print(
            f"  {label}: s={vals['s']}, p={vals['p']}, "
            f"total={vals['total']}, s/p={vals['s/p']}"
        )
