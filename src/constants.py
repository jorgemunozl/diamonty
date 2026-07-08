from pathlib import Path

BASE = Path(__file__).parent.parent

DATA = BASE / "data"
CUTOFF_DIRS = DATA / "convergence" / "cutoff"
KDENSITY_DIRS = DATA / "convergence" / "kdensity"

PBE_RELAX = DATA / "relax" / "PBE" / "RELAX_CONTCAR"
HSE06_RELAX = DATA / "relax" / "HSE06" / "RELAX_CONTCAR"

cutoff_dirs_pbe = [path for path in CUTOFF_DIRS.glob("*") if path.is_dir()]
kdensity_dirs_pbe = [path for path in KDENSITY_DIRS.glob("*") if path.is_dir()]

# ── DOS directories ────────────────────────────────────────
PBE_DOS_DIR = DATA / "dos" / "PBE"
HSE06_DOS_DIR = DATA / "dos" / "HSE06"

# ── Experimental references ─────────────────────────────────
# Diamond, cubic (Fd-3m), 2 atoms/primitive cell.
# Lattice constant: 3.56683 Å at 300 K, rounded to 3.567 Å.
#   TODO: cite source
# Indirect band gap: 5.47 eV.
#   TODO: cite source
EXP_A = 3.567  # cubic lattice constant (Å)
EXP_GAP = 5.47  # indirect band gap (eV)
