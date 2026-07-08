from pathlib import Path

BASE = Path(__file__).parent.parent

DATA = BASE / "data"
CUTOFF_DIRS = DATA / "convergence" / "cutoff"

cutoff_dirs_pbe = [path for path in CUTOFF_DIRS.glob("*") if path.is_dir()]

# ── Experimental references ─────────────────────────────────
# Diamond, cubic (Fd-3m), 2 atoms/primitive cell.
# Lattice constant: 3.56683 Å at 300 K, rounded to 3.567 Å.
#   TODO: cite source
# Indirect band gap: 5.47 eV.
#   TODO: cite source
EXP_A = 3.567  # cubic lattice constant (Å)
EXP_GAP = 5.47  # indirect band gap (eV)
