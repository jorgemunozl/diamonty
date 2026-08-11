from pathlib import Path

BASE = Path(__file__).parent.parent

PLOTS_SLIDE = BASE / "slides" / "img"

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

# ── Band structure ─────────────────────────────────────────
_VASP = BASE.parent / "diamonty" / "4.Properties_of_Diamond" / "vasp"
PBE_BAND = _VASP / "PBE" / "band" / "vasprun.xml"
HSE06_BAND = _VASP / "HSE06" / "band" / "vasprun.xml"

# ── Supercell convergence (NV center formation energy) ─────
SUPERCELL_CONV = (
    BASE.parent / "diamonty" / "5.Diagramas_Kohn-Sham" / "Supercell_convergence"
)
MU_C = -9.092944955  # elemental chemical potential, C
MU_N = -8.320881190  # elemental chemical potential, N

POINT_DEFECT = BASE.parent / "diamonty" / "5.Diagramas_Kohn-Sham" / "Point_defect"
ZPL_EXCITED_DIR = BASE.parent / "diamonty" / "7.Optical_Transitions" / "ZPL" / "N_C-V_C_-1-excited"

# ── Experimental references ─────────────────────────────────
# Diamond, cubic (Fd-3m), 2 atoms/primitive cell.
# Lattice constant: 3.56683 Å at 300 K, rounded to 3.567 Å.
#   TODO: cite source
# Indirect band gap: 5.47 eV.
#   TODO: cite source
EXP_A = 3.567  # cubic lattice constant (Å)
EXP_GAP = 5.47  # indirect band gap (eV)
