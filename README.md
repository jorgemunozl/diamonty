# Diamonty

**A DFT study of diamond for NV-center quantum technologies.**

## Motivation

The nitrogen-vacancy (NV) center in diamond is one of the most promising
platforms for room-temperature quantum computing, sensing, and
communication. Its spin state can be initialized, manipulated, and read
out optically — but designing devices around it requires a deep
understanding of the host material itself.

Before studying the defect, we must understand the *perfect crystal*.

## What this project does

This is a systematic density functional theory (DFT) study of bulk diamond
using the **VASP** plane-wave code. We characterize the ground-state
properties that set the stage for defect calculations:

1. **Convergence of computational parameters** — What plane-wave cutoff
   and k-point density are needed for converged total energies?
2. **Lattice structure** — How well do PBE (GGA) and HSE06 (hybrid)
   functionals reproduce the experimental lattice constant?
3. **Electronic structure** — Density of states (DOS), orbital-projected
   local DOS, and the fundamental band gap. How does the choice of
   exchange-correlation functional affect the predicted gap?

The project bridges *methodology* (how to converge a VASP calculation for
diamond) with *physics* (what the results tell us about bonding and
electronic structure).

## Methodology

| Detail | Value |
|---|---|
| Code | VASP 6.x (PAW method) |
| Pseudopotential | PAW_PBE (C) — same for both functionals |
| Functionals | PBE (GGA) and HSE06 (screened hybrid, 25% exact exchange) |
| ENCUT | 500 eV (production), converged to < 1 meV/atom |
| k-mesh | 10×10×10 Γ-centered (production) |
| Crystal | Diamond, Fd3̄m, 2 atoms per primitive cell |

All raw VASP output files (OUTCAR, DOSCAR, PROCAR, EIGENVAL, CONTCAR)
are included under `data/`.

## Key Findings

### Lattice constant

| Functional | a (Å) | Deviation from exp. |
|---|---|---|
| PBE | 3.574 | +0.19% |
| HSE06 | 3.548 | −0.53% |
| Experiment | 3.567 | — |

PBE slightly overestimates the lattice constant (typical GGA behavior);
HSE06 slightly underestimates it. Both are within 1% of experiment, which
is excellent for DFT.

### Band gap

| Functional | E_g (eV) | Deviation from exp. |
|---|---|---|
| PBE | 4.16 | −24.0% |
| HSE06 | 5.37 | −1.8% |
| Experiment | 5.47 | — |

The well-known *GGA band gap problem* is on full display: PBE underestimates
the indirect gap by nearly a quarter. HSE06 recovers it almost exactly —
the 25% exact Hartree-Fock exchange largely cancels the self-interaction
error that plagues semi-local functionals.

### Density of states

Both functionals produce qualitatively identical DOS: a clear gap at the
Fermi level, with the valence band centered around −5 eV and the conduction
band above +5 eV. The only difference is the *width* of the gap.

### Local DOS (orbital projection)

Integrated orbital weights from PROCAR (lm-decomposed, PAW sphere):

| Orbital | PBE | HSE06 |
|---|---|---|
| s | 0.83 | 0.84 |
| p | 1.74 | 1.76 |
| s/p ratio | 0.48 | 0.48 |

Both carbon atoms in the primitive cell are symmetry-equivalent and give
identical projections. The s/p ratio reflects mixed *sp*-bonding character
in the occupied bands. Crucially, the orbital decomposition is **invariant
under choice of functional** — the electronic structure *topology* is
robust, even when the gap is not.

## Interpretation

The results illustrate a broader principle in computational materials
science: **structural properties are cheap and reliable with GGA, but
electronic properties — especially band gaps — require hybrids.** For
the NV center project, this means:

- **PBE** is sufficient for geometry relaxation and convergence studies.
- **HSE06** (or beyond) is necessary for accurate defect levels, optical
  transitions, and spin properties.

A practical workflow: converge and relax with PBE, then compute electronic
properties with HSE06 on the PBE-relaxed geometry.

## Presentation

The `slides/` directory contains a Touying/Typst presentation
(`main.typ` → `main.pdf`) covering:

- Theoretical background (Kohn-Sham equations, pseudopotentials, SCF cycle)
- Convergence methodology (why ENCUT and k-points matter)
- All results above, with figures and tables
- Comparison of PBE vs HSE06 across every property

## Experimental references

- Lattice constant: 3.56683 Å at 300 K → rounded to 3.567 Å
  (TODO: add citation)
- Indirect band gap: 5.47 eV at 300 K
  (TODO: add citation — C.D. Clark et al., Proc. R. Soc. Lond. A 277, 312, 1964)

## Next steps

1. Band structure along high-symmetry k-path (Γ→X→W→K→Γ→L→U→W→L→K)
2. Single-atom substitutional defects (NV⁻ center)
3. Formation energies and charge-state transition levels
4. Excited-state forces for photoluminescence modeling
