#import "@preview/cetz:0.4.2"
#import "@preview/fletcher:0.5.8" as fletcher: edge, node
#import "@preview/touying:0.6.1": *
#import "@local/touying-simpl-uni:1.0.0": *

// cetz and fletcher bindings for touying
#let cetz-canvas = touying-reducer.with(reduce: cetz.canvas, cover: cetz.draw.hide.with(bounds: true))
#let fletcher-diagram = touying-reducer.with(reduce: fletcher.diagram, cover: fletcher.hide)

// ── Load analysis data ──────────────────────────────────────
#let data = json("data/analysis_data.json")
#let exp-a = data.exp.a
#let exp-g = data.exp.gap
#let dev(val, exp) = { 100 * (val - exp) / exp }

#let de-str(de) = {
  if de != none {
    str(calc.round(de, digits: 2))
  } else {
    "—"
  }
}

// Table Number


#show: ecnu-theme.with(
  // Lang and font configuration
  lang: "en",
  font: "Libertinus Serif",

  // Basic information
  config-info(
    title: [Diamonty: NV Center in Diamond],
    short-title: [Diamonty: NV Center in Diamond],
    subtitle: [A study on punctual defects on materials with a high band-gap: Research Sprint],
    author: [Munoz Laredo Jorge],
    date: datetime.today(),
    institution: [National University of Engineering],
  ),
)
#title-slide()

//#outline-slide()

= Part One: Theoretical Background

== Many Electrons problem set up

We are trying to find the wave function $Psi(r_1, r_2, ..., r_(N e))$ for a many body system in his *ground state*, with it we could obtain all the information of a system, and create any kind of material with desired property.

Analytic is extremely hard, then we use the *BO aproximation* and a Hartree product, and limit ourselves to find the *electronic density* $n_0(r)=sum_j |psi_j^(K S)(r_j)|²$. _(Hokenberg Khon Theorem)_.

You end with a functional of the energy as:

$ E[n]=T_H [n] + U_H [n] + E_(X C)[n]+U_(e E)[n] $

Where:  $E_(X C) = Delta T + Delta U$ where

$T_H [n] = sum_i ⟨ psi_j^(K S) | -frac(nabla^2, 2) | psi_j^(K S) ⟩ quad and quad U_H [n] = 1/2 integral integral frac(n(bold(r)) n(bold(r)'), |bold(r) - bold(r)'|) dif bold(r) dif bold(r)'$

== Influent terms for the Kohn Sham equations

The choosing of, effective potential

$ V_(text("eff")) = V_H (r) + V_(e E) (r) + V_(x c)(r) $

By performing optimization using the *Variational Principle* you end up with the Kohn Sham equations.

$ H_"KS"[n] psi_j^(K S)(r) = epsilon_i psi_j^(K S)(r) $

Where the Khons Sham Hamiltonian can be written as

$
  H_"KS"[n] = -frac(1, 2) nabla^2 + V_"ext"(bold(r)) + V_H[n](bold(r)) + V_(x c)[n](bold(r))
$

where each term is: $V_H[n] = integral frac(n(bold(r)'), |bold(r) - bold(r)'|) dif bold(r)'$ is the Hartree potential,
$V_(x c)[n] = frac(delta E_(x c)[n], delta n(bold(r)))$ the exchange-correlation potential,
and $V_"ext"$ the external potential from the nuclei.


== Kohn Sham is self-referential

The Kohn Sham equations are *self-referential*: the Hamiltonian $H_"KS"$ depends
on the density $n(bold(r))$, but the density is built from the KS orbitals
$n(bold(r)) = sum_i f_i |psi_i(bold(r))|^2$. Unlike the fixed Schrödinger
Hamiltonian, every $V_H[n]$ and $V_(x c)[n]$ term must be recomputed when $n$
changes.

#highlight["The Hamiltonian determines the density, and the density determines
  the Hamiltonian."] This circular dependence forces an *iterative* solution — the Self-Consistent Field (SCF) cycle:

But before of that, we need a form to write down $psi(r)$, so we need to choose a basis set.

== A Basis Set where the wave functions are expressed

To solve the KS equations numerically, we expand the orbitals $psi_i(bold(r))$ in a *basis set*. The choice of basis dramatically affects computational cost and accuracy.

*Plane Wave (PW) Basis*

In periodic systems (crystals), Bloch's theorem allows us to write:

$ psi_(n bold(k))(bold(r)) = e^(i bold(k) dot.op bold(r)) u_(n bold(k))(bold(r)) $

where $u_(n bold(k))(bold(r))$ has the periodicity of the lattice. We expand $u$ in plane waves:

$ u_(n bold(k))(bold(r)) = sum_(bold(G)) c_(n bold(k))(bold(G)) e^(i bold(G) dot.op bold(r)) $

where $bold(G)$ are reciprocal lattice vectors. In practice,
the sum is truncated at a kinetic energy cutoff, but in practice is better use PAW basis.

// == PAW Basis
//
// The *Projector Augmented Wave (PAW)* method (Blöchl, 1994) is an all-electron approach that combines efficiency of pseudopotentials with full nodal structure near nuclei.
//
// *Core Idea:* The true wavefunction $psi$ is related to a smooth pseudo-wavefunction $tilde(psi)$ via:
//
// $ psi = tilde(psi) + sum_i (phi_i - tilde(phi)_i) ⟨ tilde(p)_i | tilde(psi) ⟩ $
//
// where $tilde(psi)$ is smooth (expanded in plane waves with low $E_"cut"$), $phi_i$ are all-electron partial waves (localized, contain nodes), $tilde(phi)_i$ are pseudo partial waves, and $tilde(p)_i$ are projector functions.
//
// #highlight[Key benefit:] Accurate magnetic and core properties while keeping plane-wave efficiency. Used by VASP, ABINIT, and other modern DFT codes.

== Pseudopotentials

#grid(
  columns: (1.5fr, 1fr),
  column-gutter: 1em,
  [
    Core electrons are tightly bound and chemically inert — they don't participate in bonding. However, their rapid oscillations near the nucleus require huge plane-wave cutoffs.

    *Pseudopotential Approximation:* Replace the all-electron potential $V_"ion"$ and core electrons with an effective *pseudopotential* $V_"pseudo"$ that:

    1. Reproduces the correct valence properties (scattering, eigenvalues)
    2. Has smooth, nodeless *pseudo-wavefunctions* outside a core radius $r_c$
  ],
  [
    #figure(
      image("img/pseudo.jpg", width: 100%),
      caption: "Pseudopotential Approximation",
    )
  ],
)

== Different Flavors of Pseudopotentials

#table(
  columns: (auto, 1.8fr, 1.8fr, 1.2fr),
  inset: 8pt,
  align: (left, left, left, left),
  stroke: none,
  fill: (x, y) => if y == 0 { rgb("#3b82f6").lighten(80%) } else if calc.odd(y) { rgb("#f1f5f9") },

  table.header([*Type*], [*Key Feature*], [*Trade-off*], [*Examples*]),

  [*NC*], [Preserves wavefunction norm], [High accuracy, high $E_"cut"$], [Troullier-Martins],

  [*Ultrasoft*], [Relaxed norm constraint], [Low $E_"cut"$, less transferable], [Vanderbilt],

  [*PAW*], [All-electron reconstruction], [Best accuracy, complex], [VASP, QE],

  [*ECP*], [Analytical core replacement], [Good for heavy atoms + relativity], [Stuttgart-Dresden],
)

== PBE and HSE06 pseudopotentials

#grid(
  columns: (1fr, 1fr),
  column-gutter: 1.5em,
  [
    *PBE (Perdew-Burke-Ernzerhof)*

    - GGA-level exchange-correlation functional
    - Semi-local: depends on $n(bold(r))$ and $nabla n(bold(r))$
    - Cheap, robust, works well for structures
    - *Underestimates band gaps* (systematic error)
  ],
  [
    *HSE06 (Heyd-Scuseria-Ernzerhof)*

    - Hybrid functional: mixes 25% exact Hartree-Fock exchange with PBE
    - Screened Coulomb potential (range-separation)
    - Much more expensive (~10--100×)
    - *Much better band gaps* (corrects self-interaction error)
  ],
)

Both use the *same PAW pseudopotential* (`PAW_PBE C`). The functional is
chosen in INCAR, not in the pseudopotential file.

== The Self-Consistent Field (SCF) Method is a loop

When running a SCF calculation, e.g. quantum espresso, what the SCF cycle does is:

1. *Initial Guess*: Start with an approximate $n^((0))(bold(r))$, often a superposition of atomic charge densities.

2. *Construct $H_"KS"$*: Build the KS Hamiltonian:
  $ H_"KS" = -frac(1, 2) nabla^2 + V_"ext" + V_"Hartree"[n] + V_(x c)[n] $

3. *Diagonalize*: Solve the generalized eigenvalue problem:
  $ H c_i = epsilon_i S c_i quad "(in a plane-wave or PAW basis)" $

== The self-consistent field (SCF) cycle

4. *New Density*: Compute the output density:
  $ n_"out"(bold(r)) = sum_i f_i |Psi_i(bold(r))|^2 quad (f_i = "occupation") $

5. *Check Convergence*: If $||n_"out" - n_"in"|| < delta arrow.r$ done. Otherwise: mix $n_"in"$ and $n_"out"$, return to step 2.


= How you are supposed to interpretate the data from the calculations (example)

== VASP workflow

#figure(
  image("img/vasp_workflow.png", width: 40%),
  caption: "VASP workflow",
)


== Questions to answer

1. What plane-wave energy cutoff should we use? (Cutoff Energy Convergence)

2. How fine should our k-point grid be? (K-point Mesh Convergence)

3. What does diamond look like?  (Lattice parameters, Density of States, Band gap and band structure)

== The data that we actually have

#table(
  columns: (1.6fr, 2.5fr, 2fr, 2fr),
  inset: 6pt,
  align: (left, left, left, left),
  stroke: none,
  fill: (x, y) => if y == 0 { rgb("#3b82f6").lighten(80%) } else if calc.odd(y) { rgb("#f1f5f9") },

  table.header([*Study*], [*Folder / Type*], [*Parameters explored*], [*Key outputs*]),

  [*Convergence (ENCUT)*],
  [`convergence/cutoff/`],
  [16 values: 200 to 950 eV (step 50)],
  [OUTCAR, OSZICAR, E0 vs ENCUT],

  [*Convergence (k-points)*], [`convergence/kdensity/`], [9 meshes: 3³ to 23³], [OUTCAR, OSZICAR, KPOINTS, E0 vs Nₖ],

  [*PBE Relaxation*], [`PBE/relax/`], [ENCUT = 500 eV, 10³ mesh], [CONTCAR (lattice const.), OUTCAR],

  [*PBE DOS + LDOS*], [`PBE/dos/`], [ENCUT = 500 eV, 10³ mesh], [DOSCAR (total DOS), PROCAR (orbital proj.)],

  [*PBE Band Structure*], [`PBE/band/`], [ENCUT = 500 eV, 900 kpts, 8 bands], [EIGENVAL, PROCAR],

  [*HSE06 Relaxation*], [`HSE06/relax/`], [ENCUT = 500 eV, 10³ mesh], [CONTCAR (lattice const.), OUTCAR],

  [*HSE06 DOS + LDOS*], [`HSE06/dos/`], [ENCUT = 500 eV, 10³ mesh], [DOSCAR (total DOS), PROCAR (orbital proj.)],

  [*HSE06 Band Structure*], [`HSE06/band/`], [ENCUT = 500 eV, 47 kpts, 24 bands], [EIGENVAL, PROCAR],
)

== Convergence of the energy from SCF (example)

#figure(
  image("img/energy_convergence.png", width: 85%),
)

== Density of States obtention for different pseudopotentials (example)

#figure(
  image("img/DOS.png", width: 70%),
)

== Band properties obtention (example)

#figure(
  image("img/unit_cell_band.png", width: 80%),
)

= Plotting of the data for Diamond

== Energy Cutoff for the plane-wave basis and k mesh

#grid(
  columns: (1fr, 1fr),
  [
    #figure(
      image("img/encut_convergence.png", width: 100%),
      caption: [Energy convergence as a function of ENCUT],
    )
  ],
  [
    #figure(
      image("img/kpoint_convergence.png", width: 100%),
      caption: [Energy convergence as a function of k-point mesh],
    )
  ],
)

== Energy Cutoff

#let encut-rows = data.encut_table.map(row => ([#row.encut], [#row.e0], [#de-str(row.de)]))

#figure(
  convergence-table(
    ([*ENCUT (eV)*], [*$E_0$ (eV)*], [*$|\\Delta E|$ (meV)*]),
    encut-rows,
    highlight: (2,),
  ),
  caption: [Energy convergence as a function of ENCUT],
)

== K-point Mesh Convergence

#let kpoint-rows = data.kpoint_table.map(row => ([#row.mesh], [#row.e0], [#de-str(row.de)]))

#figure(
  convergence-table(
    ([*k-mesh*], [*$E_0$ (eV)*], [*$|\\Delta E|$ (meV)*]),
    kpoint-rows,
    highlight: (3,),
  ),
  caption: [k-point mesh convergence — criterion: 1 meV. Converged at 10×10×10.],
)


== Lattice Structure for PBE and HSE06

#prop-table("PBE Functionasdal", (
  ([Lattice constant \\(a\\)], [#data.pbe.a\ Å]),
  ([Deviation from $a_"exp"$], [#calc.round(dev(data.pbe.a, exp-a), digits: 2)\%]),
  ([Band gap $E_g$], [#data.pbe.gap\ eV]),
  ([VBM (absolute)], [#data.pbe.vbm\ eV]),
  ([CBM (absolute)], [#data.pbe.cbm\ eV]),
  ([$E_"Fermi"$ (from DOS)], [#data.pbe.fermi\ eV]),
))

== PBE vs HSE06

#comparison-table(
  ([*Property*], [*PBE*], [*HSE06*], [*Experiment*]),
  (
    ([Lattice constant \\(a\\) (Å)], [#data.pbe.a], [#data.hse06.a], [#exp-a]),
    (
      [Deviation from exp.],
      [#calc.round(dev(data.pbe.a, exp-a), digits: 2)\%],
      [#calc.round(dev(data.hse06.a, exp-a), digits: 2)\%],
      [—],
    ),
    ([Band gap $E_g$ (eV)], [#data.pbe.gap], [#data.hse06.gap], [#exp-g]),
    (
      [Gap deviation from exp.],
      [#calc.round(dev(data.pbe.gap, exp-g), digits: 1)\%],
      [#calc.round(dev(data.hse06.gap, exp-g), digits: 1)\%],
      [—],
    ),
    ([VBM (eV)], [#data.pbe.vbm], [#data.hse06.vbm], [—]),
    ([CBM (eV)], [#data.pbe.cbm], [#data.hse06.cbm], [—]),
    ([$E_"Fermi"$ (eV)], [#data.pbe.fermi], [#data.hse06.fermi], [—]),
  ),
)
