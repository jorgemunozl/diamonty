#import "@preview/cetz:0.4.2"
#import "@preview/fletcher:0.5.8" as fletcher: edge, node
#import "@preview/touying:0.6.1": *
#import "@local/touying-simpl-uni:1.0.0": *

// cetz and fletcher bindings for touying
#let cetz-canvas = touying-reducer.with(reduce: cetz.canvas, cover: cetz.draw.hide.with(bounds: true))
#let fletcher-diagram = touying-reducer.with(reduce: fletcher.diagram, cover: fletcher.hide)

// ── Hardcoded data (run `python src/diamonty.py` to update) ─
#let exp-a = 3.567
#let exp-g = 5.47
#let dev(val, exp) = { 100 * (val - exp) / exp }

#let pbe-a = 3.5737
#let pbe-vol = 11.4103
#let hse-a = 3.5481
#let hse-vol = 11.1672

// LDOS data from PROCAR (integrated orbital weights per C atom)
#let pbe-s = 0.83
#let pbe-p = 1.74
#let hse-s = 0.84
#let hse-p = 1.76

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

= Part One: Theoretical Background Needed

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


= Part Two: Primitive Cell Results

== VASP workflow

#figure(
  image("img/vasp_workflow.png", width: 40%),
  caption: "VASP workflow",
)

== Questions to answer

1. What plane-wave energy cutoff should we use? (Cutoff Energy Convergence)

2. How fine should our k-point grid be? (K-point Mesh Convergence)

3. What does diamond look like?  (Lattice parameters, Density of States, Band gap and band structure)

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

#figure(
  convergence-table(
    ([*ENCUT (eV)*], [*$E_0$ (eV)*], [*$|\\Delta E|$ (meV)*]),
    (
      ([400], [-18.194445], [27.19]),
      ([450], [-18.187092], [7.35]),
      ([500], [-18.183833], [3.26]),
      ([550], [-18.184731], [0.90]),
      ([600], [-18.188449], [3.72]),
      ([650], [-18.192154], [3.71]),
      ([700], [-18.194773], [2.62]),
      ([750], [-18.196336], [1.56]),
      ([800], [-18.197682], [1.35]),
      ([850], [-18.198601], [0.92]),
      ([900], [-18.199201], [0.60]),
      ([950], [-18.199522], [0.32]),
    ),
    highlight: (3,),
  ),
  caption: [Energy convergence — criterion: 4 meV. Converged at 500 eV.],
)

== K-point Mesh Convergence

#figure(
  convergence-table(
    ([*k-mesh*], [*$E_0$ (eV)*], [*$|\\Delta E|$ (meV)*]),
    (
      ([5×5×5], [-18.149609], [548.83]),
      ([8×8×8], [-18.183083], [33.47]),
      ([10×10×10], [-18.183833], [0.75]),
      ([13×13×13], [-18.183888], [0.05]),
      ([15×15×15], [-18.183894], [0.01]),
      ([18×18×18], [-18.183890], [0.00]),
      ([20×20×20], [-18.183889], [0.00]),
      ([23×23×23], [-18.183894], [0.00]),
    ),
    highlight: (2,),
  ),
  caption: [k-point mesh convergence — criterion: 1 meV. Converged at 10×10×10.],
)



== Lattice Structure for PBE and HSE06

#figure(
  comparison-table(
    ([*Property*], [*PBE*], [*HSE06*], [*Experiment*]),
    (
      ([Lattice constant a (Å)], [#pbe-a], [#hse-a], [#exp-a]),
      (
        [Deviation from $a_"exp"$],
        [#calc.round(dev(pbe-a, exp-a), digits: 2)\%],
        [#calc.round(dev(hse-a, exp-a), digits: 2)\%],
        [—],
      ),
      ([Primitive volume V (Å³)], [#pbe-vol], [#hse-vol], [—]),
    ),
  ),
  caption: [Lattice parameters from VASP relaxation (ENCUT = 500 eV, 10³ k-mesh)],
)

#v(0.5em)

*PBE* overestimates the lattice constant by +0.19%, while *HSE06* underestimates it by -0.53%. Both are within 1% of experiment

== Density of States

#grid(
  columns: (1fr, 1fr),
  [
    #figure(
      image("img/dos_pbe.pdf", width: 100%),
      caption: [DOS (PBE)],
    )
  ],
  [
    #figure(
      image("img/dos_hse06.pdf", width: 100%),
      caption: [DOS (HSE06)],
    )
  ],
)


== Local Density of States (Orbital Projection)

#grid(
  columns: (1fr, 1fr),
  [
    #figure(
      image("img/ldos_pbe.pdf", width: 100%),
      caption: [LDOS — PBE],
    )
  ],
  [
    #figure(
      image("img/ldos_hse06.pdf", width: 100%),
      caption: [LDOS — HSE06],
    )
  ],
)


== K paths - High Symmetry

#grid(
  columns: (1fr, 1fr),
  [
    #figure(
      image("img/kpath.jpg", width: 130%),
      caption: [K path (PBE)],
    )
  ],
  [
    #figure(
      image("img/brillouin_fcc.jpg", width: 70%),
      caption: [Defect levels (HSE06)],
    )
  ],
)


== Band Structure

#grid(
  columns: (1fr, 1fr),
  [
    #figure(
      image("img/band_PBE.pdf", width: 100%),
      caption: [PBE ($E_g = 4.12$ eV)],
    )
  ],
  [
    #figure(
      image("img/band_HSE06.pdf", width: 100%),
      caption: [HSE06 ($E_g = 5.34$ eV)],
    )
  ],
)

= Part Three: Point defects on the super cell

== Why we cant use primitive cells for point defects

#grid(
  columns: (1fr, 1fr),
  [
    #figure(
      image("img/diamond_1.png", width: 80%),
      caption: [Diamond conventional cell, *8 atoms*],
    )
  ],
  [
    #figure(
      image("img/diamond_2.png", width: 80%),
      caption: [Diamond primitive cell, *2 atoms*],
    )
  ],
)

== How obtain the energy for point defects

*How much energy does it cost to take a perfect diamond supercell and turn into a supercell containing a point *

$ E^(q=0)_("form")[D] = E^(q=0)_("def")[D] - E_("perf") - 2 mu_(C) + mu_(N) $

1. Total energy difference between the defective supercell and the pristine supercell.
2. Energy to remove a carbon atom from the supercell, $mu_(C)=-9.092944955$, _depends on choice of reservoir and growth conditions_.
3. Energy to add a nitrogen atom to the supercell, $mu_(N)=-8.32088119$

Note: You actually should remove two carbon atoms.

== NV Center — Supercell Convergence

#figure(
  image("img/formation_energy.pdf", width: 64%),
  caption: [Formation energy vs supercell size — converges to 17.79 eV],
)

#v(0.5em)

#figure(
  convergence-table(
    ([*Supercell*], [*Atoms*], [*$E_"form"$ (eV)*], [*$|\\Delta E|$ (meV)*]),
    (
      ([1×1×1], [8], [16.634], [—]),
      ([2×2×2], [64], [17.586], [952.0]),
      ([3×3×3], [216], [17.753], [167.5]),
      ([4×4×4], [512], [17.780], [26.9]),
      ([5×5×5], [1000], [17.788], [7.6]),
      ([6×6×6], [1728], [17.789], [1.8]),
      ([7×7×7], [2744], [17.790], [0.7]),
    ),
    highlight: (3,),
  ),
  caption: [Convergence at 4×4×4 (512 atoms)],
)

By 4×4×4 (512 atoms) the formation energy is within 10 meV of the isolated-defect
limit. This is the practical converged supercell size for subsequent defect
calculations.

// == Defect Levels — NV Center #figure( image("img/defect_levels.png", width: 30%), caption: [Single-particle Kohn-Sham defect levels vs charge state], ) Each column shows gap states for one charge state ($-3$ to $+2$). Filled circles = occupied, open triangles = empty, blue = spin-up, red = spin-down. VBM and CBM from the perfect 4×4×4 supercell.


== Defect Creation

#figure(
  image("img/diamond.webp", width: 110%),
  caption: [asd],
)

== Defect Number of electrons

#figure(
  image("img/five_diamond.png", width: 40%),
  caption: [In a neutral NV center, there are 5 free electrons],
)

== What we want?

I don't know

== Why KS eigenvalues are the answer

Because yes asd

== KS Eigenvalues — All Charge States

#grid(
  columns: (1fr, 1fr, 1fr),
  rows: (auto, auto, auto),
  [
    #figure(image("img/eigen_N_C-V_C_-3.pdf", width: 100%), caption: [q = −3])
  ],
  [
    #figure(image("img/eigen_N_C-V_C_-2.pdf", width: 100%), caption: [q = −2])
  ],
  [
    #figure(image("img/eigen_N_C-V_C_-1.pdf", width: 100%), caption: [q = −1])
  ],

  [
    #figure(image("img/eigen_N_C-V_C_0.pdf", width: 100%), caption: [q = 0])
  ],
  [
    #figure(image("img/eigen_N_C-V_C_1.pdf", width: 100%), caption: [q = +1])
  ],
  [
    #figure(image("img/eigen_N_C-V_C_2.pdf", width: 100%), caption: [q = +2])
  ],
)

== Why the spin is important, because yes

So yeah

== Spin state for the NV center

#figure(
  image("img/spin_states.pdf", width: 60%),
  caption: [Spin state of the NV center vs charge state],
)


== How on earth,

== Master Equation

#grid(
  columns: (1fr, 1fr),
  [
    #text(size: 0.8em)[
      For a defect with charge $q$, the _master equation_:

      $
        E^f(q, E_F) = E_(text("def"))(q) - E_(text("perf")) - sum n_i mu_i + q(E_(text("VBM")) + E_F) + E_(text("corr"))
      $

      - $n_i$ is negative when atoms are added and positive when atoms are removed
      - $mu_i$ is the chemical potential of the $i$-th band
      - $E_F$ is the Fermi level
      - $E_(text("VBM"))$ is the valence band maximum energy
      - $E_(text("corr"))$ is the correction energy.

      Plot it when:

      $ 0 <= E_F <= E_(text("g")) $

      because
    ]
  ],
  [
    #figure(
      image("img/folder-tree.png", width: 70%),
      caption: [q = 0],
    )
  ],
)


== Master Equation for NV center on diamond

$ E^f(q, E_F) = E_(text("def"))(q) - E_(text("perf")) - sum n_i mu_i + q(E_(text("VBM")) + E_F) + E_(text("corr")) $


== Formation Energy Diagram

#figure(
  image("img/formation_energy_diagram.pdf", width: 60%),
  caption: [Formation energy diagram],
)

== Optical Transitions

When a defect is optical active?
What being optical active means?
