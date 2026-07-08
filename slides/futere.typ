== Local Density of States (Orbital Projection)

Integrated orbital weights per carbon atom, summed over all occupied
bands and k-points (from PROCAR lm-decomposed).

#figure(
  comparison-table(
    ([*Atom*], [*Orbital*], [*PBE*], [*HSE06*]),
    (
      ([C₁,₂], [$s$], [#pbe-s], [#hse-s]),
      ([C₁,₂], [$p$], [#pbe-p], [#hse-p]),
      ([C₁,₂], [$s/p$ ratio], [#calc.round(pbe-s / pbe-p, digits: 2)], [#calc.round(hse-s / hse-p, digits: 2)]),
    ),
  ),
  caption: [Orbital-projected weights (PAW sphere) — both Carbons are symmetry-equivalent],
)

Both carbon atoms give identical projections (diamond $F d 3-bar m$ symmetry).
The $s/p$ ratio in the occupied bands corresponds to mixed $s$–$p$
bonding character. PBE and HSE06 yield nearly identical orbital
decompositions — the electronic structure *topology* is insensitive to
the choice of functional.

Both functionals show a clear band gap near $E_F$, with the valence band
around −5 eV and the conduction band above +5 eV — typical of an
insulator. The HSE06 gap is wider due to exact exchange.
