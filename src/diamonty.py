from dataclasses import dataclass, field
from pathlib import Path
from typing import Literal, Optional

from constants import cutoff_dirs_pbe
from utils import read_oszicar_energy


@dataclass
class Config:
    functional: Literal["PBE", "HS56"] = field(
        default="PBE",
        metadata={"description": "Functional to use for the calculation"},
    )
    cutoff_dirs: list[Path] = field(
        default_factory=lambda: cutoff_dirs_pbe,
        metadata={"description": "Directories containing cutoff energy data"},
    )


class Diamonty:
    def __init__(self, config: Config):
        self.config = config

    def encut(self):
        """
        Returns the cutoff energy for the calculation.
        """
        cutoff_data = []
        for d in self.config.cutoff_dirs:
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


config = Config(functional="PBE")
diamonty = Diamonty(config)
diamonty.encut()
