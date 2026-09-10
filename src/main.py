from diamonty import Config, Diamonty

config = Config(functional="HSE06")
# config = Config(functional="PBE")
diamonty = Diamonty(config)


def formation_energy():
    results = diamonty.formation_energy()
    diamonty.plot_formation_energy(results)


def formation_energy_diagram():
    diamonty.formation_energy_diagram()


def band_structure():
    bs, gap = diamonty.band_structure()
    diamonty.plot_band_structure(bs)


def spin_states():
    spin_data = diamonty.spin_states()
    diamonty.plot_spin_states(spin_data)


def defect_levels():
    vbm, cbm, levels = diamonty.defect_levels()
    diamonty.plot_defect_levels(vbm, cbm, levels)


def zpl():
    e_ground, e_excited, zpl = diamonty.zpl()
    diamonty.plot_zpl(e_ground, e_excited, zpl)


def dos():
    dos = diamonty.dos()
    diamonty.plot_dos(dos)


def ldos():
    ldos = diamonty.ldos()
    diamonty.plot_ldos(ldos)


if __name__ == "__main__":
    formation_energy_diagram()
