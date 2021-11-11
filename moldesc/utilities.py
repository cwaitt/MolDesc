""" Utility tools """

import numpy as np
from ase.io import read
from ase import Atoms

"""
Strip cell from atoms object and set to zero

input:
    mol - ase atoms object with unit cell

ouput:
    new_mol - ase atoms object without unitcell
"""


def strip(mole):
    atoms = mole
    positions = [i for i in atoms.positions]
    symbols = atoms.get_chemical_symbols()
    new_atoms = Atoms(symbols,positions=positions)
    pos0 = new_atoms[0].position
    trans = 0 - pos0
    new_atoms.translate(trans)

    new_mol = new_atoms
    
    return new_mol
