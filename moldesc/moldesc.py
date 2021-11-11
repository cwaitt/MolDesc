""" Molecular Properties """

import numpy as np
from ase.data import vdw_radii, chemical_symbols
from ase.neighborlist import NeighborList, natural_cutoffs
from ase.io import read
from ase.visualize import view

vdw_cutoff_list = [[i,j] for i, j in zip(chemical_symbols, vdw_radii)] # list of VDW radii for each element in ASE
                                                                       # this could be replaced with your own radii


""" 
Molecular volume using VDW radii

input:

    mol - ase atoms object (if in a unit cell center object in unit cell)

return:
    Total_Volume - molecular VDW volume in A**3

"""
def mol_vdw_volume(mole):
    
    atoms = mole
    chem_sym = atoms.get_chemical_symbols()

    n_list = natural_cutoffs(atoms, mult=1) # use covalent radii to find connectivity of atoms

    nl = NeighborList(n_list,self_interaction=False,bothways=True) # get connectivity matrix
    nl.update(atoms)
    
    Total_volume = 0
    for i in range(0,len(atoms)):
        radii_a = [l[1] for l in vdw_cutoff_list if atoms.get_chemical_symbols()[i] in l][0] # find the radii of atom A
        ind, offsets = nl.get_neighbors(i) # get atom As neighbors

        radii_cross = 0 # cross section volume
        for j in ind:
            radii_b = [l[1] for l in vdw_cutoff_list if atoms.get_chemical_symbols()[j] in l][0] # for each neighbor, get the radius of the neigbor
            dist = atoms.get_distance(i,j,mic=True) # get the distance between A and B
            total_sphere_V = (radii_a + radii_b - dist)**2
            total_org_sphere = ((2*dist*radii_a) - (3*radii_a**2) + (6*radii_a*radii_b))
            total_cro_sphere = ((dist**2) + (2*dist*radii_b) - (3*radii_b**2))
            radii_cross = (((np.pi*total_sphere_V*(total_org_sphere+total_cro_sphere))/(12*dist))*0.5) + radii_cross # compute the volumetric intersection of the two spheres

        V_a = (4/3)*np.pi * radii_a**3 # get the volume of atom A
        V_cont = V_a - radii_cross # subtract off half of the intersecting volume (as B will subtract off the other half later)

        Total_volume = V_cont + Total_volume # Total VDW volume

    return Total_volume

"""
Molecular elongation

input:

    mol - ase atoms object (user must self align OSDA)

return:
    elong_a - first elongation coeff.
    elong_b - second elongation coeff.

"""

def mol_elongation(mol):
        
    atoms = mol
    positions = [i for i in atoms.positions]

    xpos = [i[0] for i in positions]
    ypos = [i[1] for i in positions]
    zpos = [i[2] for i in positions]

    xbox = max(xpos) - min(xpos)     # get the dimensions of the box to enclose the molecule
    ybox = max(ypos) - min(ypos)
    zbox = max(zpos) - min(zpos)

    atoms.set_cell([xbox,ybox,zbox])
    atoms.center()
    view(atoms)

    ### Elongation
    if xbox > ybox and xbox > zbox:
        elong_a = ybox/xbox
        elong_b = zbox/xbox

    elif ybox > xbox and ybox > zbox:
        elong_a = xbox/ybox
        elong_b = zbox/ybox

    elif zbox > xbox and zbox > ybox:
        elong_a = xbox/zbox
        elong_b = ybox/zbox
        
    return elong_a, elong_b
