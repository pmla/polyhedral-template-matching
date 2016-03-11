"""Python module for running PTM on an Asap list of atoms.

This module runs Polyhedral Template Matching on an Atoms object from ASE/Asap.
It requires that both ASE and Asap are installed.
"""

import asap3
import ptmmodule
import numpy as np


def PTM(atoms, structures=None, alloys=False, cutoff=10.0):
    """Run Complex Hull Analysis on an Atoms object.

    Parameters:
    atoms: The atoms object

    structures=None: The list of structures to be investigated.
        The list defaults to ('sc', 'fcc', 'hcp', 'ico', 'bcc').

    alloys=False: Set to true to check for alloy structures.
    
    cutoff: A cutoff used for the neighborlist.  Must be large enough that all
        nearest neighbors are returned (second-nearest for BCC).  Using a too
        large value may impact performance, but not results.
        
    Returns:
    Well, we will see...
    """

    nblist = asap3.FullNeighborList(cutoff, atoms)
    result = np.zeros(len(atoms), int)
    for i in range(len(atoms)):
        indices, relative_positions, sqdist = nblist.get_neighbors(i)
        assert(len(indices) >= 14)
        nearest = np.argsort(sqdist)[:14]
        positions = np.zeros((15,3))
        positions[1:] = relative_positions[nearest]
        data = ptmmodule.index_structure(positions)
        #struct, alloy, rmsd, rot, nb_permut = xxx
        result[i] = data[0]
    return result

if __name__ == "__main__":
    from ase.lattice.cubic import FaceCenteredCubic
    atoms = FaceCenteredCubic("Cu", size=(7,7,7), pbc=False)
    ptm = PTM(atoms)
    print ptm
    print np.bincount(ptm)
