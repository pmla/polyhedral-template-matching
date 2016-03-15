"""Python module for running PTM on an Asap list of atoms.

This module runs Polyhedral Template Matching on an Atoms object from ASE/Asap.
It requires that both ASE and Asap are installed.
"""

import asap3
import ptmmodule
import numpy as np


def PTM(atoms, target_structures=None, calculate_strains=False, cutoff=10.0):
    """Run Complex Hull Analysis on an Atoms object.

    Parameters:
    atoms: The atoms object

    target_structures=None: A tuple of structures to be investigated.
        It defaults to ('sc', 'fcc', 'hcp', 'ico', 'bcc').
        It MUST be a tuple, not a list or other sequence (bug?).
    
    calculate_strains=False: Set to True to calculate strains.

    cutoff: A cutoff used for the neighborlist.  Must be large enough that all
        nearest neighbors are returned (second-nearest for BCC).  Using a too
        large value may impact performance, but not results.
        
    Returns:
    (structures, alloytypes, rmsds, scales, rotations, [strains])
    Each of these are a NumPy array of the same length as the number of atoms.
    For each atom, the data returned is

    structures[i]: The local crystal structure around atom i, if any.
         0 = none; 1 = SC; 2 = FCC; 3 = HCP; 4 = Icosahedral; 5 = BCC.

    alloytypes[i]: The alloy structure identified.
         0 = unidentified; 1 = pure element; 2 = L1_0;
         3 = L1_2 majority atom; 4 = L1_2 minority atom.
         (0 is returned if structures[i] != 2 or if no known alloy structure
         is recognized)

    rmsds[i]: The RMSD error in the fitting to the template, or INF if
         no structure was identified.

    scales[i]: The average distance to the nearest neighbors for
         structures 1-4; or the average distance to nearest and
         next-nearest neighbors for structure 5 (BCC); or INF if no
         structure was identified.

    rotations[i]: The rotation of the crystal lattice, expressed as a
         unit quaternion.  If no structure was found, the illegal
         value (0, 0, 0, 0) is returned.

    strains[i] (only returned if calculate_strains=True).  The strain
         tensor as a symmetric 3x3 matrix.  The trace of the matrix is
         1.0, since a hydrostatic component of the strain cannot be
         determined without a-priori knowledge of the reference
         lattice parameter.  If such knowledge is available, the
         hydrostatic component of the strain can be calculated from
         scales[i].

    """

    nblist = asap3.FullNeighborList(cutoff, atoms)
    structures = np.zeros(len(atoms), int)
    alloys = np.zeros(len(atoms), int)
    rmsds = np.zeros(len(atoms))
    scales = np.zeros(len(atoms))
    rotations = np.zeros((len(atoms), 4))
    if calculate_strains:
        strains = np.zeros((len(atoms), 3, 3))
    z_all = atoms.get_atomic_numbers()
    for i in range(len(atoms)):
        indices, relative_positions, sqdist = nblist.get_neighbors(i)
        assert(len(indices) >= 14)
        nearest = np.argsort(sqdist)[:14]
        positions = np.zeros((15,3))
        positions[1:] = relative_positions[nearest]

        z = np.zeros(15, np.int32)
        z[0] = z_all[i]
        z[1:] = z_all[indices[nearest]]
        print z.dtype
        data = ptmmodule.index_structure(positions, z, 
                                         calculate_strains=calculate_strains)
        # data = (struct, alloy, rmsd, scale, rotation)
        structures[i], alloys[i], rmsds[i], scales[i] = data[:4]
        if structures[i]:
            rotations[i] = data[4]
            if calculate_strains:
                strains[i] = data[7]
    if calculate_strains:
        return structures, alloys, rmsds, scales, rotations, strains
    else:
        return structures, alloys, rmsds, scales, rotations

if __name__ == "__main__":
    from ase.lattice.cubic import FaceCenteredCubic
    atoms = FaceCenteredCubic("Cu", size=(7,7,7), pbc=False)
    structures, alloys, rmsds, scales, rotations, strains = PTM(atoms, calculate_strains=True)
    print structures
    print np.bincount(structures)
    print 
    print alloys
    print rmsds
    print scales
    print rotations
    print

    # In the following, surface atoms will have type 0, or be detected
    # as 1 or 2 with a high rmsd and a strange rotation whereas bulk
    # atoms will be correctly identified as type 2 (fcc) with very
    # small rmsd and rotation (1,0,0,0).

    for s, a, r, sc, rot in zip(structures, alloys, rmsds, scales, 
                                   rotations):
        print s, a, r, sc, rot

