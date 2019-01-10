import numpy as np
import matplotlib.pyplot as plt
from ase.spacegroup import crystal
from ase.neighborlist import NeighborList
from ase.visualize import view



def get_distances(atoms, i):

	rs = 4 * np.ones(len(atoms))
	nl = NeighborList(rs, skin=0, self_interaction=False, bothways=True)
	nl.update(atoms)

	bonds = []
	p = atoms.positions[i]
	indices, offsets = nl.get_neighbors(i)

	nbrpos = []
	for j, offset in zip(indices, offsets):
		q = atoms.positions[j] + np.dot(offset, atoms.get_cell())
		d = np.linalg.norm(p - q)
		bonds += [(d, i, j, tuple(offset))]
		nbrpos += [q - atoms.positions[i]]

	ds = np.array([e[0] for e in bonds])
	indices = np.argsort(ds)
	ds = ds[indices]
	nbrpos = np.array(nbrpos)[indices]
	return ds, nbrpos

def go():
	a = 4.65327231
	c = 2.96920288
	rutile = crystal(['H', 'H', 'O', 'O', 'O', 'O'],
			basis=[(0, 0, 0), (0.5, 0.5, 0.5),
			(0.19542, 0.80458, 0.5), (0.80458, 0.19542, 0.5), (0.69542, 0.69542, 0), (0.30458, 0.30458, 0)],
			spacegroup=136, cellpar=[a, a, c, 90, 90, 90])

	a = 3.90
	fluorite = crystal(['H', 'F', 'F'], basis=[(0, 0, 0), (0.25, 0.25, 0.25), (0.75, 0.75, 0.75)],
			spacegroup=225, cellpar=[a, a, a, 90, 90, 90])

	diamond = crystal(['Si'] * 8, basis=[[0.  , 0.  , 0.5 ],
       [0.75, 0.25, 0.25],
       [0.  , 0.5 , 0.  ],
       [0.5 , 0.5 , 0.5 ],
       [0.25, 0.25, 0.75],
       [0.5 , 0.  , 0.  ],
       [0.75, 0.75, 0.75],
       [0.25, 0.75, 0.25]], cellpar=[a, a, a, 90, 90, 90])

	atoms = rutile
	view(atoms)
	return
	print atoms.numbers

	ds, nbrpos = get_distances(atoms, i=2)
	print ds[:4]
	#print np.round(nbrpos[:22] / nbrpos[0,0]).astype(np.int)

	ds = ds[:60]
	print [(i + 1, e / ds[0]) for i, e in enumerate(ds[1:] - ds[:-1]) if e > 1E-3]
	plt.scatter(range(1, len(ds) + 1), ds)
	plt.xlim(0, len(ds) + 1)
	plt.ylim(0, round(max(ds) + 1))
	plt.show()
go()
