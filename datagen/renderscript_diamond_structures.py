import ovito
from ovito import *
from ovito.data import *
from ovito.vis import *
import numpy as np
import math
import os
import multiprocessing
import itertools


def get_templates():

	sqrt = np.sqrt

	SC = [	[  0.            ,  0.            ,  0.             ],
		[  0.            ,  0.            , -1.             ],
		[  0.            ,  0.            ,  1.             ],
		[  0.            , -1.            ,  0.             ],
		[  0.            ,  1.            ,  0.             ],
		[ -1.            ,  0.            ,  0.             ],
		[  1.            ,  0.            ,  0.             ]	]

	FCC = [	[  0.            ,  0.            ,  0.             ],
		[  0.            ,  0.707106781187,  0.707106781187 ],
		[  0.            , -0.707106781187, -0.707106781187 ],
		[  0.            ,  0.707106781187, -0.707106781187 ],
		[  0.            , -0.707106781187,  0.707106781187 ],
		[  0.707106781187,  0.            ,  0.707106781187 ],
		[ -0.707106781187,  0.            , -0.707106781187 ],
		[  0.707106781187,  0.            , -0.707106781187 ],
		[ -0.707106781187,  0.            ,  0.707106781187 ],
		[  0.707106781187,  0.707106781187,  0.             ],
		[ -0.707106781187, -0.707106781187,  0.             ],
		[  0.707106781187, -0.707106781187,  0.             ],
		[ -0.707106781187,  0.707106781187,  0.             ]	]

	HCP = [	[  0.            ,  0.            ,  0.             ],
		[  0.707106781186,  0.            ,  0.707106781186 ],
		[ -0.235702260395, -0.942809041583, -0.235702260395 ],
		[  0.707106781186,  0.707106781186,  0.             ],
		[ -0.235702260395, -0.235702260395, -0.942809041583 ],
		[  0.            ,  0.707106781186,  0.707106781186 ],
		[ -0.942809041583, -0.235702260395, -0.235702260395 ],
		[ -0.707106781186,  0.707106781186,  0.             ],
		[  0.            ,  0.707106781186, -0.707106781186 ],
		[  0.707106781186,  0.            , -0.707106781186 ],
		[  0.707106781186, -0.707106781186,  0.             ],
		[ -0.707106781186,  0.            ,  0.707106781186 ],
		[  0.            , -0.707106781186,  0.707106781186 ]	]

	ICO = [	[  0.            ,  0.            ,  0.             ],
		[  0.            ,  0.525731112119,  0.850650808352 ],
		[  0.            , -0.525731112119, -0.850650808352 ],
		[  0.            ,  0.525731112119, -0.850650808352 ],
		[  0.            , -0.525731112119,  0.850650808352 ],
		[ -0.525731112119, -0.850650808352,  0.             ],
		[  0.525731112119,  0.850650808352,  0.             ],
		[  0.525731112119, -0.850650808352,  0.             ],
		[ -0.525731112119,  0.850650808352,  0.             ],
		[ -0.850650808352,  0.            , -0.525731112119 ],
		[  0.850650808352,  0.            ,  0.525731112119 ],
		[  0.850650808352,  0.            , -0.525731112119 ],
		[ -0.850650808352,  0.            ,  0.525731112119 ]	]

	BCC = [	[  0.            ,  0.            ,  0.             ],
		[ -0.541451884327, -0.541451884327, -0.541451884327 ],
		[  0.541451884327,  0.541451884327,  0.541451884327 ],
		[  0.541451884327, -0.541451884327, -0.541451884327 ],
		[ -0.541451884327,  0.541451884327,  0.541451884327 ],
		[ -0.541451884327,  0.541451884327, -0.541451884327 ],
		[  0.541451884327, -0.541451884327,  0.541451884327 ],
		[ -0.541451884327, -0.541451884327,  0.541451884327 ],
		[  0.541451884327,  0.541451884327, -0.541451884327 ],
		[  0.            ,  0.            , -1.082903768655 ],
		[  0.            ,  0.            ,  1.082903768655 ],
		[  0.            , -1.082903768655,  0.             ],
		[  0.            ,  1.082903768655,  0.             ],
		[ -1.082903768655,  0.            ,  0.             ],
		[  1.082903768655,  0.            ,  0.             ]	]

	DCUB = [[  0.            ,  0.            ,  0.             ],
		[ -0.391491627053,  0.391491627053,  0.391491627053 ],
		[ -0.391491627053, -0.391491627053, -0.391491627053 ],
		[  0.391491627053, -0.391491627053,  0.391491627053 ],
		[  0.391491627053,  0.391491627053, -0.391491627053 ],
		[ -0.782983254107,  0.            ,  0.782983254107 ],
		[ -0.782983254107,  0.782983254107,  0.             ],
		[  0.            ,  0.782983254107,  0.782983254107 ],
		[ -0.782983254107, -0.782983254107,  0.             ],
		[ -0.782983254107,  0.            , -0.782983254107 ],
		[  0.            , -0.782983254107, -0.782983254107 ],
		[  0.            , -0.782983254107,  0.782983254107 ],
		[  0.782983254107, -0.782983254107,  0.             ],
		[  0.782983254107,  0.            ,  0.782983254107 ],
		[  0.            ,  0.782983254107, -0.782983254107 ],
		[  0.782983254107,  0.            , -0.782983254107 ],
		[  0.782983254107,  0.782983254107,  0.             ]	]

	DHEX = [[  0.            ,  0.            ,  0.             ],
		[ -0.391491627053, -0.391491627053, -0.391491627053 ],
		[  0.391491627053, -0.391491627053,  0.391491627053 ],
		[ -0.391491627053,  0.391491627053,  0.391491627053 ],
		[  0.391491627053,  0.391491627053, -0.391491627053 ],
		[ -0.260994418036, -1.043977672142, -0.260994418036 ],
		[ -1.043977672142, -0.260994418036, -0.260994418036 ],
		[ -0.260994418036, -0.260994418036, -1.043977672142 ],
		[  0.782983254107,  0.            ,  0.782983254107 ],
		[  0.782983254107, -0.782983254107,  0.             ],
		[  0.            , -0.782983254107,  0.782983254107 ],
		[  0.            ,  0.782983254107,  0.782983254107 ],
		[ -0.782983254107,  0.782983254107,  0.             ],
		[ -0.782983254107,  0.            ,  0.782983254107 ],
		[  0.782983254107,  0.782983254107,  0.             ],
		[  0.            ,  0.782983254107, -0.782983254107 ],
		[  0.782983254107,  0.            , -0.782983254107 ]	]

	k7 = 2*sqrt(3)-1
	k8 = 6-sqrt(3)
	GRP =  np.array(
	       [[     0,     0,  0 ],
                [     0,  2*k7,  0 ],
                [    k8,  - k7,  0 ],
                [   -k8,   -k7,  0 ],
                [   -k8,  3*k7,  0 ],
                [    k8,  3*k7,  0 ],
                [  2*k8,     0,  0 ],
                [    k8, -3*k7,  0 ],
                [   -k8, -3*k7,  0 ],
                [ -2*k8,     0,  0 ]	]) * 3 / 22

	return [np.array(e) for e in [SC, FCC, HCP, ICO, BCC, DCUB, DHEX, GRP]]

def find_plane(positions):

	best = (0, None)
	cs = itertools.combinations(range(len(positions)), 3)
	for c in cs:
		qs = positions[list(c)]
		vs = qs[1:] - qs[0]
		normal = np.cross(vs[0], vs[1])
		if np.linalg.norm(normal) < 1E-3:
			continue
		normal /= np.linalg.norm(normal)
		dots = np.abs(np.dot(positions - qs[0], normal))
		count = len(np.where(dots < 1E-3)[0])
		best = max(best, (count, tuple(c)))
	print(best)
	return np.array(best[1])

def project(positions, indices):

	a, b, c = positions[indices]

	x = a / np.linalg.norm(a)

	y = b / np.linalg.norm(b)
	y -= np.dot(x, y) * x
	y /= np.linalg.norm(y)

	z = np.cross(x, y)

	assert abs(np.dot(x, y)) < 1E-9
	assert abs(np.dot(y, z)) < 1E-9
	assert abs(np.dot(z, x)) < 1E-9
	basis = np.array([x, y, z])

	return np.dot(positions, basis.T)

def render_scene(folder, name):

	nameindex = dict()
	names = ['SC', 'FCC', 'HCP', 'ICO', 'BCC', 'DCUB', 'DHEX', 'GRP']
	for i, n in enumerate(names):
		nameindex[n.lower()] = i

	index = nameindex[name]

	points = get_templates()
	positions = points[index] * 16
	num_atoms = len(positions)

	if name in [e.lower() for e in ['FCC', 'HCP', 'DCUB', 'DHEX']]:
		indices = find_plane(positions)
		positions = project(positions, indices)

	distances = [np.linalg.norm(p-q) for i, p in enumerate(positions[1:]) for j, q in enumerate(positions[1:]) if i<j]
	cutoff = 1.1 * min(distances)

	# Create the particles position property.
	pos_prop = ParticleProperty.create(ParticleProperty.Type.Position, num_atoms)
	pos_prop.mutable_array[::] = positions

	# Create the particle type property and insert two atom types.
	type_prop = ParticleProperty.create(ParticleProperty.Type.ParticleType, num_atoms)
	type_prop.type_list.append(ParticleType(id = 1, name = 'Cu', color = (1.0,1.0,1.0)))
	type_prop.type_list.append(ParticleType(id = 2, name = 'Ni', color = (31/255.,120/255.,180/255.)))
	type_prop.type_list.append(ParticleType(id = 3, name = 'Pt', color = (51/255.,160/255.,44/255.)))
	type_prop.marray[0] = 1  # First atom is Cu
	for i in range(1, num_atoms):
		type_prop.marray[i] = 2  #subsequent atoms are Ni

	indices = np.where(np.linalg.norm(positions, axis=1) > np.linalg.norm(positions[1]) * 1.02)[0]
	for i in indices:
		print(i)
		type_prop.marray[i] = 3

	w = 1.3
	r = 2
	type_prop.type_list[0].radius = r
	type_prop.type_list[1].radius = r
	type_prop.type_list[2].radius = r

	xmax = np.max(positions[:,0])
	ymax = np.max(positions[:,1])
	zmax = np.max(positions[:,2])

	# Create a simulation box.
	cell = SimulationCell()
	cell.matrix = [[xmax,0,0,0], [0,ymax,0,0], [0,0,xmax,0]]
	cell.pbc = (False, False, False)
	cell.display.line_width = 0.0

	# Create a data collection to hold the particle properties, the bonds, and the simulation cell.
	data = DataCollection()
	data.add(pos_prop)
	data.add(cell)
	data.add(type_prop)


	# Create a node and insert it into the scene.
	node = ObjectNode()
	node.source = data
	modifier = ovito.modifiers.CreateBondsModifier(mode=ovito.modifiers.CreateBondsModifier.Mode.Pairwise)
	modifier.bonds_display.color = (166/255.,206/255.,227/255.)

	if name == 'grp':
		modifier.set_pairwise_cutoff(type_prop.type_list[0].name, type_prop.type_list[0].name, cutoff)
		modifier.set_pairwise_cutoff(type_prop.type_list[0].name, type_prop.type_list[1].name, cutoff)
	else:
		modifier.set_pairwise_cutoff(type_prop.type_list[0].name, type_prop.type_list[0].name, 0)
		modifier.set_pairwise_cutoff(type_prop.type_list[0].name, type_prop.type_list[1].name, 0)

	if len(indices) == 0:
		modifier.set_pairwise_cutoff(type_prop.type_list[1].name, type_prop.type_list[1].name, cutoff)
	elif len(indices) == 6:
		modifier.set_pairwise_cutoff(type_prop.type_list[1].name, type_prop.type_list[1].name, cutoff)
		modifier.set_pairwise_cutoff(type_prop.type_list[1].name, type_prop.type_list[2].name, cutoff)
	elif len(indices) > 6:
		modifier.set_pairwise_cutoff(type_prop.type_list[1].name, type_prop.type_list[1].name, cutoff)
		modifier.set_pairwise_cutoff(type_prop.type_list[1].name, type_prop.type_list[2].name, cutoff)
		modifier.set_pairwise_cutoff(type_prop.type_list[2].name, type_prop.type_list[2].name, cutoff*2)

	modifier.bonds_display.width=w
	modifier.bonds_display.use_particle_colors=False
	node.modifiers.append(modifier)
	node.compute()

	dataset.scene_nodes.append(node)
	dataset.selected_node = node		# Select the new node and adjust viewport cameras to show everything.
	#return

	vp = Viewport()
	vp.type = Viewport.Type.ORTHO
	vp.fov = 254

	vp.camera_pos = (0, 0, 0)
	vp.camera_dir = (0.8, -0.6, -0.1667)
	if name == 'grp':
		vp.camera_dir = (0, 0, -1)
	#vp.camera_pos = (0, 0, 0)
	#vp.camera_dir = (-0.3, 0.9, -0.3)
	vp.fov = 20
	if name == 'grp':
		vp.fov = 22

	size = (126, 126)
	vp.render(RenderSettings(filename = "ptm_schematic_%s.png" % name, size=size))

import sys
def go():
	if len(sys.argv) < 2:
		raise Exception("must provide a template index")

	index = int(sys.argv[1])
	if index not in range(8):
		raise Exception("oh dear")

	folder = '.'
	names = ['SC', 'FCC', 'HCP', 'ICO', 'BCC', 'DCUB', 'DHEX', 'GRP']
	for i, name in enumerate(names):
		if i != index: continue
		name = name.lower()
		print(name)
		render_scene(folder, name)
		return
		p = multiprocessing.Process(target=render_scene, args=(i, folder, name))
		p.start()
		p.join()
go()

