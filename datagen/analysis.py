import sympy
import numpy as np


sqrt = np.sqrt
array = np.array

sqrt = sympy.sqrt
array = sympy.Matrix

def norm(x):
	return sqrt(sum([e**2 for e in x]))

def go():

	p = array([-sqrt(2)/4, sqrt(3) / sqrt(2)/6, -sqrt(3)/12])
	p = array([0, 0, sqrt(3)/4])
	p = array([0, -sqrt(3)/sqrt(2)/3, -sqrt(3)/12])
	q = array([sqrt(2)/4, -sqrt(6)/4, 0])
	q = array([0, 0, -sqrt(3)/4])
	q = array([sqrt(2.0)/4.0, -sqrt(6.0)/4.0, 0.0])
	q = array([-sqrt(2.0)/2.0, 0.0, 0.0])
	q = array([-sqrt(2.0)/4.0, sqrt(6.0)/4.0, 0.0])
	q = array([sqrt(2.0)/4.0, sqrt(6.0)/4.0, 0.0])
	q = array([sqrt(2.0)/2.0, 0.0, 0.0])
	q = array([-sqrt(2.0)/4.0, -sqrt(6.0)/4.0, 0.0])
	q = array([-sqrt(2.0)/4.0, sqrt(6.0)/12.0, -sqrt(3.0)/3.0])
	q = array([sqrt(2.0)/4.0, sqrt(6.0)/12.0, -sqrt(3.0)/3.0])
	q = array([0.0, -sqrt(6.0)/6.0, -sqrt(3.0)/3.0])
	q = array([0.0, -sqrt(6.0)/6.0, sqrt(3.0)/3.0])
	q = array([sqrt(2.0)/4.0, sqrt(6.0)/12.0, sqrt(3.0)/3.0])
	q = array([-sqrt(2.0)/4.0, sqrt(6.0)/12.0, sqrt(3.0)/3.0])

	c = sqrt(3) / norm(p)

	r = sqrt(8)/sqrt(3)
	k = (4 * norm(c * p) + 12 * norm(c * r * p)) / 16

	v = [(e).expand().simplify() for e in c * p / k]
	v = [(e).expand().simplify() for e in r * c * q / k]
	print v
	print [e.evalf() for e in v]

def go():

	q = array([4*sqrt(2)/(sqrt(3) + 6*sqrt(2)),  4*sqrt(6)/(3*(sqrt(3) + 6*sqrt(2))),  16*sqrt(3)/(3*(sqrt(3) + 6*sqrt(2)))])
	v = sqrt(norm(q)**2 - norm([4*sqrt(2)/(sqrt(3) + 6*sqrt(2))])**2)
	print (v.expand().simplify())

def go():

	sqrt = np.sqrt
	dhex = np.array([
	[                                   0,                                   0,                                   0 ],
	[      -4*sqrt(2)/(sqrt(3)+6*sqrt(2)),   4*sqrt(6)/(3*(sqrt(3)+6*sqrt(2))),    4*sqrt(3)/(3*sqrt(3)+18*sqrt(2)) ],
	[                                   0,                                   0,      -4*sqrt(3)/(sqrt(3)+6*sqrt(2)) ],
	[       4*sqrt(2)/(sqrt(3)+6*sqrt(2)),   4*sqrt(6)/(3*(sqrt(3)+6*sqrt(2))),    4*sqrt(3)/(3*sqrt(3)+18*sqrt(2)) ],
	[                                   0,   -8*sqrt(6)/(3*sqrt(3)+18*sqrt(2)),    4*sqrt(3)/(3*sqrt(3)+18*sqrt(2)) ],
	[      -8*sqrt(2)/(sqrt(3)+6*sqrt(2)),                                   0,                                   0 ],
	[      -4*sqrt(2)/(sqrt(3)+6*sqrt(2)),       4*sqrt(6)/(sqrt(3)+6*sqrt(2)),                                   0 ],
	[      -4*sqrt(2)/(sqrt(3)+6*sqrt(2)),   4*sqrt(6)/(3*(sqrt(3)+6*sqrt(2))),  16*sqrt(3)/(3*(sqrt(3)+6*sqrt(2))) ],
	[      -4*sqrt(2)/(sqrt(3)+6*sqrt(2)),   4*sqrt(6)/(3*(sqrt(3)+6*sqrt(2))), -16*sqrt(3)/(3*(sqrt(3)+6*sqrt(2))) ],
	[       4*sqrt(2)/(sqrt(3)+6*sqrt(2)),   4*sqrt(6)/(3*(sqrt(3)+6*sqrt(2))), -16*sqrt(3)/(3*(sqrt(3)+6*sqrt(2))) ],
	[                                   0,   -8*sqrt(6)/(3*sqrt(3)+18*sqrt(2)), -16*sqrt(3)/(3*(sqrt(3)+6*sqrt(2))) ],
	[       4*sqrt(2)/(sqrt(3)+6*sqrt(2)),       4*sqrt(6)/(sqrt(3)+6*sqrt(2)),                                   0 ],
	[       8*sqrt(2)/(sqrt(3)+6*sqrt(2)),                                   0,                                   0 ],
	[       4*sqrt(2)/(sqrt(3)+6*sqrt(2)),   4*sqrt(6)/(3*(sqrt(3)+6*sqrt(2))),  16*sqrt(3)/(3*(sqrt(3)+6*sqrt(2))) ],
	[       4*sqrt(2)/(sqrt(3)+6*sqrt(2)),      -4*sqrt(6)/(sqrt(3)+6*sqrt(2)),                                   0 ],
	[      -4*sqrt(2)/(sqrt(3)+6*sqrt(2)),      -4*sqrt(6)/(sqrt(3)+6*sqrt(2)),                                   0 ],
	[                                   0,   -8*sqrt(6)/(3*sqrt(3)+18*sqrt(2)),  16*sqrt(3)/(3*(sqrt(3)+6*sqrt(2))) ],
	])

	hcp = np.array([
	[          0,          0,          0 ],
	[        0.5, -sqrt(3)/2,          0 ],
	[         -1,          0,          0 ],
	[       -0.5,  sqrt(3)/6, -sqrt(6)/3 ],
	[        0.5,  sqrt(3)/6, -sqrt(6)/3 ],
	[          0, -sqrt(3)/3, -sqrt(6)/3 ],
	[       -0.5,  sqrt(3)/2,          0 ],
	[        0.5,  sqrt(3)/2,          0 ],
	[          1,          0,          0 ],
	[       -0.5, -sqrt(3)/2,          0 ],
	[          0, -sqrt(3)/3,  sqrt(6)/3 ],
	[        0.5,  sqrt(3)/6,  sqrt(6)/3 ],
	[       -0.5,  sqrt(3)/6,  sqrt(6)/3 ]
	])

	alex = np.array([
		[-sqrt(2.0)/4, sqrt(3.0/2.0)/6, -sqrt(3.0)/12],
		[0, -sqrt(3.0/2.0)/3, -sqrt(3.0)/12],
		[sqrt(2.0)/4, sqrt(3.0/2.0)/6, -sqrt(3.0)/12],
		[0, 0, sqrt(3.0)/4],
		[sqrt(2.0)/4.0, -sqrt(6.0)/4.0, 0.0],
		[-sqrt(2.0)/2.0, 0.0, 0.0],
		[-sqrt(2.0)/4.0, sqrt(6.0)/4.0, 0.0],
		[sqrt(2.0)/4.0, sqrt(6.0)/4.0, 0.0],
		[sqrt(2.0)/2.0, 0.0, 0.0],
		[-sqrt(2.0)/4.0, -sqrt(6.0)/4.0, 0.0],
		[-sqrt(2.0)/4.0, sqrt(6.0)/12.0, -sqrt(3.0)/3.0],
		[sqrt(2.0)/4.0, sqrt(6.0)/12.0, -sqrt(3.0)/3.0],
		[0.0, -sqrt(6.0)/6.0, -sqrt(3.0)/3.0],
		[0.0, -sqrt(6.0)/6.0, sqrt(3.0)/3.0],
		[sqrt(2.0)/4.0, sqrt(6.0)/12.0, sqrt(3.0)/3.0],
		[-sqrt(2.0)/4.0, sqrt(6.0)/12.0, sqrt(3.0)/3.0],
	])

	better = np.array([
[0., 0., 0.],
[-0.55365277, 0.31965157, -0.2260278],
[0., -0.63930315, -0.2260278],
[0.55365277, 0.31965157, -0.2260278],
[0., 0., 0.67808339],
[-1.10730554, 0., 0.],
[-0.55365277, 0.31965157, -0.90411119],
[-0.55365277, 0.95895472, 0.],
[0.55365277, -0.95895472, 0.],
[0., -0.63930315, -0.90411119],
[-0.55365277, -0.95895472, 0.],
[0.55365277, 0.31965157, -0.90411119],
[0.55365277, 0.95895472, 0.],
[1.10730554, 0., 0.],
[0., -0.63930315, 0.90411119],
[0.55365277, 0.31965157, 0.90411119],
[-0.55365277, 0.31965157, 0.90411119]])

	import planar_graphs
	planar_graphs.plot_points(better)
go()
