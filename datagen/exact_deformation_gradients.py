from sympy import *

fcc = [
								[          0,          0,          0 ],
								[  sqrt(2)/2,  sqrt(2)/2,          0 ],
								[          0,  sqrt(2)/2,  sqrt(2)/2 ],
								[  sqrt(2)/2,          0,  sqrt(2)/2 ],
								[ -sqrt(2)/2, -sqrt(2)/2,          0 ],
								[          0, -sqrt(2)/2, -sqrt(2)/2 ],
								[ -sqrt(2)/2,          0, -sqrt(2)/2 ],
								[ -sqrt(2)/2,  sqrt(2)/2,          0 ],
								[          0, -sqrt(2)/2,  sqrt(2)/2 ],
								[ -sqrt(2)/2,          0,  sqrt(2)/2 ],
								[  sqrt(2)/2, -sqrt(2)/2,          0 ],
								[          0,  sqrt(2)/2, -sqrt(2)/2 ],
								[  sqrt(2)/2,          0, -sqrt(2)/2 ],
]

hcp = [
								[          0,          0,          0 ],
								[        sqrt(1)/2, -sqrt(3)/2,          0 ],
								[         -1,          0,          0 ],
								[       -sqrt(1)/2,  sqrt(3)/6, -sqrt(6)/3 ],
								[        sqrt(1)/2,  sqrt(3)/6, -sqrt(6)/3 ],
								[          0, -sqrt(3)/3, -sqrt(6)/3 ],
								[       -sqrt(1)/2,  sqrt(3)/2,          0 ],
								[        sqrt(1)/2,  sqrt(3)/2,          0 ],
								[          1,          0,          0 ],
								[       -sqrt(1)/2, -sqrt(3)/2,          0 ],
								[          0, -sqrt(3)/3,  sqrt(6)/3 ],
								[        sqrt(1)/2,  sqrt(3)/6,  sqrt(6)/3 ],
								[       -sqrt(1)/2,  sqrt(3)/6,  sqrt(6)/3 ],
]

bcc = [
								[                0,                0,                0 ],
								[ 7*sqrt(3)/3-7*sqrt(1)/2, 7*sqrt(3)/3-7*sqrt(1)/2, 7*sqrt(3)/3-7*sqrt(1)/2 ],
								[ 7*sqrt(1)/2-7*sqrt(3)/3, 7*sqrt(3)/3-7*sqrt(1)/2, 7*sqrt(3)/3-7*sqrt(1)/2 ],
								[ 7*sqrt(3)/3-7*sqrt(1)/2, 7*sqrt(3)/3-7*sqrt(1)/2, 7*sqrt(1)/2-7*sqrt(3)/3 ],
								[ 7*sqrt(1)/2-7*sqrt(3)/3, 7*sqrt(1)/2-7*sqrt(3)/3, 7*sqrt(3)/3-7*sqrt(1)/2 ],
								[ 7*sqrt(3)/3-7*sqrt(1)/2, 7*sqrt(1)/2-7*sqrt(3)/3, 7*sqrt(3)/3-7*sqrt(1)/2 ],
								[ 7*sqrt(1)/2-7*sqrt(3)/3, 7*sqrt(3)/3-7*sqrt(1)/2, 7*sqrt(1)/2-7*sqrt(3)/3 ],
								[ 7*sqrt(1)/2-7*sqrt(3)/3, 7*sqrt(1)/2-7*sqrt(3)/3, 7*sqrt(1)/2-7*sqrt(3)/3 ],
								[ 7*sqrt(3)/3-7*sqrt(1)/2, 7*sqrt(1)/2-7*sqrt(3)/3, 7*sqrt(1)/2-7*sqrt(3)/3 ],
								[   14*sqrt(3)/3-7,                0,                0 ],
								[   7-14*sqrt(3)/3,                0,                0 ],
								[                0,   14*sqrt(3)/3-7,                0 ],
								[                0,   7-14*sqrt(3)/3,                0 ],
								[                0,                0,   14*sqrt(3)/3-7 ],
								[                0,                0,   7-14*sqrt(3)/3 ],
]

ico = [
								[                     0,                     0,                     0 ],
								[                     0,                     0,                     1 ],
								[                     0,                     0,                    -1 ],
								[ -sqrt((5-sqrt(5))/10),        (5+sqrt(5))/10,            -sqrt(5)/5 ],
								[  sqrt((5-sqrt(5))/10),       -(5+sqrt(5))/10,             sqrt(5)/5 ],
								[                     0,          -2*sqrt(5)/5,            -sqrt(5)/5 ],
								[                     0,           2*sqrt(5)/5,             sqrt(5)/5 ],
								[  sqrt((5+sqrt(5))/10),       -(5-sqrt(5))/10,            -sqrt(5)/5 ],
								[ -sqrt((5+sqrt(5))/10),        (5-sqrt(5))/10,             sqrt(5)/5 ],
								[ -sqrt((5+sqrt(5))/10),       -(5-sqrt(5))/10,            -sqrt(5)/5 ],
								[  sqrt((5+sqrt(5))/10),        (5-sqrt(5))/10,             sqrt(5)/5 ],
								[  sqrt((5-sqrt(5))/10),        (5+sqrt(5))/10,            -sqrt(5)/5 ],
								[ -sqrt((5-sqrt(5))/10),       -(5+sqrt(5))/10,             sqrt(5)/5 ],
]

sc = [
								[  0,  0,  0 ],
								[  0,  0, -sqrt(1) ],
								[  0,  0,  sqrt(1) ],
								[  0, -sqrt(1),  0 ],
								[  0,  sqrt(1),  0 ],
								[ -sqrt(1),  0,  0 ],
								[  sqrt(1),  0,  0 ],
]

dcub = [
								[                      0,                      0,                      0 ],
								[  4/(sqrt(3)+6*sqrt(2)),  4/(sqrt(3)+6*sqrt(2)),  4/(sqrt(3)+6*sqrt(2)) ],
								[  4/(sqrt(3)+6*sqrt(2)), -4/(sqrt(3)+6*sqrt(2)), -4/(sqrt(3)+6*sqrt(2)) ],
								[ -4/(sqrt(3)+6*sqrt(2)), -4/(sqrt(3)+6*sqrt(2)),  4/(sqrt(3)+6*sqrt(2)) ],
								[ -4/(sqrt(3)+6*sqrt(2)),  4/(sqrt(3)+6*sqrt(2)), -4/(sqrt(3)+6*sqrt(2)) ],
								[  8/(sqrt(3)+6*sqrt(2)),  8/(sqrt(3)+6*sqrt(2)),                      0 ],
								[                      0,  8/(sqrt(3)+6*sqrt(2)),  8/(sqrt(3)+6*sqrt(2)) ],
								[  8/(sqrt(3)+6*sqrt(2)),                      0,  8/(sqrt(3)+6*sqrt(2)) ],
								[                      0, -8/(sqrt(3)+6*sqrt(2)), -8/(sqrt(3)+6*sqrt(2)) ],
								[  8/(sqrt(3)+6*sqrt(2)), -8/(sqrt(3)+6*sqrt(2)),                      0 ],
								[  8/(sqrt(3)+6*sqrt(2)),                      0, -8/(sqrt(3)+6*sqrt(2)) ],
								[ -8/(sqrt(3)+6*sqrt(2)), -8/(sqrt(3)+6*sqrt(2)),                      0 ],
								[                      0, -8/(sqrt(3)+6*sqrt(2)),  8/(sqrt(3)+6*sqrt(2)) ],
								[ -8/(sqrt(3)+6*sqrt(2)),                      0,  8/(sqrt(3)+6*sqrt(2)) ],
								[ -8/(sqrt(3)+6*sqrt(2)),                      0, -8/(sqrt(3)+6*sqrt(2)) ],
								[ -8/(sqrt(3)+6*sqrt(2)),  8/(sqrt(3)+6*sqrt(2)),                      0 ],
								[                      0,  8/(sqrt(3)+6*sqrt(2)), -8/(sqrt(3)+6*sqrt(2)) ],
]

dhex = [
								[                                   0,                                   0,                                   0 ],
								[      -4*sqrt(2)/(sqrt(3)+6*sqrt(2)),   4*sqrt(6)/(3*(sqrt(3)+6*sqrt(2))),   -4*sqrt(3)/(3*sqrt(3)+18*sqrt(2)) ],
								[                                   0,   -8*sqrt(6)/(3*sqrt(3)+18*sqrt(2)),   -4*sqrt(3)/(3*sqrt(3)+18*sqrt(2)) ],
								[       4*sqrt(2)/(sqrt(3)+6*sqrt(2)),   4*sqrt(6)/(3*(sqrt(3)+6*sqrt(2))),   -4*sqrt(3)/(3*sqrt(3)+18*sqrt(2)) ],
								[                                   0,                                   0,       4*sqrt(3)/(sqrt(3)+6*sqrt(2)) ],
								[      -8*sqrt(2)/(sqrt(3)+6*sqrt(2)),                                   0,                                   0 ],
								[      -4*sqrt(2)/(sqrt(3)+6*sqrt(2)),   4*sqrt(6)/(3*(sqrt(3)+6*sqrt(2))), -16*sqrt(3)/(3*(sqrt(3)+6*sqrt(2))) ],
								[      -4*sqrt(2)/(sqrt(3)+6*sqrt(2)),       4*sqrt(6)/(sqrt(3)+6*sqrt(2)),                                   0 ],
								[       4*sqrt(2)/(sqrt(3)+6*sqrt(2)),      -4*sqrt(6)/(sqrt(3)+6*sqrt(2)),                                   0 ],
								[                                   0,   -8*sqrt(6)/(3*sqrt(3)+18*sqrt(2)), -16*sqrt(3)/(3*(sqrt(3)+6*sqrt(2))) ],
								[      -4*sqrt(2)/(sqrt(3)+6*sqrt(2)),      -4*sqrt(6)/(sqrt(3)+6*sqrt(2)),                                   0 ],
								[       4*sqrt(2)/(sqrt(3)+6*sqrt(2)),   4*sqrt(6)/(3*(sqrt(3)+6*sqrt(2))), -16*sqrt(3)/(3*(sqrt(3)+6*sqrt(2))) ],
								[       4*sqrt(2)/(sqrt(3)+6*sqrt(2)),       4*sqrt(6)/(sqrt(3)+6*sqrt(2)),                                   0 ],
								[       8*sqrt(2)/(sqrt(3)+6*sqrt(2)),                                   0,                                   0 ],
								[                                   0,   -8*sqrt(6)/(3*sqrt(3)+18*sqrt(2)),  16*sqrt(3)/(3*(sqrt(3)+6*sqrt(2))) ],
								[       4*sqrt(2)/(sqrt(3)+6*sqrt(2)),   4*sqrt(6)/(3*(sqrt(3)+6*sqrt(2))),  16*sqrt(3)/(3*(sqrt(3)+6*sqrt(2))) ],
								[      -4*sqrt(2)/(sqrt(3)+6*sqrt(2)),   4*sqrt(6)/(3*(sqrt(3)+6*sqrt(2))),  16*sqrt(3)/(3*(sqrt(3)+6*sqrt(2))) ],
]


def format(n, penrose):

	data = [[str(penrose[j,i]) for i in range(3)] for j in range(n)]
	nmax = max([max([len(e) for e in row]) for row in data])

	output = []
	for row in data:
		print "{ " + ", ".join([e.rjust(nmax) for e in row]) + " },"

def go():

	for structure in [sc, fcc, hcp, ico, bcc, dcub, dhex]:
		#M = np.dot(structure.T, structure)
		#mpi_scale = np.trace(M) / 3
		#inv = structure / mpi_scale

		m = 0
		n = len(structure)
		m = sum([structure[j][0]**2 for j in range(n)])
		m = m.expand().simplify()

		structure = Matrix(structure)
		penrose = structure / m
		for j in range(n):
			for i in range(3):
				penrose[j,i] = penrose[j,i].expand().simplify()
		format(n, penrose)
go()
