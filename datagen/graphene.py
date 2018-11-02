import numpy as np
from sympy import *
import matplotlib.pyplot as plt

def gen_template():

	first = [(np.sin(i * 2 * np.pi / 3), np.cos(i * 2 * np.pi / 3), 0) for i in range(3)]
	first = np.array(first)
	second = []
	for a in first:
		for d in first:
			b = a + (d[0], -d[1], d[2])
			if np.linalg.norm(b) > 1E-3:
				second += [b]
	angles = [np.arctan2(y, x) for x, y, _ in second]
	indices = np.argsort(angles)[::-1]
	indices = np.roll(indices, -1)
	second = np.array(second)[indices]

	template = [(0, 0, 0)] + list(first) + list(second)
	template = np.array(template)

	if 1:
		xs, ys, _ = zip(*template)
		plt.scatter(xs, ys)
		for i, (x, y) in enumerate(zip(xs, ys)):
			plt.text(x, y, str(i))
		plt.show()
		return

	#print template
	scale = np.mean(np.linalg.norm(template[1:], axis=1))
	#print scale
	scale = (3 + 6 * sqrt(3)) / 9
	#print scale.evalf()
	#print scale


	sym = Matrix([
		[         0,            0,  0],
		[         0,            1,  0],
		[ sqrt(3)/2,   -sqrt(1)/2,  0],
		[-sqrt(3)/2,   -sqrt(1)/2,  0],
		[-sqrt(3)/2,  3*sqrt(1)/2,  0],
		[ sqrt(3)/2,  3*sqrt(1)/2,  0],
		[   sqrt(3),            0,  0],
		[ sqrt(3)/2, -3*sqrt(1)/2,  0],
		[-sqrt(3)/2, -3*sqrt(1)/2,  0],
		[  -sqrt(3),            0,  0]])
	sym = sym / scale

	for i in range(len(sym)):
		sym[i] = sym[i].expand().simplify()
	print sym


	#return template / scale

gen_template()
