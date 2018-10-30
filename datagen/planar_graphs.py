import math
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection


def plot_points(points, colours=None, labels=None):

	fig = plt.figure()#figsize=(16,14))
	fig.set_tight_layout(True)
	ax = fig.add_subplot(111, projection='3d', proj_type='ortho')

	if colours is not None:
		for (index, c) in zip([0, 1, 2, 3], ['C0', 'C1', 'C2', 'C3']):
			indices = np.where(colours == index)[0]
			(xs, ys, zs) = zip(*points[indices])
			ax.scatter(xs, ys, zs, c=c)
	else:
		(xs, ys, zs) = zip(*points)
		ax.scatter(xs, ys, zs)
		#for i, e in enumerate(points):
		#	c = 'rb'[i < 3]
		#	plt.plot([0, e[0]], [0, e[1]], [0, e[2]], c=c)

	if labels is not None:
		for p, l in zip(points, labels):
			ax.text(p[0], p[1], p[2], l, size=30, color='k')

	(xs, ys, zs) = zip(*points)
	lim = max([abs(e) for e in xs+ys+zs])
	ax.set_xlim(-lim, lim)
	ax.set_ylim(-lim, lim)
	ax.set_zlim(-lim, lim)
	plt.show()

def _plot_hull(points, simplices, labels=None):

	triangles = points[simplices]

	fig = plt.figure()#figsize=(16,14))
	fig.set_tight_layout(True)
	ax = fig.add_subplot(111, projection='3d')

	if 1:
		color = (0.0, 0.0, 1., 0.2)

		tri = Poly3DCollection(triangles)
		tri.set_facecolor(color)
		ax.add_collection3d(tri)

		for t in triangles:
			(xs, ys, zs) = zip(*[t[0], t[1], t[2], t[0]])
			ax.plot(xs, ys, zs, c='k')

	if labels is not None:
		for p, l in zip(points, labels):
			ax.text(p[0], p[1], p[2], l, size=30, color='k')

	(xs, ys, zs) = zip(*points)
	lim = max([abs(e) for e in xs+ys+zs])
	#ax.set_xlim(-lim, lim)
	#ax.set_ylim(-lim, lim)
	#ax.set_zlim(-lim, lim)
	plt.show()
