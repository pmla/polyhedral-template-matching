from sympy import *

m = Matrix([
[sqrt(2)/4, -sqrt(6)/4, 0],
[-sqrt(2)/2, 0, 0],
[-sqrt(2)/4, sqrt(6)/12, -sqrt(3)/3],
[sqrt(2)/4, sqrt(6)/12, -sqrt(3)/3],
[0, -sqrt(6)/6, -sqrt(3)/3],
[-sqrt(2)/4, sqrt(6)/4, 0],
[sqrt(2)/4, sqrt(6)/4, 0],
[sqrt(2)/2, 0, 0],
[-sqrt(2)/4, -sqrt(6)/4, 0],
[0, -sqrt(6)/6, sqrt(3)/3],
[sqrt(2)/4, sqrt(6)/12, sqrt(3)/3],
[-sqrt(2)/4, sqrt(6)/12, sqrt(3)/3],

[0, sqrt(6)/6, sqrt(3)/3],
[-sqrt(2)/4, -sqrt(6)/12, -sqrt(3)/3],
[sqrt(2)/4, -sqrt(6)/12, sqrt(3)/3],
[0, sqrt(6)/6, -sqrt(3)/3],
[sqrt(2)/4, -sqrt(6)/12, -sqrt(3)/3],
[-sqrt(2)/4, -sqrt(6)/12, sqrt(3)/3],
])

m = 12 * m
for j in range(12):

	print "[", ",".join([str(m[j,i]).rjust(12) for i in range(3)]), "],"
