import numpy as np

generator_laue_C1 = [ [1.0, 0, 0, 0] ]

generator_laue_C2 = [ [1.0, 0, 0, 0],
                      [0, 0, 0, 1.0] ]

generator_laue_C3 = [ [1.0, 0, 0, 0],
                      [0.5, 0, 0, np.sqrt(3)/2],
                      [0.5, 0, 0, -np.sqrt(3)/2] ]

generator_laue_C4 = [ [1.0, 0, 0, 0],
                      [np.sqrt(2)/2, 0, 0, np.sqrt(2)/2],
                      [0, 0, 0, 1.0],
                      [np.sqrt(2)/2, 0, 0, -np.sqrt(2)/2] ]

generator_laue_C6 = [ [1.0, 0, 0, 0],
                      [np.sqrt(3)/2, 0, 0, 0.5],
                      [0.5, 0, 0, np.sqrt(3)/2],
                      [0, 0, 0, 1.0],
                      [0.5, 0, 0, -np.sqrt(3)/2],
                      [np.sqrt(3)/2, 0, 0, -0.5] ]

generator_laue_D2 = [ [1.0, 0, 0, 0],
                      [0, 0, 0, 1.0],
                      [0, 1.0, 0, 0],
                      [0, 0, 1.0, 0] ]

generator_laue_D3 = [ [1.0, 0, 0, 0],
                      [0.5, 0, 0, np.sqrt(3)/2],
                      [0, 1.0, 0, 0],
                      [0.5, 0, 0, -np.sqrt(3)/2],
                      [0, -0.5, np.sqrt(3)/2, 0],
                      [0, 0.5, np.sqrt(3)/2, 0] ]

generator_laue_D4 = [ [1.0, 0, 0, 0],
                      [np.sqrt(2)/2, 0, 0, np.sqrt(2)/2],
                      [0, 1.0, 0, 0],
                      [0, 0, 0, 1.0],
                      [0, -np.sqrt(2)/2, np.sqrt(2)/2, 0],
                      [np.sqrt(2)/2, 0, 0, -np.sqrt(2)/2],
                      [0, 0, 1.0, 0],
                      [0, np.sqrt(2)/2, np.sqrt(2)/2, 0] ]

generator_laue_D6 = [ [1.0, 0, 0, 0],
                      [np.sqrt(3)/2, 0, 0, 0.5],
                      [0, 1.0, 0, 0],
                      [0.5, 0, 0, np.sqrt(3)/2],
                      [0, -np.sqrt(3)/2, 0.5, 0],
                      [0, 0, 0, 1.0],
                      [0, -0.5, np.sqrt(3)/2, 0],
                      [0.5, 0, 0, -np.sqrt(3)/2],
                      [0, 0, 1.0, 0],
                      [np.sqrt(3)/2, 0, 0, -0.5],
                      [0, 0.5, np.sqrt(3)/2, 0],
                      [0, np.sqrt(3)/2, 0.5, 0] ]

generator_laue_T = [ [1.0, 0, 0, 0],
                     [0, 0, 0, 1.0],
                     [0.5, 0.5, 0.5, 0.5],
                     [0.5, -0.5, 0.5, -0.5],
                     [0.5, -0.5, -0.5, -0.5],
                     [0.5, 0.5, 0.5, -0.5],
                     [0, 1.0, 0, 0],
                     [0.5, 0.5, -0.5, 0.5],
                     [0, 0, 1.0, 0],
                     [0.5, 0.5, -0.5, -0.5],
                     [0.5, -0.5, 0.5, 0.5],
                     [0.5, -0.5, -0.5, 0.5] ]

generator_cubic = [ [  1,   0,   0,   0],
[  np.sqrt(2)/2,   np.sqrt(2)/2,   0,   0],
[  np.sqrt(2)/2,   0,   np.sqrt(2)/2,   0],
[  np.sqrt(2)/2,   0,   0,   np.sqrt(2)/2],
[  np.sqrt(2)/2,   0,   0,  -np.sqrt(2)/2],
[  np.sqrt(2)/2,   0,  -np.sqrt(2)/2,   0],
[  np.sqrt(2)/2,  -np.sqrt(2)/2,  -0,  -0],
[  0.5,   0.5,   0.5,   0.5],
[  0.5,   0.5,   0.5,  -0.5],
[  0.5,   0.5,  -0.5,   0.5],
[  0.5,   0.5,  -0.5,  -0.5],
[  0.5,  -0.5,   0.5,   0.5],
[  0.5,  -0.5,   0.5,  -0.5],
[  0.5,  -0.5,  -0.5,   0.5],
[  0.5,  -0.5,  -0.5,  -0.5],
[  0,   1,   0,   0],
[  0,   np.sqrt(2)/2,   np.sqrt(2)/2,   0],
[  0,   np.sqrt(2)/2,   0,   np.sqrt(2)/2],
[  0,   np.sqrt(2)/2,   0,  -np.sqrt(2)/2],
[  0,   np.sqrt(2)/2,  -np.sqrt(2)/2,   0],
[  0,   0,   1,   0],
[  0,   0,   np.sqrt(2)/2,   np.sqrt(2)/2],
[  0,   0,   np.sqrt(2)/2,  -np.sqrt(2)/2],
[  0,   0,   0,   1]]

generator_hcp = [
	[ 1.00000000000000,  0.00000000000000,  0.00000000000000,  0.00000000000000],
	[ 0.00000000000000,  0.00000000000000,  1.00000000000000,  0.00000000000000],
	[ 0.50000000000000,  0.00000000000000,  0.00000000000000,  0.86602540378444],
	[-0.50000000000000,  0.00000000000000,  0.00000000000000,  0.86602540378444],
	[ 0.00000000000000,  0.86602540378444,  0.50000000000000,  0.00000000000000],
	[ 0.00000000000000,  0.86602540378444, -0.50000000000000,  0.00000000000000],
]

generator_ico = [
	[  1.00000000000000,   0.00000000000000,   0.00000000000000,   0.00000000000000],
	[  0.80901699437495,   0.50000000000000,   0.16245984811645,   0.26286555605957],
	[  0.80901699437495,   0.50000000000000,  -0.16245984811645,  -0.26286555605957],
	[  0.80901699437495,   0.30901699437495,   0.42532540417602,  -0.26286555605957],
	[  0.80901699437495,   0.30901699437495,  -0.42532540417602,   0.26286555605957],
	[  0.80901699437495,   0.00000000000000,   0.52573111211913,   0.26286555605957],
	[  0.80901699437495,   0.00000000000000,   0.00000000000000,   0.58778525229247],
	[  0.80901699437495,   0.00000000000000,   0.00000000000000,  -0.58778525229247],
	[  0.80901699437495,   0.00000000000000,  -0.52573111211913,  -0.26286555605957],
	[  0.80901699437495,  -0.30901699437495,   0.42532540417602,  -0.26286555605957],
	[  0.80901699437495,  -0.30901699437495,  -0.42532540417602,   0.26286555605957],
	[  0.80901699437495,  -0.50000000000000,   0.16245984811645,   0.26286555605957],
	[  0.80901699437495,  -0.50000000000000,  -0.16245984811645,  -0.26286555605957],
	[  0.50000000000000,   0.80901699437495,   0.26286555605957,  -0.16245984811645],
	[  0.50000000000000,   0.80901699437495,  -0.26286555605957,   0.16245984811645],
	[  0.50000000000000,   0.50000000000000,   0.68819096023559,   0.16245984811645],
	[  0.50000000000000,   0.50000000000000,   0.16245984811645,  -0.68819096023559],
	[  0.50000000000000,   0.50000000000000,  -0.16245984811645,   0.68819096023559],
	[  0.50000000000000,   0.50000000000000,  -0.68819096023559,  -0.16245984811645],
	[  0.50000000000000,   0.30901699437495,   0.42532540417602,   0.68819096023559],
	[  0.50000000000000,   0.30901699437495,  -0.42532540417602,  -0.68819096023559],
	[  0.50000000000000,   0.00000000000000,   0.85065080835204,  -0.16245984811645],
	[  0.50000000000000,  -0.00000000000000,   0.52573111211913,  -0.68819096023559],
	[  0.50000000000000,   0.00000000000000,  -0.52573111211913,   0.68819096023559],
	[  0.50000000000000,  -0.00000000000000,  -0.85065080835204,   0.16245984811645],
	[  0.50000000000000,  -0.30901699437495,   0.42532540417602,   0.68819096023559],
	[  0.50000000000000,  -0.30901699437495,  -0.42532540417602,  -0.68819096023559],
	[  0.50000000000000,  -0.50000000000000,   0.68819096023559,   0.16245984811645],
	[  0.50000000000000,  -0.50000000000000,   0.16245984811645,  -0.68819096023559],
	[  0.50000000000000,  -0.50000000000000,  -0.16245984811645,   0.68819096023559],
	[  0.50000000000000,  -0.50000000000000,  -0.68819096023559,  -0.16245984811645],
	[  0.50000000000000,  -0.80901699437495,   0.26286555605957,  -0.16245984811645],
	[  0.50000000000000,  -0.80901699437495,  -0.26286555605957,   0.16245984811645],
	[  0.30901699437495,   0.80901699437495,   0.26286555605957,   0.42532540417602],
	[  0.30901699437495,   0.80901699437495,  -0.26286555605957,  -0.42532540417602],
	[  0.30901699437495,   0.50000000000000,   0.68819096023559,  -0.42532540417602],
	[  0.30901699437495,   0.50000000000000,  -0.68819096023559,   0.42532540417602],
	[  0.30901699437495,   0.00000000000000,   0.85065080835204,   0.42532540417602],
	[  0.30901699437495,   0.00000000000000,   0.00000000000000,   0.95105651629515],
	[  0.30901699437495,  -0.00000000000000,  -0.00000000000000,  -0.95105651629515],
	[  0.30901699437495,  -0.00000000000000,  -0.85065080835204,  -0.42532540417602],
	[  0.30901699437495,  -0.50000000000000,   0.68819096023559,  -0.42532540417602],
	[  0.30901699437495,  -0.50000000000000,  -0.68819096023559,   0.42532540417602],
	[  0.30901699437495,  -0.80901699437495,   0.26286555605957,   0.42532540417602],
	[  0.30901699437495,  -0.80901699437495,  -0.26286555605957,  -0.42532540417602],
	[  0.00000000000000,   1.00000000000000,   0.00000000000000,   0.00000000000000],
	[  0.00000000000000,   0.80901699437495,   0.58778525229247,   0.00000000000000],
	[  0.00000000000000,   0.80901699437495,   0.26286555605957,  -0.52573111211913],
	[  0.00000000000000,   0.80901699437495,  -0.26286555605957,   0.52573111211913],
	[  0.00000000000000,   0.80901699437495,  -0.58778525229247,   0.00000000000000],
	[  0.00000000000000,   0.50000000000000,   0.68819096023559,   0.52573111211913],
	[  0.00000000000000,   0.50000000000000,   0.16245984811645,   0.85065080835204],
	[ -0.00000000000000,   0.50000000000000,  -0.16245984811645,  -0.85065080835204],
	[ -0.00000000000000,   0.50000000000000,  -0.68819096023559,  -0.52573111211913],
	[  0.00000000000000,   0.30901699437495,   0.95105651629515,   0.00000000000000],
	[ -0.00000000000000,   0.30901699437495,   0.42532540417602,  -0.85065080835204],
	[  0.00000000000000,   0.30901699437495,  -0.42532540417602,   0.85065080835204],
	[ -0.00000000000000,   0.30901699437495,  -0.95105651629515,  -0.00000000000000],
	[  0.00000000000000,   0.00000000000000,   0.85065080835204,  -0.52573111211913],
	[  0.00000000000000,   0.00000000000000,   0.52573111211913,   0.85065080835204],
]

generator_dcub = [
	[  1.00000000000000,   0.00000000000000,   0.00000000000000,   0.00000000000000],
	[  0.50000000000000,   0.50000000000000,   0.50000000000000,   0.50000000000000],
	[  0.50000000000000,   0.50000000000000,   0.50000000000000,  -0.50000000000000],
	[  0.50000000000000,   0.50000000000000,  -0.50000000000000,   0.50000000000000],
	[  0.50000000000000,   0.50000000000000,  -0.50000000000000,  -0.50000000000000],
	[  0.50000000000000,  -0.50000000000000,   0.50000000000000,   0.50000000000000],
	[  0.50000000000000,  -0.50000000000000,   0.50000000000000,  -0.50000000000000],
	[  0.50000000000000,  -0.50000000000000,  -0.50000000000000,   0.50000000000000],
	[  0.50000000000000,  -0.50000000000000,  -0.50000000000000,  -0.50000000000000],
	[  0.00000000000000,   1.00000000000000,   0.00000000000000,   0.00000000000000],
	[  0.00000000000000,   0.00000000000000,   1.00000000000000,   0.00000000000000],
	[  0.00000000000000,   0.00000000000000,   0.00000000000000,   1.00000000000000],
]

generator_dhex = [
	[  1.00000000000000,   0.00000000000000,   0.00000000000000,   0.00000000000000],
	[  0.50000000000000,   0.00000000000000,   0.00000000000000,   0.86602540378444],
	[  0.50000000000000,  -0.00000000000000,  -0.00000000000000,  -0.86602540378444],
]
