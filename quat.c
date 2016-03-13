#include <string.h>
#include <math.h>


#define SQRT_2         1.4142135623730951454746218587388284504414
#define HALF_SQRT_2    0.7071067811865474617150084668537601828575

#define PHI            1.6180339887498949025257388711906969547272
#define HALF_PHI       0.8090169943749474512628694355953484773636

#define INV_PHI        0.6180339887498947915034364086750429123640
#define HALF_INV_PHI   0.3090169943749473957517182043375214561820

#define SQRT_5_        2.23606797749978969640917366873127623544061835961152572427089

double generator_cubic[24][4] = {	{1,	0,	0,	0	},
					{0,	1,	0,	0	},
					{0,	0,	1,	0	},
					{0,	0,	0,	1	},
					{0.5,	0.5,	0.5,	0.5	},
					{0.5,	0.5,	-0.5,	0.5	},
					{0.5,	-0.5,	0.5,	0.5	},
					{0.5,	-0.5,	-0.5,	0.5	},
					{-0.5,	0.5,	0.5,	0.5	},
					{-0.5,	0.5,	-0.5,	0.5	},
					{-0.5,	-0.5,	0.5,	0.5	},
					{-0.5,	-0.5,	-0.5,	0.5	},
					{HALF_SQRT_2,	HALF_SQRT_2,	0,	0	},
					{HALF_SQRT_2,	0,	HALF_SQRT_2,	0	},
					{HALF_SQRT_2,	0,	0,	HALF_SQRT_2	},
					{-HALF_SQRT_2,	HALF_SQRT_2,	0,	0	},
					{-HALF_SQRT_2,	0,	HALF_SQRT_2,	0	},
					{-HALF_SQRT_2,	0,	0,	HALF_SQRT_2	},
					{0,	HALF_SQRT_2,	HALF_SQRT_2,	0	},
					{0,	HALF_SQRT_2,	0,	HALF_SQRT_2	},
					{0,	0,	HALF_SQRT_2,	HALF_SQRT_2	},
					{0,	-HALF_SQRT_2,	HALF_SQRT_2,	0	},
					{0,	-HALF_SQRT_2,	0,	HALF_SQRT_2	},
					{0,	0,	-HALF_SQRT_2,	HALF_SQRT_2	}	};

static void quat_rot(double* r, double* a, double* b)
{
	b[0] = (r[0] * a[0] - r[1] * a[1] - r[2] * a[2] - r[3] * a[3]);
	b[1] = (r[0] * a[1] + r[1] * a[0] + r[2] * a[3] - r[3] * a[2]);
	b[2] = (r[0] * a[2] - r[1] * a[3] + r[2] * a[0] + r[3] * a[1]);
	b[3] = (r[0] * a[3] + r[1] * a[2] - r[2] * a[1] + r[3] * a[0]);
}

void rotate_quaternion_into_cubic_fundamental_zone(double* q)
{
	double max = 0.0;
	int i = 0, bi = -1;
	for (i=0;i<24;i++)
	{
		double* g = generator_cubic[i];
		double t = fabs(q[0] * g[0] - q[1] * g[1] - q[2] * g[2] - q[3] * g[3]);
		if (t > max)
		{
			max = t;
			bi = i;
		}
	}

	double f[4];
	quat_rot(q, generator_cubic[bi], f);
	memcpy(q, &f, 4 * sizeof(double));
}

