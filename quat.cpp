#include <cstring>
#include <cmath>


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
					{sqrt(2) / 2,	sqrt(2) / 2,	0,	0	},
					{sqrt(2) / 2,	0,	sqrt(2) / 2,	0	},
					{sqrt(2) / 2,	0,	0,	sqrt(2) / 2	},
					{-sqrt(2) / 2,	sqrt(2) / 2,	0,	0	},
					{-sqrt(2) / 2,	0,	sqrt(2) / 2,	0	},
					{-sqrt(2) / 2,	0,	0,	sqrt(2) / 2	},
					{0,	sqrt(2) / 2,	sqrt(2) / 2,	0	},
					{0,	sqrt(2) / 2,	0,	sqrt(2) / 2	},
					{0,	0,	sqrt(2) / 2,	sqrt(2) / 2	},
					{0,	-sqrt(2) / 2,	sqrt(2) / 2,	0	},
					{0,	-sqrt(2) / 2,	0,	sqrt(2) / 2	},
					{0,	0,	-sqrt(2) / 2,	sqrt(2) / 2	}	};

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

