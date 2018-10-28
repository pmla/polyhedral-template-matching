#include <string.h>
#include <cmath>
#include <cfloat>


namespace ptm {

#define SIGN(x) (x >= 0 ? 1 : -1)
#define MIN(X, Y) (((X) < (Y)) ? (X) : (Y))
#define MAX(X, Y) (((X) > (Y)) ? (X) : (Y))


#define SQRT_2         1.4142135623730951454746218587388284504414
#define HALF_SQRT_2    0.7071067811865474617150084668537601828575

#define PHI            1.6180339887498949025257388711906969547272
#define HALF_PHI       0.8090169943749474512628694355953484773636

#define INV_PHI        0.6180339887498947915034364086750429123640
#define HALF_INV_PHI   0.3090169943749473957517182043375214561820

#define SQRT_2_3       0.8164965809277260344600790631375275552273
#define SQRT_1_6       0.4082482904638630172300395315687637776136

#define SQRT_3_4       0.8660254037844385965883020617184229195118
#define SQRT_1_3       0.5773502691896257310588680411456152796745
#define SQRT_1_12      0.2886751345948128655294340205728076398373



double generator_cubic[24][4] = {		{1,	0,	0,	0	},
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

double generator_diamond_cubic[12][4] = {	{1,	0,	0,	0	},
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
						{-0.5,	-0.5,	-0.5,	0.5	}	};

double generator_hcp[6][4] = {			{1, 0, 0, 0},
						{0.5, 0.5, 0.5, 0.5},
						{0.5, -0.5, -0.5, -0.5},
						{0, SQRT_2_3, -SQRT_1_6, -SQRT_1_6},
						{0, SQRT_1_6, -SQRT_2_3, SQRT_1_6},
						{0, SQRT_1_6, SQRT_1_6, -SQRT_2_3}	};

double generator_hcp_crystalline[12][4] = {	{1, 0, 0, 0},
						{0.5, 0.5, 0.5, 0.5},
						{0.5, -0.5, -0.5, -0.5},
						{0, SQRT_2_3, -SQRT_1_6, -SQRT_1_6},
						{0, SQRT_1_6, -SQRT_2_3, SQRT_1_6},
						{0, SQRT_1_6, SQRT_1_6, -SQRT_2_3},
						{0, -SQRT_1_3, -SQRT_1_3, -SQRT_1_3},
						{SQRT_3_4, SQRT_1_12, SQRT_1_12, SQRT_1_12},
						{SQRT_3_4, -SQRT_1_12, -SQRT_1_12, -SQRT_1_12},
						{0, 1 / SQRT_2, -1 / SQRT_2, 0},
						{0, 0, 1 / SQRT_2, -1 / SQRT_2},
						{0, 1 / SQRT_2, -0, -1 / SQRT_2}	};

double generator_diamond_hexagonal[3][4] = {	{1, 0, 0, 0},
						{0.5, 0.5, 0.5, 0.5},
						{0.5, -0.5, -0.5, -0.5}	};

double generator_icosahedral[60][4] = {		{1, 0, 0, 0},
						{HALF_PHI, -HALF_INV_PHI, -0.5, 0},
						{HALF_PHI, 0, -HALF_INV_PHI, -0.5},
						{HALF_PHI, -0.5, 0, -HALF_INV_PHI},
						{HALF_PHI, HALF_INV_PHI, -0.5, 0},
						{HALF_PHI, 0, HALF_INV_PHI, -0.5},
						{HALF_PHI, -0.5, 0, HALF_INV_PHI},
						{HALF_PHI, 0.5, 0, -HALF_INV_PHI},
						{HALF_PHI, 0, -HALF_INV_PHI, 0.5},
						{HALF_PHI, -HALF_INV_PHI, 0.5, 0},
						{HALF_PHI, 0, HALF_INV_PHI, 0.5},
						{HALF_PHI, HALF_INV_PHI, 0.5, 0},
						{HALF_PHI, 0.5, 0, HALF_INV_PHI},
						{0.5, HALF_PHI, -HALF_INV_PHI, 0},
						{0.5, HALF_PHI, HALF_INV_PHI, 0},
						{0.5, 0.5, 0.5, 0.5},
						{0.5, 0.5, 0.5, -0.5},
						{0.5, 0.5, -0.5, 0.5},
						{0.5, 0.5, -0.5, -0.5},
						{0.5, HALF_INV_PHI, 0, HALF_PHI},
						{0.5, HALF_INV_PHI, 0, -HALF_PHI},
						{0.5, 0, HALF_PHI, -HALF_INV_PHI},
						{0.5, 0, HALF_PHI, HALF_INV_PHI},
						{0.5, 0, -HALF_PHI, -HALF_INV_PHI},
						{0.5, 0, -HALF_PHI, HALF_INV_PHI},
						{0.5, -HALF_INV_PHI, 0, HALF_PHI},
						{0.5, -HALF_INV_PHI, 0, -HALF_PHI},
						{0.5, -0.5, 0.5, 0.5},
						{0.5, -0.5, 0.5, -0.5},
						{0.5, -0.5, -0.5, 0.5},
						{0.5, -0.5, -0.5, -0.5},
						{0.5, -HALF_PHI, -HALF_INV_PHI, 0},
						{0.5, -HALF_PHI, HALF_INV_PHI, 0},
						{HALF_INV_PHI, -HALF_PHI, 0, -0.5},
						{HALF_INV_PHI, 0, -0.5, -HALF_PHI},
						{HALF_INV_PHI, -0.5, -HALF_PHI, 0},
						{HALF_INV_PHI, 0, 0.5, -HALF_PHI},
						{HALF_INV_PHI, -HALF_PHI, 0, 0.5},
						{HALF_INV_PHI, 0.5, -HALF_PHI, 0},
						{HALF_INV_PHI, HALF_PHI, 0, -0.5},
						{HALF_INV_PHI, -0.5, HALF_PHI, 0},
						{HALF_INV_PHI, 0, -0.5, HALF_PHI},
						{HALF_INV_PHI, HALF_PHI, 0, 0.5},
						{HALF_INV_PHI, 0, 0.5, HALF_PHI},
						{HALF_INV_PHI, 0.5, HALF_PHI, 0},
						{0, 1, 0, 0},
						{0, HALF_PHI, -0.5, HALF_INV_PHI},
						{0, HALF_PHI, -0.5, -HALF_INV_PHI},
						{0, HALF_PHI, 0.5, HALF_INV_PHI},
						{0, HALF_PHI, 0.5, -HALF_INV_PHI},
						{0, 0.5, HALF_INV_PHI, -HALF_PHI},
						{0, 0.5, HALF_INV_PHI, HALF_PHI},
						{0, 0.5, -HALF_INV_PHI, -HALF_PHI},
						{0, 0.5, -HALF_INV_PHI, HALF_PHI},
						{0, HALF_INV_PHI, -HALF_PHI, 0.5},
						{0, HALF_INV_PHI, -HALF_PHI, -0.5},
						{0, HALF_INV_PHI, HALF_PHI, 0.5},
						{0, HALF_INV_PHI, HALF_PHI, -0.5},
						{0, 0, 1, 0},
						{0, 0, 0, 1}	};

static void quat_rot(double* r, double* a, double* b)
{
	b[0] = (r[0] * a[0] - r[1] * a[1] - r[2] * a[2] - r[3] * a[3]);
	b[1] = (r[0] * a[1] + r[1] * a[0] + r[2] * a[3] - r[3] * a[2]);
	b[2] = (r[0] * a[2] - r[1] * a[3] + r[2] * a[0] + r[3] * a[1]);
	b[3] = (r[0] * a[3] + r[1] * a[2] - r[2] * a[1] + r[3] * a[0]);
}

static int rotate_quaternion_into_fundamental_zone(int num_generators, double (*generator)[4], double* q)
{
	double max = 0.0;
	int i = 0, bi = -1;
	for (i=0;i<num_generators;i++)
	{
		double* g = generator[i];
		double t = fabs(q[0] * g[0] - q[1] * g[1] - q[2] * g[2] - q[3] * g[3]);
		if (t > max)
		{
			max = t;
			bi = i;
		}
	}

	double f[4];
	quat_rot(q, generator[bi], f);
	memcpy(q, &f, 4 * sizeof(double));
	if (q[0] < 0)
	{
		q[0] = -q[0];
		q[1] = -q[1];
		q[2] = -q[2];
		q[3] = -q[3];
	}

	return bi;
}

int rotate_quaternion_into_cubic_fundamental_zone(double* q)
{
	return rotate_quaternion_into_fundamental_zone(24, generator_cubic, q);
}

int rotate_quaternion_into_diamond_cubic_fundamental_zone(double* q)
{
	return rotate_quaternion_into_fundamental_zone(12, generator_diamond_cubic, q);
}

int rotate_quaternion_into_icosahedral_fundamental_zone(double* q)
{
	return rotate_quaternion_into_fundamental_zone(60, generator_icosahedral, q);
}

int rotate_quaternion_into_hcp_fundamental_zone(double* q)
{
	return rotate_quaternion_into_fundamental_zone(6, generator_hcp, q);
}

int rotate_quaternion_into_hcp_crystalline_fundamental_zone(double* q)
{
	return rotate_quaternion_into_fundamental_zone(12, generator_hcp_crystalline, q);
}

int rotate_quaternion_into_diamond_hexagonal_fundamental_zone(double* q)
{
	return rotate_quaternion_into_fundamental_zone(3, generator_diamond_hexagonal, q);
}

void quaternion_to_rotation_matrix(double* q, double* u)
{
	double a = q[0];
	double b = q[1];
	double c = q[2];
	double d = q[3];

	u[0] = a*a + b*b - c*c - d*d;
	u[1] = 2*b*c - 2*a*d;
	u[2] = 2*b*d + 2*a*c;

	u[3] = 2*b*c + 2*a*d;
	u[4] = a*a - b*b + c*c - d*d;
	u[5] = 2*c*d - 2*a*b;

	u[6] = 2*b*d - 2*a*c;
	u[7] = 2*c*d + 2*a*b;
	u[8] = a*a - b*b - c*c + d*d;
}

}

