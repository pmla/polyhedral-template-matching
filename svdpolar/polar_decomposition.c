//#####################################################################
// Copyright (c) 2009-2011, Eftychios Sifakis.
//
// Modified by PM Larsen for use in Polyhedral Template Matching
//
// Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
//   * Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
//   * Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or
//     other materials provided with the distribution.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING,
// BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT
// SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
// INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//#####################################################################

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define USE_SCALAR_IMPLEMENTATION
#define COMPUTE_V_AS_MATRIX
#define COMPUTE_U_AS_MATRIX

#define ENABLE_SCALAR_IMPLEMENTATION(X) X
#define ENABLE_SSE_IMPLEMENTATION(X)
#define ENABLE_AVX_IMPLEMENTATION(X)

static float rsqrt(const float f)
{
	return 1 / sqrt(f);
}

static void matmul(double* A, double* B, double* C)
{
	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			double acc = 0;
			for (int k = 0; k < 3; k++)
			{
				acc += A[i*3 + k] * B[k*3 + j];
			}

			C[i*3 + j] = acc;
		}
	}
}

static void transpose_matrix(double* V)
{
	double temp;

	temp = V[1];
	V[1] = V[3];
	V[3] = temp;

	temp = V[2];
	V[2] = V[6];
	V[6] = temp;

	temp = V[5];
	V[5] = V[7];
	V[7] = temp;
}

/*static void print_matrix(char c, double* V)
{
	printf("%c:\n", c);
	printf("%-8f %-8f %-8f\n", V[0], V[1], V[2]);
	printf("%-8f %-8f %-8f\n", V[3], V[4], V[5]);
	printf("%-8f %-8f %-8f\n", V[6], V[7], V[8]);
}*/

static void svd_3x3(double* A, double* U, double* S, double* V)
{
	#include "Singular_Value_Decomposition_Kernel_Declarations.h"

	ENABLE_SCALAR_IMPLEMENTATION(Sa11.f=A[0];)
	ENABLE_SCALAR_IMPLEMENTATION(Sa21.f=A[3];)
	ENABLE_SCALAR_IMPLEMENTATION(Sa31.f=A[6];)
	ENABLE_SCALAR_IMPLEMENTATION(Sa12.f=A[1];)
	ENABLE_SCALAR_IMPLEMENTATION(Sa22.f=A[4];)
	ENABLE_SCALAR_IMPLEMENTATION(Sa32.f=A[7];)
	ENABLE_SCALAR_IMPLEMENTATION(Sa13.f=A[2];)
	ENABLE_SCALAR_IMPLEMENTATION(Sa23.f=A[5];)
	ENABLE_SCALAR_IMPLEMENTATION(Sa33.f=A[8];)

	#include "Singular_Value_Decomposition_Main_Kernel_Body.h"
}

void left_sided_polar_decomposition_3x3(double* A, double* P, double* U)
{
	double S[9], W[9], Vt[9], t[9];

	svd_3x3(A, W, S, Vt);
	transpose_matrix(Vt);
	//W S Vt = A

	//Calculate U
	matmul(W, Vt, U);

	//Calculate P
	matmul(W, S, t);
	transpose_matrix(W);
	matmul(t, W, P);

	/*for (int i = 0;i<3;i++)
	{
		double s0 = U[i*3+0]*U[i*3+0] + U[i*3+1]*U[i*3+1] + U[i*3+2]*U[i*3+2];
		double s1 = U[i+0]*U[i+0] + U[i+3]*U[i+3] + U[i+6]*U[i+6];
		//printf("%f\t%f\n", s0, s1);

		//if (fabs(s0 - 1) > 1E-5) printf("s0: %f\n", s0);
		//if (fabs(s1 - 1) > 1E-5) printf("s1: %f\n", s1);

		//assert(fabs(s0 - 1) < 1E-6);
		//assert(fabs(s1 - 1) < 1E-6);
	}*/
}

