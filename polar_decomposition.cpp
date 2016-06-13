/*******************************************************************************
 *  -/_|:|_|_\- 
 *
 *  This code is a modification of Theobald's QCP rotation code.
 *  It has been adapted to calculate the polar decomposition of a 3x3 matrix
 *  Adaption by PM Larsen
 *
 *  Original Author(s):	  Douglas L. Theobald
 *				  Department of Biochemistry
 *				  MS 009
 *				  Brandeis University
 *				  415 South St
 *				  Waltham, MA  02453
 *				  USA
 *
 *				  dtheobald@brandeis.edu
 *				  
 *				  Pu Liu
 *				  Johnson & Johnson Pharmaceutical Research and Development, L.L.C.
 *				  665 Stockton Drive
 *				  Exton, PA  19341
 *				  USA
 *
 *				  pliu24@its.jnj.com
 * 
 *
 *	If you use this QCP rotation calculation method in a publication, please
 *	reference:
 *
 *	  Douglas L. Theobald (2005)
 *	  "Rapid calculation of RMSD using a quaternion-based characteristic
 *	  polynomial."
 *	  Acta Crystallographica A 61(4):478-480.
 *
 *	  Pu Liu, Dmitris K. Agrafiotis, and Douglas L. Theobald (2009)
 *	  "Fast determination of the optimal rotational matrix for macromolecular 
 *	  superpositions."
 *	  Journal of Computational Chemistry 31(7):1561-1563.
 *
 *
 *  Copyright (c) 2009-2013 Pu Liu and Douglas L. Theobald
 *  All rights reserved.
 *
 *  Redistribution and use in source and binary forms, with or without modification, are permitted
 *  provided that the following conditions are met:
 *
 *  * Redistributions of source code must retain the above copyright notice, this list of
 *	conditions and the following disclaimer.
 *  * Redistributions in binary form must reproduce the above copyright notice, this list
 *	of conditions and the following disclaimer in the documentation and/or other materials
 *	provided with the distribution.
 *  * Neither the name of the <ORGANIZATION> nor the names of its contributors may be used to
 *	endorse or promote products derived from this software without specific prior written
 *	permission.
 *
 *  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 *  "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 *  LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
 *  PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 *  HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 *  SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 *  LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 *  DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 *  THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 *  (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 *  OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE. 
 *
 *  Source:		 started anew.
 *
 *  Change History:
 *	2009/04/13	  Started source
 *	2010/03/28	  Modified FastCalcRMSDAndRotation() to handle tiny qsqr
 *					If trying all rows of the adjoint still gives too small
 *					qsqr, then just return identity matrix. (DLT)
 *	2010/06/30	  Fixed prob in assigning A[9] = 0 in InnerProduct()
 *					invalid mem access
 *	2011/02/21	  Made CenterCoords use weights
 *	2011/05/02	  Finally changed CenterCoords declaration in qcprot.h
 *					Also changed some functions to static
 *	2011/07/08	  put in fabs() to fix taking sqrt of small neg numbers, fp error
 *	2012/07/26	  minor changes to comments and main.c, more info (v.1.4)
 *
 *      2016/05/29        QCP method adapted for polar decomposition of a 3x3 matrix.  For use in Polyhedral Template Matching.
 *  
 ******************************************************************************/

#include <cstdbool>
#include <cmath>


static void matmul(double* A, double* x, double* b)
{
	b[0] = A[0] * x[0] + A[1] * x[3] + A[2] * x[6];
	b[3] = A[3] * x[0] + A[4] * x[3] + A[5] * x[6];
	b[6] = A[6] * x[0] + A[7] * x[3] + A[8] * x[6];

	b[1] = A[0] * x[1] + A[1] * x[4] + A[2] * x[7];
	b[4] = A[3] * x[1] + A[4] * x[4] + A[5] * x[7];
	b[7] = A[6] * x[1] + A[7] * x[4] + A[8] * x[7];

	b[2] = A[0] * x[2] + A[1] * x[5] + A[2] * x[8];
	b[5] = A[3] * x[2] + A[4] * x[5] + A[5] * x[8];
	b[8] = A[6] * x[2] + A[7] * x[5] + A[8] * x[8];
}

static void quaternion_to_rotation_matrix(double* q, double* u)
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

int polar_decomposition_3x3(double* _A, bool right_sided, double* U, double* P)
{
	const double evecprec = 1e-6;
	const double evalprec = 1e-11;

	double A[9] = {_A[0], _A[1], _A[2], _A[3], _A[4], _A[5], _A[6], _A[7], _A[8]};
	double det = A[0] * (A[4]*A[8] - A[5]*A[7]) - A[1] * (A[3]*A[8] - A[5]*A[6]) + A[2] * (A[3]*A[7] - A[4]*A[6]);
	if (det < 0)
	{
		for (int i=0;i<9;i++)
			A[i] = -A[i];
	}

	double	Sxx = A[0], Sxy = A[1], Sxz = A[2],
		Syx = A[3], Syy = A[4], Syz = A[5],
		Szx = A[6], Szy = A[7], Szz = A[8];

	double	Sxx2 = Sxx * Sxx, Syy2 = Syy * Syy, Szz2 = Szz * Szz,
		Sxy2 = Sxy * Sxy, Syz2 = Syz * Syz, Sxz2 = Sxz * Sxz,
		Syx2 = Syx * Syx, Szy2 = Szy * Szy, Szx2 = Szx * Szx;

	double SyzSzymSyySzz2 = 2.0*(Syz*Szy - Syy*Szz);
	double Sxx2Syy2Szz2Syz2Szy2 = Syy2 + Szz2 - Sxx2 + Syz2 + Szy2;
	double SxzpSzx = Sxz + Szx;
	double SyzpSzy = Syz + Szy;
	double SxypSyx = Sxy + Syx;
	double SyzmSzy = Syz - Szy;
	double SxzmSzx = Sxz - Szx;
	double SxymSyx = Sxy - Syx;
	double SxxpSyy = Sxx + Syy;
	double SxxmSyy = Sxx - Syy;
	double Sxy2Sxz2Syx2Szx2 = Sxy2 + Sxz2 - Syx2 - Szx2;

	double C[3];
	C[0] = Sxy2Sxz2Syx2Szx2 * Sxy2Sxz2Syx2Szx2
		 + (Sxx2Syy2Szz2Syz2Szy2 + SyzSzymSyySzz2) * (Sxx2Syy2Szz2Syz2Szy2 - SyzSzymSyySzz2)
		 + (-(SxzpSzx)*(SyzmSzy)+(SxymSyx)*(SxxmSyy-Szz)) * (-(SxzmSzx)*(SyzpSzy)+(SxymSyx)*(SxxmSyy+Szz))
		 + (-(SxzpSzx)*(SyzpSzy)-(SxypSyx)*(SxxpSyy-Szz)) * (-(SxzmSzx)*(SyzmSzy)-(SxypSyx)*(SxxpSyy+Szz))
		 + (+(SxypSyx)*(SyzpSzy)+(SxzpSzx)*(SxxmSyy+Szz)) * (-(SxymSyx)*(SyzmSzy)+(SxzpSzx)*(SxxpSyy+Szz))
		 + (+(SxypSyx)*(SyzmSzy)+(SxzmSzx)*(SxxmSyy-Szz)) * (-(SxymSyx)*(SyzpSzy)+(SxzmSzx)*(SxxpSyy-Szz));

	C[1] = 8.0 * (Sxx*Syz*Szy + Syy*Szx*Sxz + Szz*Sxy*Syx - Sxx*Syy*Szz - Syz*Szx*Sxy - Szy*Syx*Sxz);

	C[2] = -2.0 * (Sxx2 + Syy2 + Szz2 + Sxy2 + Syx2 + Sxz2 + Szx2 + Syz2 + Szy2);

	double fnorm_squared = 0.0;
	for (int i=0;i<9;i++)
		fnorm_squared += A[i]*A[i];

	//Newton-Raphson
	double mxEigenV = sqrt(3 * fnorm_squared);
	if (mxEigenV > evalprec)
	{
		for (int i=0;i<50;i++)
		{
			double oldg = mxEigenV;
			double x2 = mxEigenV*mxEigenV;
			double b = (x2 + C[2])*mxEigenV;
			double a = b + C[1];
			double delta = ((a * mxEigenV + C[0]) / (2 * x2 * mxEigenV + b + a));
			mxEigenV -= delta;
			if (fabs(mxEigenV - oldg) < fabs(evalprec * mxEigenV))
				break;
		}
	}
	else
	{
		mxEigenV = 0.0;
	}

	double a11 = SxxpSyy + Szz - mxEigenV;
	double a12 = SyzmSzy;
	double a13 = -SxzmSzx;
	double a14 = SxymSyx;

	double a21 = SyzmSzy;
	double a22 = SxxmSyy - Szz  -mxEigenV;
	double a23 = SxypSyx;
	double a24 = SxzpSzx;

	double a31 = a13;
	double a32 = a23;
	double a33 = Syy - Sxx - Szz - mxEigenV;
	double a34 = SyzpSzy;

	double a41 = a14;
	double a42 = a24;
	double a43 = a34;
	double a44 = Szz - SxxpSyy - mxEigenV;

	double a3344_4334 = a33 * a44 - a43 * a34;
	double a3244_4234 = a32 * a44 - a42 * a34;
	double a3243_4233 = a32 * a43 - a42 * a33;
	double a3143_4133 = a31 * a43 - a41 * a33;
	double a3144_4134 = a31 * a44 - a41 * a34;
	double a3142_4132 = a31 * a42 - a41 * a32;

	double q1 =  a22*a3344_4334-a23*a3244_4234+a24*a3243_4233;
	double q2 = -a21*a3344_4334+a23*a3144_4134-a24*a3143_4133;
	double q3 =  a21*a3244_4234-a22*a3144_4134+a24*a3142_4132;
	double q4 = -a21*a3243_4233+a22*a3143_4133-a23*a3142_4132;

	double qsqr = q1 * q1 + q2 * q2 + q3 * q3 + q4 * q4;
	double q[4];

	bool too_small = false;
	//The following code tries to calculate another column in the adjoint matrix when the norm of the 
	//current column is too small.
	//Usually this block will never be activated.  To be absolutely safe this should be
	//uncommented, but it is most likely unnecessary.
	if (qsqr < evecprec)
	{
		q1 =  a12*a3344_4334 - a13*a3244_4234 + a14*a3243_4233;
		q2 = -a11*a3344_4334 + a13*a3144_4134 - a14*a3143_4133;
		q3 =  a11*a3244_4234 - a12*a3144_4134 + a14*a3142_4132;
		q4 = -a11*a3243_4233 + a12*a3143_4133 - a13*a3142_4132;
		qsqr = q1*q1 + q2*q2 + q3*q3 + q4*q4;

		if (qsqr < evecprec)
		{
			double a1324_1423 = a13 * a24 - a14 * a23, a1224_1422 = a12 * a24 - a14 * a22;
			double a1223_1322 = a12 * a23 - a13 * a22, a1124_1421 = a11 * a24 - a14 * a21;
			double a1123_1321 = a11 * a23 - a13 * a21, a1122_1221 = a11 * a22 - a12 * a21;

			q1 =  a42 * a1324_1423 - a43 * a1224_1422 + a44 * a1223_1322;
			q2 = -a41 * a1324_1423 + a43 * a1124_1421 - a44 * a1123_1321;
			q3 =  a41 * a1224_1422 - a42 * a1124_1421 + a44 * a1122_1221;
			q4 = -a41 * a1223_1322 + a42 * a1123_1321 - a43 * a1122_1221;
			qsqr = q1*q1 + q2*q2 + q3*q3 + q4*q4;

			if (qsqr < evecprec)
			{
				q1 =  a32 * a1324_1423 - a33 * a1224_1422 + a34 * a1223_1322;
				q2 = -a31 * a1324_1423 + a33 * a1124_1421 - a34 * a1123_1321;
				q3 =  a31 * a1224_1422 - a32 * a1124_1421 + a34 * a1122_1221;
				q4 = -a31 * a1223_1322 + a32 * a1123_1321 - a33 * a1122_1221;
				qsqr = q1*q1 + q2*q2 + q3*q3 + q4*q4;
				
				if (qsqr < evecprec)
				{
					//if qsqr is still too small, return the identity matrix.
					q[0] = 1.0;
					q[1] = 0.0;
					q[2] = 0.0;
					q[3] = 0.0;
					U[0] = U[4] = U[8] = 1.0;
					U[1] = U[2] = U[3] = U[5] = U[6] = U[7] = 0.0;
					too_small = true;
				}
			}
		}
	}

	if (!too_small)
	{
		double normq = sqrt(qsqr);
		q1 /= normq;
		q2 /= normq;
		q3 /= normq;
		q4 /= normq;
		q[0] = -q1;
		q[1] = q2;
		q[2] = q3;
		q[3] = q4;
		quaternion_to_rotation_matrix(q, U);
	}

	if (det < 0)
	{
		for (int i=0;i<9;i++)
			U[i] = -U[i];
	}

	double invU[9] = {U[0], U[3], U[6], U[1], U[4], U[7], U[2], U[5], U[8]};

	if (right_sided)
		matmul(invU, _A, P);
	else
		matmul(_A, invU, P);

	return !too_small;
}

