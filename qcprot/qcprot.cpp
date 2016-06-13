/*******************************************************************************
 *  -/_|:|_|_\- 
 *
 *  File:		   qcprot.cpp
 *  Version:		1.4
 *
 *  Function:	   Rapid calculation of the least-squares rotation using a 
 *				  quaternion-based characteristic polynomial and 
 *				  a cofactor matrix
 *
 *  Author(s):	  Douglas L. Theobald
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
 *  Modified by PM Larsen for use in Polyhedral Template Matching
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
 *      2016/02/20        modified for use in Polyhedral Template Matching.  InnerProduct function now takes permutation array.
 *  
 ******************************************************************************/
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include "qcprot.hpp"
#include "quat.hpp"

int FastCalcRMSDAndRotation(double *q, double *A, double *rmsd, double E0, int len, double minScore, double* rot)
{
	double 	SyzSzymSyySzz2, Sxx2Syy2Szz2Syz2Szy2, Sxy2Sxz2Syx2Szx2,
		   SxzpSzx, SyzpSzy, SxypSyx, SyzmSzy,
		   SxzmSzx, SxymSyx, SxxpSyy, SxxmSyy;
	double C[4];
	double mxEigenV; 
	double b, a, delta, rms, qsqr;
	double q1, q2, q3, q4, normq;
	double a11, a12, a13, a14, a21, a22, a23, a24;
	double a31, a32, a33, a34, a41, a42, a43, a44;
	double x2;
	double a3344_4334, a3244_4234, a3243_4233, a3143_4133,a3144_4134, a3142_4132; 
	double evecprec = 1e-6;
	double evalprec = 1e-11;

	double	Sxx = A[0], Sxy = A[1], Sxz = A[2],
		Syx = A[3], Syy = A[4], Syz = A[5],
		Szx = A[6], Szy = A[7], Szz = A[8];

	double	Sxx2 = Sxx * Sxx, Syy2 = Syy * Syy, Szz2 = Szz * Szz,
		Sxy2 = Sxy * Sxy, Syz2 = Syz * Syz, Sxz2 = Sxz * Sxz,
		Syx2 = Syx * Syx, Szy2 = Szy * Szy, Szx2 = Szx * Szx;

	SyzSzymSyySzz2 = 2.0*(Syz*Szy - Syy*Szz);
	Sxx2Syy2Szz2Syz2Szy2 = Syy2 + Szz2 - Sxx2 + Syz2 + Szy2;

	C[2] = -2.0 * (Sxx2 + Syy2 + Szz2 + Sxy2 + Syx2 + Sxz2 + Szx2 + Syz2 + Szy2);
	C[1] = 8.0 * (Sxx*Syz*Szy + Syy*Szx*Sxz + Szz*Sxy*Syx - Sxx*Syy*Szz - Syz*Szx*Sxy - Szy*Syx*Sxz);

	SxzpSzx = Sxz + Szx;
	SyzpSzy = Syz + Szy;
	SxypSyx = Sxy + Syx;
	SyzmSzy = Syz - Szy;
	SxzmSzx = Sxz - Szx;
	SxymSyx = Sxy - Syx;
	SxxpSyy = Sxx + Syy;
	SxxmSyy = Sxx - Syy;
	Sxy2Sxz2Syx2Szx2 = Sxy2 + Sxz2 - Syx2 - Szx2;

	C[0] = Sxy2Sxz2Syx2Szx2 * Sxy2Sxz2Syx2Szx2
		 + (Sxx2Syy2Szz2Syz2Szy2 + SyzSzymSyySzz2) * (Sxx2Syy2Szz2Syz2Szy2 - SyzSzymSyySzz2)
		 + (-(SxzpSzx)*(SyzmSzy)+(SxymSyx)*(SxxmSyy-Szz)) * (-(SxzmSzx)*(SyzpSzy)+(SxymSyx)*(SxxmSyy+Szz))
		 + (-(SxzpSzx)*(SyzpSzy)-(SxypSyx)*(SxxpSyy-Szz)) * (-(SxzmSzx)*(SyzmSzy)-(SxypSyx)*(SxxpSyy+Szz))
		 + (+(SxypSyx)*(SyzpSzy)+(SxzpSzx)*(SxxmSyy+Szz)) * (-(SxymSyx)*(SyzmSzy)+(SxzpSzx)*(SxxpSyy+Szz))
		 + (+(SxypSyx)*(SyzmSzy)+(SxzmSzx)*(SxxmSyy-Szz)) * (-(SxymSyx)*(SyzpSzy)+(SxzmSzx)*(SxxpSyy-Szz));

	//Newton-Raphson
	mxEigenV = E0;
	int i = 0;
	for (i = 0; i < 50; ++i)
	{
		double oldg = mxEigenV;
		x2 = mxEigenV*mxEigenV;
		b = (x2 + C[2])*mxEigenV;
		a = b + C[1];
		delta = ((a*mxEigenV + C[0])/(2.0*x2*mxEigenV + b + a));
		mxEigenV -= delta;
		//printf("\n diff[%3d]: %16g %16g %16g", i, mxEigenV - oldg, evalprec*mxEigenV, mxEigenV);
		if (fabs(mxEigenV - oldg) < fabs(evalprec*mxEigenV))
			break;
	}

	//if (i == 50) 
	//   fprintf(stderr,"\nMore than %d iterations needed!\n", i);

	//the fabs() is to guard against extremely small, but *negative* numbers due to floating point error
	rms = sqrt(fabs(2.0 * (E0 - mxEigenV)/len));
	(*rmsd) = rms;
	//printf("\n\n %16g %16g %16g \n", rms, E0, 2.0 * (E0 - mxEigenV)/len);

	if (minScore > 0) 
		if (rms < minScore)
			return -1; // Don't bother with rotation. 

	a11 = SxxpSyy + Szz - mxEigenV;
	a12 = SyzmSzy;
	a13 = -SxzmSzx;
	a14 = SxymSyx;

	a21 = SyzmSzy;
	a22 = SxxmSyy - Szz  -mxEigenV;
	a23 = SxypSyx;
	a24 = SxzpSzx;

	a31 = a13;
	a32 = a23;
	a33 = Syy - Sxx - Szz - mxEigenV;
	a34 = SyzpSzy;

	a41 = a14;
	a42 = a24;
	a43 = a34;
	a44 = Szz - SxxpSyy - mxEigenV;

	a3344_4334 = a33 * a44 - a43 * a34;
	a3244_4234 = a32 * a44-a42*a34;
	a3243_4233 = a32 * a43 - a42 * a33;
	a3143_4133 = a31 * a43-a41*a33;
	a3144_4134 = a31 * a44 - a41 * a34;
	a3142_4132 = a31 * a42-a41*a32;

	q1 =  a22*a3344_4334-a23*a3244_4234+a24*a3243_4233;
	q2 = -a21*a3344_4334+a23*a3144_4134-a24*a3143_4133;
	q3 =  a21*a3244_4234-a22*a3144_4134+a24*a3142_4132;
	q4 = -a21*a3243_4233+a22*a3143_4133-a23*a3142_4132;

	qsqr = q1 * q1 + q2 * q2 + q3 * q3 + q4 * q4;

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
		qsqr = q1*q1 + q2 *q2 + q3*q3+q4*q4;

		if (qsqr < evecprec)
		{
			double a1324_1423 = a13 * a24 - a14 * a23, a1224_1422 = a12 * a24 - a14 * a22;
			double a1223_1322 = a12 * a23 - a13 * a22, a1124_1421 = a11 * a24 - a14 * a21;
			double a1123_1321 = a11 * a23 - a13 * a21, a1122_1221 = a11 * a22 - a12 * a21;

			q1 =  a42 * a1324_1423 - a43 * a1224_1422 + a44 * a1223_1322;
			q2 = -a41 * a1324_1423 + a43 * a1124_1421 - a44 * a1123_1321;
			q3 =  a41 * a1224_1422 - a42 * a1124_1421 + a44 * a1122_1221;
			q4 = -a41 * a1223_1322 + a42 * a1123_1321 - a43 * a1122_1221;
			qsqr = q1*q1 + q2 *q2 + q3*q3+q4*q4;

			if (qsqr < evecprec)
			{
				q1 =  a32 * a1324_1423 - a33 * a1224_1422 + a34 * a1223_1322;
				q2 = -a31 * a1324_1423 + a33 * a1124_1421 - a34 * a1123_1321;
				q3 =  a31 * a1224_1422 - a32 * a1124_1421 + a34 * a1122_1221;
				q4 = -a31 * a1223_1322 + a32 * a1123_1321 - a33 * a1122_1221;
				qsqr = q1*q1 + q2 *q2 + q3*q3 + q4*q4;
				
				if (qsqr < evecprec)
				{
					//if qsqr is still too small, return the identity matrix.
					q[0] = 1.0;
					q[1] = 0.0;
					q[2] = 0.0;
					q[3] = 0.0;
					//rot[0] = rot[4] = rot[8] = 1.0;
					//rot[1] = rot[2] = rot[3] = rot[5] = rot[6] = rot[7] = 0.0;
					return 0;
				}
			}
		}
	}

	normq = sqrt(qsqr);
	q1 /= normq;
	q2 /= normq;
	q3 /= normq;
	q4 /= normq;
	q[0] = q1;
	q[1] = q2;
	q[2] = q3;
	q[3] = q4;

	quaternion_to_rotation_matrix(q, rot);
	return 1;
}

void InnerProduct(double *A, int num, const double (*coords1)[3], double (*coords2)[3], int8_t* permutation)
{
	A[0] = A[1] = A[2] = A[3] = A[4] = A[5] = A[6] = A[7] = A[8] = 0.0;

	for (int i = 0; i < num; ++i)
	{
		double x1 = coords1[i][0];
		double y1 = coords1[i][1];
		double z1 = coords1[i][2];

		double x2 = coords2[permutation[i]][0];
		double y2 = coords2[permutation[i]][1];
		double z2 = coords2[permutation[i]][2];

		A[0] += x1 * x2;
		A[1] += x1 * y2;
		A[2] += x1 * z2;

		A[3] += y1 * x2;
		A[4] += y1 * y2;
		A[5] += y1 * z2;

		A[6] += z1 * x2;
		A[7] += z1 * y2;
		A[8] += z1 * z2;  
	}
}

