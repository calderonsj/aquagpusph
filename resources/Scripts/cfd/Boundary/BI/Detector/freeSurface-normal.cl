/*
 *  This file is part of AQUAgpusph, a free CFD program based on SPH.
 *  Copyright (C) 2012  Jose Luis Cercos Pita <jl.cercos@upm.es>
 *
 *  AQUAgpusph is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  AQUAgpusph is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with AQUAgpusph.  If not, see <http://www.gnu.org/licenses/>.
 */

/** @addtogroup basic
 * @{
 */

/** @file
 * @brief delta-SPH methods, including the correction terms
 */

#if defined(LOCAL_MEM_SIZE) && defined(NO_LOCAL_MEM)
    #error NO_LOCAL_MEM has been set.
#endif

#include "resources/Scripts/types/types.h"
#include "resources/Scripts/KernelFunctions/Kernel.h"

/** @brief Eigenvalues computation.
 *
 * Here the eigenvalues of the MLS transformation matrix, including BI terms, are computed.
 *
 * @param imove Moving flags.
 *   - imove > 0 for regular fluid particles.
 *   - imove = 0 for sensors.
 *   - imove < 0 for boundary elements/particles.
 * @param r Position \f$ \mathbf{r} \f$.
 * @param mls_bi Kernel MLS transformation matrix with BI \f$ L \f$.
 * @param lambda_bi Eigenvalues of the MLS transformation matrix with BI.
 * @param N Number of particles.
 */
__kernel void lambdabi(const __global int* imove,
                   const __global matrix* mls_bi,
                   __global float* lambda_bi,
                   usize N)
{
    const usize i = get_global_id(0);
    if(i >= N)
        return;
	if (imove[i] != 1)
		return;
    #ifndef HAVE_3D

	float a = mls_bi[i].s0;
	float b = mls_bi[i].s1;
	float c = mls_bi[i].s2;
	float d = mls_bi[i].s3;

	const float lambdabi1 = 1.f / ( 0.5f * (a + d) + 0.5f * sqrt(fmax(0.f, 4.f * b * c + (a - d) * (a - d) ) ) );
	const float lambdabi2 = 1.f / ( 0.5f * (a + d) - 0.5f * sqrt(fmax(0.f, 4.f * b * c + (a - d) * (a - d) ) ) );
	
    lambda_bi[i] = fmin(lambdabi1,lambdabi2);

	#else

	const float c1 = MATRIX_TRACE(mls_bi[i]) * MATRIX_TRACE(mls_bi[i]);
	const float16 c2 = MATRIX_MUL(mls_bi[i], mls_bi[i]);

	const float a = - MATRIX_TRACE(mls_bi[i]) ;
	const float b =   0.5f * (c1 - MATRIX_TRACE(c2));
	const float c = - det(mls_bi[i]);

	const float p = (3.f * b - a*a) / 3.f;
	const float q = (2.f*a*a*a-9.f*a*b+27.f*c)/27.f;
	const float delta = pow( q/2.f , 2.f) + pow( p/3.f , 3.f);
	const float phi = acos( (-q/2.f) / sqrt(-pow(p/3.f , 3.f) ) );
	const float pi = 3.14159f;

	if(delta == 0.00f){
		if((p == 0.f)  && (q == 0.f)){
			const float x1 = - a/3.f;
			const float x2 = - a/3.f;
			const float x3 = - a/3.f;
	
	const float inter = min(1.f/x1,1.f/x2);

	lambda_bi[i] = min(inter, 1.f/x3);
		}
		else {
 			const float x1 = - (3.f*q)/(2.f*p) - a/3.f;
			const float x2 = - (3.f*q)/(2.f*p) - a/3.f;
			const float x3 = - (4.f*p*p)/(9.f*q) - a/3.f;
	
	const float inter = min(1.f/x1,1.f/x2);

	lambda_bi[i] = min(inter, 1.f/x3); 
		}
	}
	else if (delta > 0.00f){
		const float xi = -q/2.f - sqrt(delta);
		const float yi = -q/2.f + sqrt(delta);

		const float x1 = copysign(pow( fabs(yi) , 1.f/3.f ), yi) + copysign(pow( fabs(xi) , 1.f/3.f ), xi) - a/3.f;

	lambda_bi[i] = 1.f/x1;
	}
	else {
		const float x1 = 2.f * sqrt(-p/3.f) * cos((phi + 2.f*0*pi)/3.f) - a/3.f;
		const float x2 = 2.f * sqrt(-p/3.f) * cos((phi + 2.f*1*pi)/3.f) - a/3.f;
		const float x3 = 2.f * sqrt(-p/3.f) * cos((phi + 2.f*2*pi)/3.f) - a/3.f;
	
	const float inter = min(1.f/x1,1.f/x2);

	lambda_bi[i] = min(inter, 1.f/x3);

	}

    #endif
}

/** @brief Free surface normals computation.
 *
 * Here the normals of the fluid particles are computed, taking into account
 * the influence of BI boundaries.
 *
 * @param imove Moving flags.
 *   - imove > 0 for regular fluid particles.
 *   - imove = 0 for sensors.
 *   - imove < 0 for boundary elements/particles.
 * @param r Position \f$ \mathbf{r} \f$.
 * @param rho Density \f$ \rho \f$.
 * @param m Mass \f$ m \f$.
 * @param mls_bi Kernel MLS transformation matrix with BI \f$ L \f$.
 * @param lambda_bi Eigenvalues of the MLS transformation matrix with BI.
 * @param inormal Normal of the fluid particles.
 * @param jhoc Head and tail of chains for each cell.
 * @param N Number of particles.
 */

__kernel void normal(const __global int* imove,
                    const __global matrix* mls_bi,
                    const __global float* lambda_bi,
                    const __global float* shepard,
		    const __global float* m,
		    const __global float* rho,
		    const __global vec* r,
		    const __global vec* normal,
		    __global vec* inormal,
                    // Link-list data
                    const __global svec2* jhoc,
                    // Simulation data
                    usize N)
{
    const usize i = get_global_id(0);
    const uint it = get_local_id(0);
    if(i >= N)
        return;
	if(imove[i] != 1){
	    inormal[i].XYZ = VEC_ZERO.XYZ;
	    return;
	}

	const vec_xyz r_i = r[i].XYZ;

    // Initialize the output
    #ifndef LOCAL_MEM_SIZE
        #define _INORMAL_ inormal[i].XYZ
    #else
        #define _INORMAL_ inormal_l[it]
        __local vec_xyz inormal_l[LOCAL_MEM_SIZE];
	__INORMAL_ = VEC_ZERO.XYZ;
    #endif

	FOR_NEIGHS(N, jhoc){
        if(i == j){
            continue;
        }
        if((imove[j] != 1) && (imove[j] != -3)){
            continue;
        }
        const vec_xyz r_ij = r[j].XYZ - r_i;
        const float q = length(r_ij) / H;
        if(q >= SUPPORT)
        {
            continue;
        }
		const vec_xyz n_j = normal[j].XYZ; // Assumed outwarding oriented
		const vec_xyz nabla_wij = kernelF(q) * CONF * r_ij;
		const vec_xyz nwj_Wij_Sj = n_j * kernelW(q) * CONW * m[j];
		if(lambda_bi[i] > 0.75f){
			if(imove[j] == 1){
				_INORMAL_ -= (lambda_bi[j] - lambda_bi[i]) * MATRIX_DOT(mls_bi[i], nabla_wij).XYZ * m[j] / rho[j] ;
			}
			if(imove[j] == -3){
				const vec_xyz n_j = normal[j].XYZ; // Assumed outwarding oriented
				_INORMAL_ -= (1.f - lambda_bi[i]) * MATRIX_DOT(mls_bi[i], nwj_Wij_Sj).XYZ; 
				//the 1.f would be lambda_bi[j] at a boundary term, but before changing all we just change this.
			}
		}
		else{
			if(imove[j] == 1){
				_INORMAL_ -= (lambda_bi[j]) * MATRIX_DOT(mls_bi[i], nabla_wij).XYZ  * m[j] / rho[j];
			}
			if(imove[j] == -3){
				const vec_xyz n_j = normal[j].XYZ; // Assumed outwarding oriented
				_INORMAL_ -= (1.f) * MATRIX_DOT(mls_bi[i], nwj_Wij_Sj).XYZ; 
				//the 1.f would be lambda_bi[j] at a boundary term, but before changing all we just change this.
			}
		}

	}END_FOR_NEIGHS()

	_INORMAL_ /= shepard[i];

	if(length(_INORMAL_) > 0.1f * lambda_bi[i]/H){
		_INORMAL_ = fast_normalize(_INORMAL_);
	}
	else{
		_INORMAL_ = VEC_ZERO.XYZ;
	}

	#ifdef LOCAL_MEM_SIZE
	    inormal[i].XYZ = _INORMAL_;
	#endif

}