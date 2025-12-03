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
 * Here the eigenvalues of the MLS fluid part are computed.
 *
 * @param imove Moving flags.
 *   - imove > 0 for regular fluid particles.
 *   - imove = 0 for sensors.
 *   - imove < 0 for boundary elements/particles.
 * @param r Position \f$ \mathbf{r} \f$.
 * @param mls_fluid Kernel MLS fluid transformation matrix \f$ L \f$.
 * @param lambda Eigenvalues of the MLS fluid transformation matrix.
 * @param N Number of particles.
 */
__kernel void lambda(const __global int* imove,
                   const __global matrix* mls_fluid,
                   __global float* lambda,
                   usize N)
{
    const usize i = get_global_id(0);
    if(i >= N)
        return;
	if (imove[i] != 1)
		return;
    #ifndef HAVE_3D

	float a = mls_fluid[i].s0;
	float b = mls_fluid[i].s1;
	float c = mls_fluid[i].s2;
	float d = mls_fluid[i].s3;

	const float lambda1 = 1.f / ( ( 0.5f * (a + d) + 0.5f * sqrt(fmax(0.f, 4.f * b * c + (a - d) * (a - d) ) ) ) );
	const float lambda2 = 1.f / ( ( 0.5f * (a + d) - 0.5f * sqrt(fmax(0.f, 4.f * b * c + (a - d) * (a - d) ) ) ) );
	
    lambda[i] = fmin(lambda1,lambda2);

	#else

	const float c1 = MATRIX_TRACE(mls_fluid[i]) * MATRIX_TRACE(mls_fluid[i]);
	const float16 c2 = MATRIX_MUL(mls_fluid[i], mls_fluid[i]);

	const float a = - MATRIX_TRACE(mls_fluid[i]) ;
	const float b =   0.5f * (c1 - MATRIX_TRACE(c2));
	const float c = - det(mls_fluid[i]);

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

	lambda[i] = fmin(inter, 1.f/x3);
		}
		else {
 			const float x1 = - (3.f*q)/(2.f*p) - a/3.f;
			const float x2 = - (3.f*q)/(2.f*p) - a/3.f;
			const float x3 = - (4.f*p*p)/(9.f*q) - a/3.f;
	
	const float inter = min(1.f/x1,1.f/x2);

	lambda[i] = fmin(inter, 1.f/x3); 
		}
	}
	else if (delta > 0.00f){
		const float xi = -q/2.f - sqrt(delta);
		const float yi = -q/2.f + sqrt(delta);

		const float x1 = copysign(pow( fabs(yi) , 1.f/3.f ), yi) + copysign(pow( fabs(xi) , 1.f/3.f ), xi) - a/3.f;

	lambda[i] = 1.f/x1;
	}
	else {
		const float x1 = 2.f * sqrt(-p/3.f) * cos((phi + 2.f*0*pi)/3.f) - a/3.f;
		const float x2 = 2.f * sqrt(-p/3.f) * cos((phi + 2.f*1*pi)/3.f) - a/3.f;
		const float x3 = 2.f * sqrt(-p/3.f) * cos((phi + 2.f*2*pi)/3.f) - a/3.f;
	
	const float inter = min(1.f/x1,1.f/x2);

	lambda[i] = fmin(inter, 1.f/x3);

	}

    #endif
}

/** @brief Free surface normals computation.
 *
 * Here the normals of the fluid particles are computed.
 *
 * @param imove Moving flags.
 *   - imove > 0 for regular fluid particles.
 *   - imove = 0 for sensors.
 *   - imove < 0 for boundary elements/particles.
 * @param r Position \f$ \mathbf{r} \f$.
 * @param rho Density \f$ \rho \f$.
 * @param m Mass \f$ m \f$.
 * @param mls_fluid Kernel MLS fluid transformation matrix \f$ L \f$.
 * @param lambda Eigenvalues of the MLS fluid transformation matrix.
 * @param inormal Normal of the fluid particles.
 * @param jhoc Head and tail of chains for each cell.
 * @param N Number of particles.
 */

__kernel void normal(const __global int* imove,
                    const __global matrix* mls_fluid,
                    const __global float* lambda,
		    const __global float* m,
		    const __global float* rho,
		    const __global vec* r,
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
    if(imove[i] != 1)
	return;

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
        if(imove[j] != 1){
            continue;
        }
        const vec_xyz r_ij = r[j].XYZ - r_i;
        const float q = length(r_ij) / H;
        if(q >= SUPPORT)
        {
            continue;
        }
		const vec_xyz nabla_wij = kernelF(q) * CONF * r_ij;
		if(lambda[i] > 0.75f){
			_INORMAL_ -= (lambda[j] - lambda[i]) * MATRIX_DOT(mls_fluid[i], nabla_wij).XYZ * m[j] / rho[j] ;
		}
		else{
			_INORMAL_ -= (lambda[j]) * MATRIX_DOT(mls_fluid[i], nabla_wij).XYZ  * m[j] / rho[j];
		}
	}END_FOR_NEIGHS()

	if(length(_INORMAL_) > 0.1f * lambda[i]/H){
	    _INORMAL_ = fast_normalize(_INORMAL_);
	}
	else{
	    _INORMAL_ = VEC_ZERO.XYZ;
	}
	#ifdef LOCAL_MEM_SIZE
	    inormal[i].XYZ = _INORMAL_;
	#endif

}