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

/** @brief Free surface detector step 1.
 *
 * Here the first step of the free surface detector algorithm is computed.
 *
 * @param imove Moving flags.
 *   - imove > 0 for regular fluid particles.
 *   - imove = 0 for sensors.
 *   - imove < 0 for boundary elements/particles.
 * @param frees Free surface flag.
 *   - frees = 1 for particles belonging to the free surface.
 *   - frees = 0 for particles not belonging to the free surface.
 * @param lambda Eigenvalues of the MLS fluid transformation matrix.
 * @param N Number of particles.
 */
__kernel void first(const __global int* imove,          
		    __global unsigned int* frees,
                    const __global float* lambda,
                    usize N)
{
    const usize i = get_global_id(0);
    if(i >= N)
        return;
	if (imove[i] != 1)
		return;

	if (lambda[i] < 0.2f){
	    frees[i] = 1;
	    return;
	}
}

/** @brief Free surface detector step 2.
 *
 * Here the second step of the free surface detector algorithm is computed.
 *
 * @param imove Moving flags.
 *   - imove > 0 for regular fluid particles.
 *   - imove = 0 for sensors.
 *   - imove < 0 for boundary elements/particles.
 * @param frees Free surface flag.
 *   - frees = 1 for particles belonging to the free surface.
 *   - frees = 0 for particles not belonging to the free surface.
 * @param lambda Eigenvalues of the MLS fluid transformation matrix.
 * @param inormal Normal of the fluid particles.
 * @param N Number of particles.
 */
__kernel void second(const __global int* imove,
                    const __global vec* r,
                    const __global vec* inormal,
		    const __global float* lambda,
		    __global unsigned int* frees,
                    // Link-list data
                    const __global svec2* jhoc,
                    // Simulation data
                    usize N)
{
    const usize i = get_global_id(0);
    const usize it = get_local_id(0);
    if(i >= N)
        return;
    if(imove[i] != 1)
        return;
    if(frees[i] == 1)
        return;

    const vec_xyz r_i = r[i].XYZ;
    const vec_xyz n_i = inormal[i].XYZ;

    const vec_xyz dt_i = H * n_i;
    const vec_xyz t_i = r_i + dt_i;

    #ifndef HAVE_3D
    vec_xyz tau_i = VEC_ZERO;
    	tau_i.x = -n_i.y;
	tau_i.y = n_i.x;
    #else
    vec_xyz tau_i = cross(r_i, n_i);
    #endif
 
    // Initialize the output
    #ifndef LOCAL_MEM_SIZE
        #define _FREES_ frees[i]
    #else
        #define _FREES_ frees_l[it]
        __local unsigned int frees_l[LOCAL_MEM_SIZE];
    #endif

    _FREES_ = 1;

    if(lambda[i] > 0.75f){
	frees[i] = 0;
	return;
    }

    FOR_NEIGHS(N, jhoc){

        const vec_xyz r_ij = r[j].XYZ - r_i;
	const vec_xyz r_jT = r[j].XYZ - t_i;
        const float q = length(r_ij) / H;

	if(q >= SUPPORT)
        {
            continue;
        }
	if(imove[j] != 1){
		continue;
	}
	#ifndef HAVE_3D
	{
		if((length(r_ij) >= sqrt(2.f) * H) && (length(r_jT) < H)){
			_FREES_ = 0;
		}
		if((length(r_ij) < sqrt(2.f) * H) && (fabs(dot(n_i, r_jT)) + fabs(dot(tau_i, r_jT)) < 0.99f * H) ){
			_FREES_ = 0;
		}
	}
	#else
	{
		if((length(r_ij) >= sqrt(2.f) * H) && (length(r_jT) < H)){ 
			_FREES_ = 0;
		}
		if((length(r_ij) < sqrt(2.f) * H) && (acos( dot(inormal[i].XYZ, r_ij) / length(r_ij) ) < 0.4f * 3.14159265359f )){
			_FREES_ = 0;
		}
	}
	#endif
    }END_FOR_NEIGHS()

    #ifdef LOCAL_MEM_SIZE
	frees[i] = _FREES_;
    #endif
}