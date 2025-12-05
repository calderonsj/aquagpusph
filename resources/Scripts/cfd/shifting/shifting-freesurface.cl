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

#if defined(LOCAL_MEM_SIZE) && defined(NO_LOCAL_MEM)
    #error NO_LOCAL_MEM has been set.
#endif

#include "resources/Scripts/types/types.h"
#include "resources/Scripts/KernelFunctions/Kernel.h"

/** @brief Kappa factor computation.
 *
 * @param imove Moving flags.
 *   - imove > 0 for regular fluid particles.
 *   - imove = 0 for sensors.
 *   - imove < 0 for boundary elements/particles.
 * @param r Position \f$ \mathbf{r} \f$.
 * @param kappa Curvature correction factor at free surface.
 * @param inormal Normal of the fluid particles.
 * @param jhoc Head of chain for each cell (first particle found).
 * @param N Number of particles.
 */

__kernel void threshold(const __global int* imove,
                    const __global vec* r,
		    __global int* kappa,
		    const __global vec* inormal,
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

	const vec_xyz r_i = r[i].XYZ;

    // Initialize the output
    #ifndef LOCAL_MEM_SIZE
        #define _KAPPA_ kappa[i]
    #else
        #define _KAPPA_ kappa_l[it]
        __local int kappa_l[LOCAL_MEM_SIZE];
    #endif

    _KAPPA_ = 1;

    const vec_xyz n_i = inormal[i].XYZ;
    if(length(n_i) < 1e-6f){
        return;
    }

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

        const vec_xyz n_j = inormal[j].XYZ;
        if(length(n_j) < 1e-6f){
            continue;
        }

	    float norma = dot(fast_normalize(n_i), fast_normalize(n_j));
	    norma = clamp(norma, -1.0f, 1.0f);

	if(acos(norma) > 0.2617f)
        {
            _KAPPA_ = 0;
        }
    }END_FOR_NEIGHS()

    #ifdef LOCAL_MEM_SIZE
	kappa[i] = _KAPPA_;
    #endif
}

/** @brief Definition of region around free surface.
 *
 * @param imove Moving flags.
 *   - imove > 0 for regular fluid particles.
 *   - imove = 0 for sensors.
 *   - imove < 0 for boundary elements/particles.
 * @param r Position \f$ \mathbf{r} \f$.
 * @param frees Free surface flag.
 *   - frees = 1 for particles belonging to the free surface.
 *   - frees = 0 for particles not belonging to the free surface.
 * @param region Area around free surface (up to 2h).
 * @param jhoc Head of chain for each cell (first particle found).
 * @param N Number of particles.
 */
__kernel void region(const __global int* imove,
                    const __global vec* r,
		    const __global unsigned int* frees,
		    __global unsigned int* region,
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
 
    const vec_xyz r_i = r[i].XYZ;

    // Initialize the output
    #ifndef LOCAL_MEM_SIZE
        #define _REGION_ region[i]
    #else
        #define _REGION_ region_l[it]
        __local unsigned int region_l[LOCAL_MEM_SIZE];
    #endif
        _REGION_ = 0;

    FOR_NEIGHS(N, jhoc){
        if((imove[j] != 1)){
            continue;
        }

        const vec_xyz r_ij = r[j].XYZ - r_i;
        const float q = length(r_ij) / H;

        if(q >= SUPPORT)
        {
            continue;
        }

        {
	    if( (frees[j] == 1) ){
		_REGION_ = 1;
		break;
	    }
        }
    }END_FOR_NEIGHS()

    #ifdef LOCAL_MEM_SIZE
        region[i] = _REGION_;
    #endif
}

/** @brief Shifting correction at the free surface.
 *
 * @param imove Moving flags.
 *   - imove > 0 for regular fluid particles.
 *   - imove = 0 for sensors.
 *   - imove < 0 for boundary elements/particles.
 * @param lambda_bi Eigenvalues of the MLS total transformation matrix.
 * @param kappa Curvature correction factor at free surface.
 * @param inormal Normal of the fluid particles.
 * @param region Area around free surface (up to 2h).
 * @param N Number of particles.
 */
__kernel void entry(const __global int* imove,
		    const __global float* lambda_bi,
		    const __global int* kappa,
		    const __global vec* inormal,
		    const __global int* frees, 
		    const __global unsigned int* region, 
		    __global vec* du_shift,
                    usize N)
{
    const usize i = get_global_id(0);
    if(i >= N)
        return;
    if(imove[i] != 1)
        return;

    if(region[i] == 0)
	return;

    const vec_xyz n_i = inormal[i].XYZ;    
        
    //2019 PN Sun based
    // Conditions are applied based on the free surface region, not FS itself!
    // If particle is outside the free surface region, keep shifting.
    if(region[i] == 0){
        du_shift[i] = du_shift[i] ;
    }
    // If particle is inside the free surface region, conditions from PN Sun.
    else{
        if(lambda_bi[i] < 0.55f){
            du_shift[i] = VEC_ZERO;
        }
        else{
            const matrix out = outer(inormal[i].XYZ, inormal[i].XYZ);
            const matrix I = MAT_EYE - out;
            const vec R = MATRIX_DOT( I, du_shift[i]);

            if (dot(inormal[i].XYZ, du_shift[i].XYZ) >= 0.f){
                du_shift[i] = kappa[i] * R;
            }
            else{
                du_shift[i] = du_shift[i] ;
            }
        }
    }       
}
