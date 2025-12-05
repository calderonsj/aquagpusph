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

/** @brief Gradient of shifting formula computation.
 *
 * The gradient of the shifting formula is computed for
 * boundary integrals particle influence.
 *
 * @param imove Moving flags.
 *   - imove > 0 for regular fluid particles.
 *   - imove = 0 for sensors.
 *   - imove < 0 for boundary elements/particles.
 * @param r Position \f$ \mathbf{r} \f$.
 * @param rho Density \f$ \rho \f$.
 * @param m Mass \f$ m \f$.
 * @param normal Normal \f$ \mathbf{n} \f$.
 * @param shepard Shepard term
 * \f$ \gamma(\mathbf{x}) = \int_{\Omega}
 *     W(\mathbf{y} - \mathbf{x}) \mathrm{d}\mathbf{y} \f$.
 * @param gradC Gradient of the shifting formula.
 * @param dr Characteristic distance between particles.
 * @param jhoc Head of chain for each cell (first particle found).
 * @param N Number of particles.
 */
__kernel void drbound(const __global int* imove,
                    const __global vec* r,
                    const __global float* rho,
                    const __global float* m,
                    __global vec* gradC,
		    const __global vec* normal,
		    const __global float* shepard,
		    float dr,
                    // Link-list data
                    const __global svec2* jhoc,
                    // Simulation data
                    const usize N)
{
    const uint i = get_global_id(0);
    const uint it = get_local_id(0);
    if(i >= N)
        return;
    if(imove[i] != 1)
        return;

    const vec_xyz r_i = r[i].XYZ;

    // Initialize the output
    #ifndef LOCAL_MEM_SIZE
        #define _GRADC_ gradC[i]
    #else
        #define _GRADC_ gradC_l[it]
        __local vec gradC_l[LOCAL_MEM_SIZE];
        _GRADC_ = gradC[i];
    #endif

    FOR_NEIGHS(N, jhoc){
        if(imove[j] != -3){
            continue;
        }

        const vec_xyz r_ij = r[j].XYZ - r_i;
        const float q = length(r_ij) / H;

        if(q >= SUPPORT)
        {
            continue;
        }

        {
	    const float w_ij = kernelW(q) * CONW * m[j];
	    const vec_xyz n_j = normal[j].XYZ;
	    // No anti-clumping term for the boundary elements contribution.
	    _GRADC_ += w_ij * n_j;
	}
    }END_FOR_NEIGHS()

    const float shepard_i = shepard[i];

    #ifdef LOCAL_MEM_SIZE
	gradC[i] = _DUSHIFT_ / shepard_i;
    #endif
}