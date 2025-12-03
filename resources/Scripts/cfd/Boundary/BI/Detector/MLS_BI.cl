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
 * @brief MLS kernel transformation matrix computation
 */

#if defined(LOCAL_MEM_SIZE) && defined(NO_LOCAL_MEM)
    #error NO_LOCAL_MEM has been set.
#endif

#include "resources/Scripts/types/types.h"
#include "resources/Scripts/KernelFunctions/Kernel.h"

/** @brief Compute the MLS transformation matrix inverse, \f$ L_i^{-1} \f$.
 * 
 * Such transformation matrix can be multiplied by the kernel gradient to
 * produce a new kernel gradient,
 * \f$ \nabla W^{L}_{ij} = L_i \cdot \nabla W_{ij} \f$, such that the lienar
 * fields differential operators are consistently computed.
 *
 * @param imove Moving flags.
 *   - imove > 0 for regular fluid/solid particles.
 *   - imove = 0 for sensors.
 *   - imove < 0 for boundary elements/particles.
 * @param r Position \f$ \mathbf{r}_{n+1} \f$.
 * @param rho Density \f$ \rho \f$.
 * @param m Mass \f$ m \f$.
 * @param mls Kernel MLS transformation matrix \f$ L \f$.
 * @param mls_bi Kernel MLS transformation matrix \f$ L \f$, including BI.
 * @param jhoc Head of chain for each cell (first particle found).
 * @param N Number of particles.
 * @param n_cells Number of cells in each direction
 * @param mls_imove Type of particles affected
 * @note The MLS kernel transformation will be computed just for the particles
 * with the moving flag mls_imove, and using just the information of the
 * particles with the moving flag mls_imove
 */
__kernel void entry(const __global int* imove,
                    const __global vec* r,
                    const __global float* rho,
                    const __global float* m,
		    const __global float* shepard,
		    const __global vec* normal,
		    const __global matrix* mls,
                    __global matrix* mls_bi,
                    const __global svec2* jhoc,
                    usize N,
                    uint mls_imove)
{
    const usize i = get_global_id(0);
    const usize it = get_local_id(0);
    if(i >= N)
        return;
    if(imove[i] != mls_imove){
        return;
    }

    const vec_xyz r_i = r[i].XYZ;

    // Initialize the output
    #ifndef LOCAL_MEM_SIZE
        #define _MLS_BI_ mls_bi[i]
    #else
        #define _MLS_BI_ mls_bi_l[it]
        __local matrix mls_bi_l[LOCAL_MEM_SIZE];
    #endif
    _MLS_BI_ = mls[i];

    FOR_NEIGHS(N, jhoc){
        if((i == j) || (imove[j] != -3)){
            continue;
        }
        const vec_xyz r_ij = r[j].XYZ - r_i;
        const float q = length(r_ij) / H;
		const vec_xyz n_j = normal[j].XYZ;  // Assumed outwarding oriented
        if(q >= SUPPORT)
        {
            continue;
        }
        {
            _MLS_BI_ += outer(r_ij, kernelW(q) * CONW * n_j * m[j]);
        }
    }END_FOR_NEIGHS()

    #ifdef LOCAL_MEM_SIZE
        mls_bi[i] = _MLS_BI_;
    #endif
}

/** @brief Invert the matrix computed in entry() to get the final MLS
 * transformation matrix, \f$ L_i \f$.
 * 
 * Such transformation matrix can be multiplied by the kernel gradient to
 * produce a new kernel gradient,
 * \f$ \nabla W^{L}_{ij} = L_i \cdot \nabla W_{ij} \f$, such that the lienar
 * fields differential operators are consistently computed.
 *
 * @param imove Moving flags.
 *   - imove > 0 for regular fluid/solid particles.
 *   - imove = 0 for sensors.
 *   - imove < 0 for boundary elements/particles.
 * @param mls_bi Kernel MLS transformation matrix \f$ L \f$, including BI.
 * @param N Number of particles.
 * @param mls_imove Type of particles affected
 */
__kernel void mls_bi_inv(const __global int* imove,
                      __global matrix* mls_bi,
                      const __global float* shepard,
                      usize N,
                      uint mls_imove)
{
    const usize i = get_global_id(0);
    const usize it = get_local_id(0);
    if(i >= N)
        return;
    if(imove[i] != mls_imove){
        return;
    }

    float shepard_i = shepard[i];
    if(shepard_i < 1.0E-6f){
        // It will be considered that there are not enough
        // particles to interpolate
        shepard_i = 1.f;
    }
	mls_bi[i] /= shepard_i;

        mls_bi[i] = MATRIX_INV(mls_bi[i]);
}

/*
 * @}
 */
