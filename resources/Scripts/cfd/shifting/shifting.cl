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

#ifndef __AC_FACTOR__
    /** @def __AC_FACTOR__
     * @brief Anticlumping term factor (R) for shifting.
     * \f$ R \left(\f$, 
     * \frac{W(q)}{W\!\left(\dfrac{\Delta r}{H}\right)}\right)^{4}
     */
    #define __AC_FACTOR__ 0.f
#endif

/** @brief Gradient of shifting formula computation.
 *
 * The gradient of the shifting formula is computed for
 * fluid particles.
 *
 * @param imove Moving flags.
 *   - imove > 0 for regular fluid particles.
 *   - imove = 0 for sensors.
 *   - imove < 0 for boundary elements/particles.
 * @param r Position \f$ \mathbf{r} \f$.
 * @param rho Density \f$ \rho \f$.
 * @param m Mass \f$ m \f$.
 * @param gradC Gradient of the shifting formula.
 * @param dr Characteristic distance between particles.
 * @param jhoc Head of chain for each cell (first particle found).
 * @param N Number of particles.
 */
__kernel void drfluid(const __global int* imove,
                    const __global vec* r,
                    const __global float* rho,
                    const __global float* m,
                    __global vec* gradC,
  		    float dr,
                    // Link-list data
                    const __global svec2* jhoc,
                    // Simulation data
                    const usize N)
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
        #define _GRADC_ gradC[i]
    #else
        #define _GRADC_ gradC_l[it]
        __local vec gradC_l[LOCAL_MEM_SIZE];
        _GRADC_ = VEC_ZERO;
    #endif

    const float dr_over_H = dr / H;
    const float Wref = kernelW(dr_over_H);

    FOR_NEIGHS(N, jhoc){
        if(imove[j] != 1){
            continue;
        }
        const vec_xyz r_ij = r[j].XYZ - r_i;
        const float q = length(r_ij) / H;
		
        if(q >= SUPPORT)
        {
            continue;
        }

        {
		const float f_ij = kernelF(q) * CONF * m[j];
		// --- Anticlumping factor ---
        	const float Wq    = kernelW(q);
        	const float ratio = Wq / Wref;          // W(q) / W(dr/H)

        	// ratio^4 without pow()
        	const float ratio2 = ratio * ratio;
        	const float ratio4 = ratio2 * ratio2;

		// Anticlumping term factor R reduced to R <= 0.05 to keep stability
        	const float Aac = 1.0f + __AC_FACTOR__ * ratio4; // 1 + R * (W(q)/W(dr/H))^4

        	// Apply to shifting gradient contribution
        	_GRADC_ += Aac * r_ij * f_ij / rho[j];
	}
    }END_FOR_NEIGHS()

    #ifdef LOCAL_MEM_SIZE
	gradC[i] = _GRADC_ ;
    #endif
}

/** @brief Final shifting formula computation.
 *
 * The shifting formula is computed including 
 * corrections and the maximum variation.
 *
 * @param imove Moving flags.
 *   - imove > 0 for regular fluid particles.
 *   - imove = 0 for sensors.
 *   - imove < 0 for boundary elements/particles.
 * @param gradC Gradient of the shifting formula.
 * @param du_shift Shifting velocity.
 * @param Uref Characteristic velocity factor for shifting.
 * @param jhoc Head of chain for each cell (first particle found).
 * @param N Number of particles.
 */

__kernel void du(__global int* imove,
                    __global vec* du_shift,
                    const __global vec* gradC,
		    float Uref,
                    usize N)
{
    usize i = get_global_id(0);
    if(i >= N)
        return;

    if(imove[i] != 1)
        return;

	#if __SHIFTING_FORMULATION__ == __SUN_SHIFTING__
            const vec_xyz shift = - Uref * (2.f * H) * gradC[i];
            const float max_shift = 0.5f * Uref;
        #else
            #error Unknown Shifting formulation: __SHIFTING_FORMULATION__
        #endif
	
    	const float norm_du = length(shift);
	const float M = min(max_shift, norm_du);

	du_shift[i] = M * fast_normalize(shift);
}

/** @brief Corrector to include shifting.
 *
 * The shifting formula is added to the position.
 *
 * @param imove Moving flags.
 *   - imove > 0 for regular fluid particles.
 *   - imove = 0 for sensors.
 *   - imove < 0 for boundary elements/particles.
 * @param r Position \f$ \mathbf{r} \f$.
 * @param du_shift Shifting velocity.
 * @param dt Time step.
 * @param N Number of particles.
 */

__kernel void correct(__global int* imove,
                    __global vec* r,
		    const __global vec* du_shift,
		    const float dt,
                    usize N)
{
    usize i = get_global_id(0);
    if(i >= N)
        return;

    if(imove[i] != 1)
        return;

    float DT = dt;

    r[i] += DT * du_shift[i];
}

