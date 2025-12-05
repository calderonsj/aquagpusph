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

/** @file
 * @brief Extra viscous terms computed for the acoustic damper method.
 */

#if defined(LOCAL_MEM_SIZE) && defined(NO_LOCAL_MEM)
    #error NO_LOCAL_MEM has been set.
#endif

#include "resources/Scripts/types/types.h"
#include "resources/Scripts/KernelFunctions/Kernel.h"

/** @brief Compute the acoustic damper viscous force term:
 *           F^{ad} = \lambda^{ad} \nabla (\nabla \cdot \vec{u}(\vec{r}))
 *
 * Compute the acoustic damper viscous force term.
 *
 * @param imove Moving flags.
 *   - imove = 1 for regular fluid particles.
 *   - imove = 0 for sensors.
 *   - imove = -3 for boundary elements/particles (BI method).
 * @param iset Set of particles index.
 * @param r Position \f$ \mathbf{r} \f$.
 * @param rho Density \f$ \rho \f$.
 * @param m Mass \f$ m \f$.
 * @param fi_ad Acoustic damper force \f$ \mathbf{F}^{ad} \f$.
 * @param div_u Velocity divergence \f$ \nabla \cdot \mathbf{u} \f$.
 * @param lambda_ad Acoustic damper coefficient \f$ \lambda^{ad} \f$.
 * @param jhoc Head of chain for each cell (first particle found).
 * @param N Number of particles.
 */
__kernel void acousticTerm(const __global int* imove,
                    const __global int* iset,
                    const __global vec* r,
                    const __global float* rho,
                    const __global float* m,
                    __global vec* fi_ad,
                    const __global float* div_u,
                    __constant float* lambda_ad,
                    // Link-list data
                    const __global svec2* jhoc,
                    // Simulation data
                    const usize N)
{
    const usize i = get_global_id(0);
    const usize it = get_local_id(0);
    if(i >= N)
        return;
    if(imove[i] != 1){
        return;
    }

    const vec_xyz r_i = r[i].XYZ;
    const float div_u_i = div_u[i];
    const float rho_i = rho[i];

    // Initialize the output
    #ifndef LOCAL_MEM_SIZE
        #define _FIAD_ fi_ad[i].XYZ
    #else
        #define _FIAD_ fi_ad_l[it]
        __local vec_xyz fi_ad_l[LOCAL_MEM_SIZE];
        _FIAD_ = VEC_ZERO.XYZ;
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
        {
            const float rho_j = rho[j];
            const float div_u_j = div_u[j];
            const float f_ij = kernelF(q) * CONF * m[j];

	        _FIAD_ += lambda_ad[iset[i]] * (div_u_j + div_u_i) * f_ij * r_ij / (rho_i * rho_j);	

        }
    }END_FOR_NEIGHS()

    #ifdef LOCAL_MEM_SIZE
        fi_ad[i].XYZ = _FIAD_;
    #endif
}

/** @brief Variation rates for acoustic damper term.
 *
 * @param iset Set of particles index.
 * @param imove Moving flags.
 *   - imove = 1 for regular fluid particles.
 *   - imove = 0 for sensors (ignored by this preset).
 *   - imove < 0 for boundary elements/particles.
 * @param rho Density \f$ \rho \f$.
 * @param fi_ad Acoustic damper force \f$ \mathbf{F}^{ad} \f$.
 * @param shepard Shepard renormalization factor.
 * @param dudt Velocity rate of change \f$ \frac{d \mathbf{u}}{d t} \f$.
 * @param refd Density of reference of the fluid \f$ \rho_0 \f$.
 * @param N Number of particles.
 * @param h Smoothing length \f$ h \f$.
 * @param cs Speed of sound \f$ c_s \f$.
 */
__kernel void rates(const __global unsigned int* iset,
                    const __global int* imove,
                    const __global float* rho,
                    const __global vec* fi_ad,
                    const __global float* shepard,
                    __global vec* dudt,
                    __constant float* refd,
                    const usize N,
					float h,
					float cs)
{
    const usize i = get_global_id(0);
    if(i >= N)
        return;
    if(imove[i] != 1)
        return;

    const uint set_i = iset[i];
    const float rho_i = rho[i];

	dudt[i] += fi_ad[i] / (rho_i * shepard[i]);
}

/*
 * @}
 */