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

/** @brief Fluid particles interactions computation for extra du terms.
 *
 * Compute the differential operators involved in the delta-plus numerical scheme, taking
 * into account just the fluid-fluid interactions.
 *
 * @param imove Moving flags.
 *   - imove > 0 for regular fluid particles.
 *   - imove = 0 for sensors.
 *   - imove < 0 for boundary elements/particles.
 * @param r Position \f$ \mathbf{r} \f$.
 * @param u Velocity \f$ \mathbf{u} \f$.
 * @param rho Density \f$ \rho \f$.
 * @param m Mass \f$ m \f$.
 * @param jhoc Head and tail of chains for each cell.
 * @param du_shift Shifting velocity.
 * @param div_du Divergence of shifting velocity.
 * @param div_drho Divergence of the product of density and shifting velocity.
 * @param div_mu Divergence of the outer product of velocity and shifting velocity.
 * @param N Number of particles.
 */
__kernel void extra(const __global int* imove,
                    const __global vec* r,
                    const __global float* rho,
                    const __global float* m,
			const __global vec* u,
			const __global vec* du_shift,
                    __global float* div_du,
                    __global float* div_drho,
                    __global vec* div_mu,
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
	const vec_xyz u_i = u[i].XYZ;
	const vec_xyz du_i = du_shift[i].XYZ;
    const float rho_i = rho[i];

    // Initialize the output
    #ifndef LOCAL_MEM_SIZE
        #define _DIVDU_ div_du[i]
        #define _DIVDR_ div_drho[i]
        #define _DIVMU_ div_mu[i].XYZ
    #else
        #define _DIVDU_ div_du_l[it]
        #define _DIVDR_ div_drho_l[it]
        #define _DIVMU_ div_mu_l[it]
        __local float div_du_l[LOCAL_MEM_SIZE];
        __local float div_drho_l[LOCAL_MEM_SIZE];
        __local vec_xyz div_mu_l[LOCAL_MEM_SIZE];
        _DIVDU_ = 0.f;
        _DIVDR_ = 0.f;
        _DIVMU_ = VEC_ZERO;
    #endif

    FOR_NEIGHS(N, jhoc){
        if(imove[j] != 1){
            j++;
            continue;
        }

        const vec_xyz r_ij = r[j].XYZ - r_i;
    	const float rho_j = rho[j];
        const float q = length(r_ij) / H;
        const float f_ij = kernelF(q) * CONF * m[j];
		const vec_xyz u_j = u[j].XYZ;
		const vec_xyz du_j = du_shift[j].XYZ;
		const float dudr = dot(du_j - du_i, r_ij);
		
        if(q >= SUPPORT)
        {
            j++;
            continue;
        }

        {
		_DIVDU_ += dudr * f_ij / rho_j;
		_DIVDR_ += dot((rho_j * du_j + rho_i * du_i), r_ij) * f_ij / rho_j;
		_DIVMU_ += MATRIX_DOT((outer(u_j, du_j) + outer(u_i, du_i)) , f_ij * r_ij / rho_j).XYZ;
        }
    }END_FOR_NEIGHS()

    #ifdef LOCAL_MEM_SIZE
        div_du[i] = _DIVDU_;
        div_drho[i] = _DIVDR_;
        div_mu[i] = _DIVMU_;
    #endif

}

/** @brief Adding shifting extra terms to governing equations
 *
 * @param imove Moving flags.
 *   - imove > 0 for regular fluid particles.
 *   - imove = 0 for sensors.
 *   - imove < 0 for boundary elements/particles.
 * @param u Velocity \f$ \mathbf{u} \f$.
 * @param rho Density \f$ \rho \f$.
 * @param div_du Divergence of shifting velocity.
 * @param div_drho Divergence of the product of density and shifting velocity.
 * @param div_mu Divergence of the outer product of velocity and shifting velocity.
 * @param dudt Velocity rate of change
 * \f$ \left. \frac{d \mathbf{u}}{d t} \right\vert_{n+1} \f$.
 * @param drhodt Density rate of change
 * \f$ \left. \frac{d \rho}{d t} \right\vert_{n+1} \f$.
 * @param N Number of particles.
 */
__kernel void rates(const __global uint* iset,
                    const __global int* imove,
                    const __global vec* u,
                    const __global float* rho,
                    const __global float* div_du,
                    const __global float* div_drho,
                    const __global float* div_mu,
                    __global float* drhodt,
		    __global vec* dudt,
                    unsigned int N)
{
    unsigned int i = get_global_id(0);
    if(i >= N)
        return;
    if(imove[i] != 1)
        return;

 	   drhodt[i] += - rho[i] * div_du[i] + div_drho[i];
	   dudt[i]   += div_mu[i] - u[i] * div_du[i];


}