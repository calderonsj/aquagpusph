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

#ifndef HAVE_3D
    #include "../types/2D.h"
    #include "../KernelFunctions/Wendland2D.hcl"
#else
    #include "../types/3D.h"
    #include "../KernelFunctions/Wendland3D.hcl"
#endif

/** @brief Gradient of gamma uneveness function computation.
 *
 * \f[ \nabla \Gamma(\mathbf{x}) = \int_{\Omega}
 *     \nabla W(\mathbf{y} - \mathbf{x}) \mathrm{d}\mathbf{x} \f]
 *
 * The gradient of the gamma function is applied to account for the 
 * uneveness of the particles and achieve a regular initial distribution
 * that does not end up in spurious movements.
 *
 * In the gradient of the gamma function computation the fluid extension
 * particles are not taken into account.
 *
 * @param imove Moving flags.
 *   - imove > 0 for regular fluid particles.
 *   - imove = 0 for sensors.
 *   - imove < 0 for boundary elements/particles.
 * @param r Position \f$ \mathbf{r} \f$.
 * @param m Mass \f$ m \f$.
 * @param grad_gamma Gradient of the gamma term.
 * \f$ \nabla \gamma(\mathbf{x}) = \int_{\Omega}
 *    \nabla W(\mathbf{y} - \mathbf{x}) \mathrm{d}\mathbf{x} \f$.
 * @param icell Cell where each particle is located.
 * @param ihoc Head of chain for each cell (first particle found).
 * @param N Number of particles.
 * @param n_cells Number of cells in each direction
 */
__kernel void entry(const __global int* imove,
                    const __global vec* r,
                    const __global float* m,
                    const __global float* rho,
                    __global vec* grad_gamma,
                    // Link-list data
                    const __global uint *icell,
                    const __global uint *ihoc,
                    // Simulation data
                    uint N,
                    uivec4 n_cells)
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
        #define _DGAMMA_ grad_gamma[i]
    #else
        #define _DGAMMA_ grad_gamma_l[it]
        __local vec grad_gamma_l[LOCAL_MEM_SIZE];
        _DGAMMA_ = VEC_ZERO;
    #endif

    BEGIN_LOOP_OVER_NEIGHS(){
        if(imove[j] != 1){
            j++;
            continue;
        }

        const vec_xyz r_ij = r[j].XYZ - r_i;
        const float q = length(r_ij) / H;
        if(q >= SUPPORT)
        {
            j++;
            continue;
        }
	// It is NOT renormalized by H
        {
            _DGAMMA_.XYZ += r_ij * kernelF(q) * CONF * m[j] / rho[j];
        }
    }END_LOOP_OVER_NEIGHS()

    #ifdef LOCAL_MEM_SIZE
        grad_gamma[i] = _DGAMMA_;
    #endif
}
