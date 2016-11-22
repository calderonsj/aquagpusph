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
 * @brief Boundary integral term computation over the fluid particles.
 */

#if defined(LOCAL_MEM_SIZE) && defined(NO_LOCAL_MEM)
    #error NO_LOCAL_MEM has been set.
#endif

#ifndef HAVE_3D
    #include "/home/calderon/SPH/aquagpusph.cmake/resources/Scripts/types/2D.h"
    #include "/home/calderon/SPH/aquagpusph.cmake/resources/Scripts/KernelFunctions/Wendland2D.hcl"
#else
    #include "/home/calderon/SPH/aquagpusph.cmake/resources/Scripts/types/3D.h"
    #include "/home/calderon/SPH/aquagpusph.cmake/resources/Scripts/KernelFunctions/Wendland3D.hcl"
#endif

/** @brief Performs the boundary effect on the fluid particles.
 * @param iset Set of particles index.
 * @param imove Moving flags.
 *   - imove > 0 for regular fluid particles.
 *   - imove = 0 for sensors.
 *   - imove < 0 for boundary elements/particles.
 * @param r Position \f$ \mathbf{r} \f$.
 * @param normal Normal \f$ \mathbf{n} \f$.
 * @param rho Density \f$ \rho \f$.
 * @param m Area of the boundary element \f$ s \f$.
 * @param grad_p Pressure gradient \f$ \nabla p \f$.
 * @param icell Cell where each particle is located.
 * @param ihoc Head of chain for each cell (first particle found).
 * @param N Number of particles.
 * @param n_cells Number of cells in each direction
 */
__kernel void deltashepardSymBI(const __global int* imove,
                                const __global int* imirrored,
                                const __global vec* r,
                                const __global vec* rmirrored,
                                const __global vec* normal,
                                const __global float* rho,
                                const __global float* m,
                                __global vec* dshepard,
                                // Link-list data
                                __global uint *icell,
                                __global uint *ihoc,
                                // Simulation data
                                uint N,
                                uivec4 n_cells)

{
    const uint i = get_global_id(0);
    const uint it = get_local_id(0);
    if(i >= N)
        return;
    if((!imirrored[i]) || (imove[i] != 1))
        return;

    const vec_xyz r_i = r[i].XYZ;
    const float rho_i = rho[i];

    // Initialize the output
    #ifndef LOCAL_MEM_SIZE
        #define _DSHEPARD_ dshepard[i]
    #else
        #define _DSHEPARD_ dshepard_l[it]
        __local vec dshepard_l[LOCAL_MEM_SIZE];
        _DSHEPARD_ = dshepard[i];
    #endif

    BEGIN_LOOP_OVER_NEIGHS(){
		if((!imirrored[j]) || ((imove[j] != -3))){ // && (imove[j] != -1))){
			j++;
			continue;
        }
		const vec_xyz r_ij = rmirrored[j].XYZ - r_i;
        const float q = length(r_ij) / H;
        if(q >= SUPPORT)
        {
            j++;
            continue;
        }

        {
            const vec_xyz n_j = normal[j].XYZ;  // Assumed outwarding oriented

            _DSHEPARD_.XYZ -= kernelW(q) * CONW * H * n_j * m[j];
        }
    }END_LOOP_OVER_NEIGHS()

    #ifdef LOCAL_MEM_SIZE
        dshepard[i] = _DSHEPARD_;
    #endif
}
