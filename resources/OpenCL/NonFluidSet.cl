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

#ifndef HAVE_3D
    #include "types/2D.h"
#else
    #include "types/3D.h"
#endif

#ifndef HAVE_3D
	#ifndef NEIGH_CELLS
		/** @def NEIGH_CELLS Number of neighbour cells. In 2D case 8,
		 * and the main cells must be computed, but in 3D 27 cells,
		 * must be computed.
		 */ 
		#define NEIGH_CELLS 9
	#endif
#else
	#ifndef NEIGH_CELLS
		#define NEIGH_CELLS 27
	#endif
#endif

/** To set the values of pressure and density two approaches may be considered:
 *   - _ALL_INTERPOLATED_: Proposed by Ferrand, the density and pressure are interpolated from the neighbour particles.
 *   - _DENSITY_BASED_: Just the density is interpolated, computing the pressure using the EOS.
 *   - _PRESSURE_BASED_: The pressure value is interpolated, and the density is computed the EOS.
 * Pressure field is corrected (in all the cases) by the hydrostatic field.
 * Since density is not corrected, if you are using weakly compressible
 * (where the fluid pressure and density are related) just using the
 * _PRESSURE_BASED_ algorithm may return a consistent density field.
 */
#define _PRESSURE_BASED_


__kernel void main(const __global uint* iset, const __global int* imove,
                   __global vec* grad_p, __global float* div_u,
                   __global float* rho, __global float* p,
                   const __global float* shepard, __constant float* refd,
                   __constant float* gamma, uint N, float cs)
{
    uint i = get_global_id(0);
    if(i >= N)
        return;

    if((imove == -1) || (imove > 0)){
        return;
    }

    float shepard_i = shepard[i];
    if(shepard_i < 0.01f){
        // It will be considered that there are not enough
        // particles to interpolate
        shepard_i = 1.f;
        grad_p[i] = VEC_ZERO;
        div_u[i] = 0.f;
    }

    #ifdef _ALL_INTERPOLATED_
        rho[i] = div_u[i] / shepard_i;
        p[i] = (grad_p[i].x + grad_p[i].y) / shepard_i;
    #elif defined _DENSITY_BASED_
        rho[i] = div_u[i] / shepard_i;
        // Batchelor 1967
        const float rdenf = refd[iset[i]];
        const float gammf = gamma[iset[i]];
        const float ddenf = rho[i] / rdenf;
        const float prb = cs * cs * rdenf / gammf;
        p[i] = max(prb * (pow(ddenf, gammf) - 1.f),
                   grad_p[i].x / shepard_i) + grad_p[i].y / shepard_i;
    #elif defined _PRESSURE_BASED_
        p[i] = max(0.f, (grad_p[i].x + grad_p[i].y) / shepard_i);
        // Reversed Batchelor 1967
        const float rdenf = refd[iset[i]];
        const float gammf = gamma[iset[i]];
        const float iprb = gammf / (cs * cs * rdenf);
        rho[i] = rdenf * pow(1.f + iprb * p[i], 1.f / gammf);
    #else
        #error "Unknow boundary elements field computation algorithm"
    #endif
    // Rinitializate output variables
    grad_p[i] = VEC_ZERO;
    div_u[i] = 0.f;
}
