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
 * @brief Boundary element - Fluid particle interaction.
 * (See DeLeffe.cl for details)
 *
 * It is prefearable to use a header to be included instead of generating a
 * function for thye particles interaction, which imply more registries
 * consumption.
 */

// Artificial viscosity factor
#ifndef __CLEARY__
    #ifndef HAVE_3D
        #define __CLEARY__ 8.f
    #else
        #define __CLEARY__ 10.f
    #endif
#endif

// ------------------------------------------------------------------
// face properties
// ------------------------------------------------------------------
const vec_xyz n_j = nmirrored[j].XYZ;  // Assumed outwarding oriented
const float rho_j = rho[j];
if(rho_j <= 0.01f * refd_i){
    j++;
    continue;
}
const float area_j = m[j];

// ------------------------------------------------------------------
// Boundary element computation
// ------------------------------------------------------------------
{
    const float p_j = p[j];
    const vec_xyz du = umirrored[j].XYZ - u_i;
    const float udn = rho_j * dot(du, n_j);
    //---------------------------------------------------------------
    //       calculate the kernel w_ij
    //---------------------------------------------------------------
    const float w_ij = kernelW(q) * CONW * area_j;
    //---------------------------------------------------------------
    //       calculate the pressure factor
    //---------------------------------------------------------------
    const vec_xyz prfac = rho_j * (prfac_i + p_j / (rho_j * rho_j)) * n_j;
    //---------------------------------------------------------------
    //       calculate viscosity terms
    //---------------------------------------------------------------
    // const float r2 = (q * q + 0.01f) * H * H;
    // const vec_xyz lapufac = __CLEARY__ * vdn / (r2 * rho_i * rho_j) * n_j;
    //---------------------------------------------------------------
    //     Momentum equation (grad(p)/rho and lap(u)/rho)
    //---------------------------------------------------------------
    _GRADP_ += prfac * w_ij;
    // _LAPU_ += lapufac * w_ij;
    //---------------------------------------------------------------
    //     Continuity equation (rho*div(u))
    //---------------------------------------------------------------
    _DIVU_ += udn * w_ij;
    //---------------------------------------------------------------
    //     Density diffusion term (lap(p))
    //---------------------------------------------------------------
    // const float ndr = - rho_j * dot(r_ij, n_j) / r2;
    // const float drfac = (p_j - p_i) - refd_i * dot(g, r_ij);
    // _LAPP_ -= drfac * ndr * w_ij;
}
