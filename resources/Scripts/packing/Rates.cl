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
 * @brief Velocity variation rate computation for the particle packing algorithm.
 */

#ifndef HAVE_3D
    #include "../types/2D.h"
#else
    #include "../types/3D.h"
#endif


/** @brief Velocity variation rate computation for the particle packing algorithm.
 *
 * The mass conservation is neglected and velocity rate equation is applied from the already
 * computed differential operators:
 *
 *    \f$ \frac{\mathrm{d} \mathbf{u}}{\mathrm{d} t} =
 *     - \beta \, \nabla \gamma_i
 *     - \zeta \, \mathbf{u}_i
 *   - \f$ \frac{\mathrm{d} \rho}{\mathrm{d} t} =
 *     0.0f$
 *
 * @param imove Moving flags.
 *   - imove > 0 for regular fluid particles.
 *   - imove = 0 for sensors.
 *   - imove < 0 for boundary elements/particles.
 * @param u Velocity \f$ \mathbf{u} \f$.
 * @param shepard Shepard term
 * \f$ \gamma(\mathbf{x}) = \int_{\Omega}
 *     W(\mathbf{y} - \mathbf{x}) \mathrm{d}\mathbf{x} \f$.
 * @param grad_gamma Gradient of the gamma term.
 * \f$ \nabla \gamma(\mathbf{x}) = \int_{\Omega}
 *    \nabla W(\mathbf{y} - \mathbf{x}) \mathrm{d}\mathbf{x} \f$.
 * @param dudt Velocity rate of change
 * \f$ \left. \frac{d \mathbf{u}}{d t} \right\vert_{n+1} \f$.
 * @param beta: pressure gradient constant factor.
 * \f$ \beta = \frac{2 p_0}{\rho_0} \f$.
 * @param zeta: artificial dissipation constant factor.
 * \f$ \zeta = \alpha \frac{\sqrt{\beta}}{V_{0}^{1 / DIMS}}
 * @param N Number of particles.
 */
__kernel void entry(const __global int* imove,
					const __global vec* u,
					const __global float* shepard,
                    const __global vec* grad_gamma,
                    __global vec* dudt,
                    float beta,
                    float zeta,
                    unsigned int N)
{
    unsigned int i = get_global_id(0);
    if(i >= N)
        return;
    if(imove[i] != 1)
        return;
    dudt[i] = -beta * grad_gamma[i] / shepard[i] - zeta * u[i];
}