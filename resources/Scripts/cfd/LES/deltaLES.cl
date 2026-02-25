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

#ifndef EXCLUDED_PARTICLE
    /** @brief Condition to exclude a particle from the delta-SPH model
     * 
     * By default all the boundary elements are excluded. Even though it is
     * enough for simulation where fluid and solid mechanics are not combined,
     * it is strongly recommended to conveniently overload this macro. 
     * @note Redefining this macro this OpenCL script can be recicled
     */
    #define EXCLUDED_PARTICLE(index) imove[index] <= 0
#endif

#if defined(LOCAL_MEM_SIZE) && defined(NO_LOCAL_MEM)
    #error NO_LOCAL_MEM has been set.
#endif

#include "resources/Scripts/types/types.h"
#include "resources/Scripts/KernelFunctions/Kernel.h"

#if __LAP_FORMULATION__ == __LAP_MONAGHAN__
    #ifndef HAVE_3D
        #define __CLEARY__ 8.f
    #else
        #define __CLEARY__ 10.f
    #endif
#endif

/** @brief Strain-rate tensor computation.
 *
 * @param imove Moving flags.
 *   - imove > 0 for regular fluid particles.
 *   - imove = 0 for sensors.
 *   - imove < 0 for boundary elements/particles.
 * @param r Position \f$ \mathbf{r} \f$.
 * @param rho Density \f$ \rho \f$.
 * @param m Mass \f$ m \f$.
 * @param u Velocity \f$ \mathbf{u} \f$.
 * @param mls Kernel MLS transformation matrix \f$ L \f$.
 * @param D strain-rate tensor matrix
 * @param jhoc Head and tail of chains for each cell.
 * @param N Number of particles.
 */
__kernel void diff(const __global int* restrict imove,
                   const __global vec* restrict r,
                   const __global float* restrict rho,
                   const __global float* restrict m,
                   const __global vec* restrict u,
		   const __global matrix* restrict mls,
                   __global matrix* restrict D,
                   const __global svec2* restrict jhoc,
                   usize N)
{
    const usize i = get_global_id(0);
    const usize it = get_local_id(0);
    if(i >= N)
        return;
    if(EXCLUDED_PARTICLE(i)){
        return;
    }

    const vec_xyz r_i = r[i].XYZ;
    const vec_xyz u_i = u[i].XYZ;

    // Initialize the output
    #ifndef LOCAL_MEM_SIZE
        #define _D_ D[i]
    #else
        #define _D_ D_l[it]
        __local matrix D_l[LOCAL_MEM_SIZE];
        _D_ = MAT_ZERO;
    #endif

    FOR_NEIGHS(N, jhoc){
        if( (i == j) || (EXCLUDED_PARTICLE(j))){
            continue;
        }

        const vec_xyz r_ij = r[j].XYZ - r_i;
	const vec_xyz u_ij = u[j].XYZ - u_i;
        const float q = length(r_ij) / H;
	const float V_j = m[j] / rho[j];

        if(q >= SUPPORT)
        {
            continue;
        }
        {
            const float f_ij = kernelF(q) * CONF;
	    const vec_xyz mls_i = MATRIX_DOT(mls[i], 0.5f * f_ij * r_ij * V_j);

	    _D_ += outer(u_ij, mls_i) + outer(mls_i, u_ij);
        }
    }END_FOR_NEIGHS()

    #ifdef LOCAL_MEM_SIZE
        D[i] = _D_;
    #endif
}

/** @brief Norm of strain-rate tensor.
 *
 * @param imove Moving flags.
 *   - imove > 0 for regular fluid particles.
 *   - imove = 0 for sensors.
 *   - imove < 0 for boundary elements/particles.
 * @param D strain-rate tensor matrix.
 * @param absD norm of the strain-rate tensor matrix.
 * @param N Number of particles.
 */
__kernel void absDiff(const __global int* restrict imove,
                      const __global matrix* restrict D,
		      __global float* restrict absD,
                      usize N)
{
    const usize i = get_global_id(0);

    if(i >= N)
        return;
    if(EXCLUDED_PARTICLE(i)){
        return;
    }

    const float Diff_i = doubledot(D[i],D[i]);
	
    absD[i] = sqrt(2.0 * Diff_i);

}

/** @brief Local delta and alpha values for every particle.
 *	These are computed from characteristic filtering length $L_{LES}$
 *
 * @param imove Moving flags.
 *   - imove > 0 for regular fluid particles.
 *   - imove = 0 for sensors.
 *   - imove < 0 for boundary elements/particles.
 * @param absD norm of the strain-rate tensor matrix.
 * @param Ldelta local value of delta.
 * @param Lalpha local value of alpha.
 * @param LLES filtering length.
 * @param cs speed of sound.
 * @param N Number of particles.
 */
__kernel void deltaalpha(const __global int* restrict imove,
			 const __global float* restrict absD,
                   	 usize N,
			 __global float* restrict Ldelta,
			 __global float* restrict Lalpha,
			 const float LLes,
			 const float cs)
{
    const usize i = get_global_id(0);
    
    if(i >= N)
        return;
    if(EXCLUDED_PARTICLE(i)){
        return;
    }

	const float mu_d = (1.5 * LLes) * (1.5 * LLes) * absD[i];
	const float mu_a = (0.12 * LLes) * (0.12 * LLes) * absD[i];

	// Limit values for impact flows. More relevant for coarse resolutions.
	Ldelta[i] = fmin(0.5f, mu_d / (cs * H));
	Lalpha[i] = fmin(0.05f, mu_a / (cs * H));
}

/** @brief Fluid particles interactions computation.
 *
 * Compute the differential operators involved in the numerical scheme, taking
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
 * @param p Pressure \f$ p \f$.
 * @param grad_p Pressure gradient \f$ \frac{\nabla p}{rho} \f$.
 * @param lap_u Velocity laplacian \f$ \frac{\Delta \mathbf{u}}{rho} \f$.
 * @param div_u Velocity divergence \f$ \rho \nabla \cdot \mathbf{u} \f$.
 * @param icell Cell where each particle is located.
 * @param ihoc Head of chain for each cell (first particle found).
 * @param N Number of particles.
 * @param n_cells Number of cells in each direction
 */
__kernel void alphaij(const __global uint* restrict iset,
		    const __global int* restrict imove,
                    const __global vec* restrict r,
                    const __global vec* restrict u,
                    const __global float* restrict rho,
                    const __global float* restrict m,
                    const __global float* restrict p,
		    const __global float* restrict Lalpha,
                    __global vec* restrict lap_u_turb,
		    float cs,
		    __constant float* restrict refd,
		    // Link-list data
                    const __global svec2* restrict jhoc,
                    // Simulation data
                    usize N)
{
    const usize i = get_global_id(0);
    const usize it = get_local_id(0);
    if(i >= N)
        return;
    if(imove[i] != 1){
        return;
    }

    const vec_xyz r_i = r[i].XYZ;
    const vec_xyz u_i = u[i].XYZ;
    const float rho_i = rho[i];

    // Initialize the output
    #ifndef LOCAL_MEM_SIZE
        #define _LAPUTURB_ lap_u_turb[i].XYZ
    #else
        #define _LAPUTURB_ lap_u_turb_l[it]
        __local vec_xyz lap_u_turb_l[LOCAL_MEM_SIZE];
        _LAPUTURB_ = VEC_ZERO.XYZ;
    #endif

    FOR_NEIGHS(N, jhoc){
        if(i == j){
            j++;
            continue;
        }
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
        {
            const float rho_j = rho[j];
            const float udr = dot(u[j].XYZ - u_i, r_ij);
            const float f_ij = kernelF(q) * CONF * m[j];		
	    
	    float alpha_ij = 2.f * Lalpha[i] * Lalpha[j] / (Lalpha[i] + Lalpha[j] + 1e-08f);
	    float mu_t = alpha_ij * H * cs * refd[iset[i]];

            #if __LAP_FORMULATION__ == __LAP_MONAGHAN__
                const float r2 = (q * q + 0.01f) * H * H;
                _LAPUTURB_ += mu_t * __CLEARY__ * f_ij * udr / (r2 * rho_i * rho_j) * r_ij;
            #elif __LAP_FORMULATION__ == __LAP_MORRIS__
                _LAPUTURB_ += mu_t * f_ij * 2.f / (rho_i * rho_j) * (u[j].XYZ - u_i);
            #else
                #error Unknown Laplacian formulation: __LAP_FORMULATION__
            #endif
        }
    }END_FOR_NEIGHS()

    #ifdef LOCAL_MEM_SIZE
        lap_u_turb[i].XYZ = _LAPUTURB_;
    #endif
}

/** @brief Laplacian of the pressure computation.
 *
 * @param imove Moving flags.
 *   - imove > 0 for regular fluid particles.
 *   - imove = 0 for sensors.
 *   - imove < 0 for boundary elements/particles.
 * @param r Position \f$ \mathbf{r} \f$.
 * @param rho Density \f$ \rho \f$.
 * @param m Mass \f$ m \f$.
 * @param jhoc Head and tail of chains for each cell.
 * @param lap_p Pressure laplacian \f$ \Delta p \f$.
 * @param N Number of particles.
 */
__kernel void deltaij1(const __global int* restrict imove,
                   const __global vec* restrict r,
                   const __global float* restrict rho,
                   const __global float* restrict m,
                   __global float* restrict lap_p,
		   const __global float* restrict Ldelta,
                   const __global svec2* restrict jhoc,
                   usize N)
{
    const usize i = get_global_id(0);
    const usize it = get_local_id(0);
    if(i >= N)
        return;
    if(EXCLUDED_PARTICLE(i)){
        return;
    }

    const vec_xyz r_i = r[i].XYZ;
    const float rho_i = rho[i];

    // Initialize the output
    #ifndef LOCAL_MEM_SIZE
        #define _LAPP_ lap_p[i]
    #else
        #define _LAPP_ lap_p_l[it]
        __local float lap_p_l[LOCAL_MEM_SIZE];
        _LAPP_ = 0.f;
    #endif

    FOR_NEIGHS(N, jhoc){
        if( (i == j) || (EXCLUDED_PARTICLE(j))){
            continue;
        }
        const vec_xyz r_ij = r[j].XYZ - r_i;
        const float q = length(r_ij) / H;
        if(q >= SUPPORT)
        {
            continue;
        }
        {
            const float f_ij = kernelF(q) * CONF * m[j] / rho[j];
	    const float r2 = (q * q + 0.01f) * H * H;

	    const float delta_ij = 2.f * Ldelta[i] * Ldelta[j] / (Ldelta[i] + Ldelta[j] + 1e-08f);

            _LAPP_ += 2.f * delta_ij * (rho[j] - rho[i]) * dot(r_ij, f_ij * r_ij) / r2;
        }
    }END_FOR_NEIGHS()

    #ifdef LOCAL_MEM_SIZE
        lap_p[i] = _LAPP_;
    #endif
}

__kernel void deltaij2(const __global int* restrict imove,
                        const __global vec* restrict r,
                        const __global float* restrict rho,
                        const __global float* restrict m,
                        const __global vec* restrict lap_p_corr,
                        __global float* restrict lap_p,
			const __global float* restrict Ldelta,
                        const __global svec2* restrict jhoc,
                        usize N)
{
    const usize i = get_global_id(0);
    const usize it = get_local_id(0);
    if(i >= N)
        return;
    if(EXCLUDED_PARTICLE(i)){
        return;
    }

    const vec_xyz r_i = r[i].XYZ;
    const vec_xyz gradp_i = lap_p_corr[i].XYZ;

    // Initialize the output
    #ifndef LOCAL_MEM_SIZE
        #define _LAPP_ lap_p[i]
    #else
        #define _LAPP_ lap_p_l[it]
        __local float lap_p_l[LOCAL_MEM_SIZE];
        _LAPP_ = lap_p[i];
    #endif

    FOR_NEIGHS(N, jhoc){
        if( (i == j) || (EXCLUDED_PARTICLE(j))){
            continue;
        }
        const vec_xyz r_ij = r[j].XYZ - r_i;
        const float q = length(r_ij) / H;
        if(q >= SUPPORT)
        {
            continue;
        }
        {
            const vec_xyz gradp_ij = lap_p_corr[j].XYZ + gradp_i;
            const float f_ij = kernelF(q) * CONF * m[j] / rho[j];
	    const float r2 = (q * q + 0.01f) * H * H;

	    const float delta_ij = 2.f * Ldelta[i] * Ldelta[j] / (Ldelta[i] + Ldelta[j] + 1e-08f);

            _LAPP_ -= delta_ij * dot(gradp_ij, r_ij) * dot(r_ij, f_ij * r_ij) / r2;
        }
    }END_FOR_NEIGHS()

    #ifdef LOCAL_MEM_SIZE
        lap_p[i] = _LAPP_;
    #endif
}

/** @brief Adding turbulent LES viscosity to rates computation
 *
 * The momentum equation is applied from the already
 * computed differential operators:
 *
 *   - \f$ \frac{\mathrm{d} \mathbf{u}}{\mathrm{d} t} =
 *     - \frac{\nabla p}{rho}
 *     + \frac{\mu}{rho} \Delta \mathbf{u}
 *     + \mathbf{g}\f$
 *   - \f$ \frac{\mathrm{d} \rho}{\mathrm{d} t} =
 *     - \rho \nabla \cdot \mathbf{u}
 *     + \delta \Delta t \frac{\rho_a}{\rho_0} \Delta p\f$
 *
 * @param iset Set of particles index.
 * @param imove Moving flags.
 *   - imove > 0 for regular fluid particles.
 *   - imove = 0 for sensors.
 *   - imove < 0 for boundary elements/particles.
 * @param rho Density \f$ \rho_{n+1} \f$.
 * @param grad_p Pressure gradient \f$ \frac{\nabla p}{rho} \f$.
 * @param lap_u Velocity laplacian \f$ \frac{\Delta \mathbf{u}}{rho} \f$.
 * @param div_u Velocity divergence \f$ \rho \nabla \cdot \mathbf{u} \f$.
 * @param dudt Velocity rate of change
 * \f$ \left. \frac{d \mathbf{u}}{d t} \right\vert_{n+1} \f$.
 * @param drhodt Density rate of change
 * \f$ \left. \frac{d \rho}{d t} \right\vert_{n+1} \f$.
 * @param visc_dyn Dynamic viscosity \f$ \mu \f$.
 * @param N Number of particles.
 * @param g Gravity acceleration \f$ \mathbf{g} \f$.
 */
__kernel void ratesLES(const __global uint* restrict iset,
                    const __global int* restrict imove,
                    const __global vec* restrict lap_u_turb,
                    __global vec* restrict dudt,
                    usize N)
{
    const usize i = get_global_id(0);
    if(i >= N)
        return;
    if(imove[i] != 1)
        return;

    // Momentum equation
    dudt[i] += lap_u_turb[i];
}

__kernel void deltaSPH(const __global int* restrict imove,
                    const __global float* restrict lap_p,
                    __global float* restrict drhodt,
                    usize N,
		    float cs)
{
    const usize i = get_global_id(0);
    if(i >= N)
        return;
    if(EXCLUDED_PARTICLE(i))
        return;

        const float delta_f = H * cs;
	drhodt[i] += delta_f * lap_p[i];
    
}

/*
 * @}
 */