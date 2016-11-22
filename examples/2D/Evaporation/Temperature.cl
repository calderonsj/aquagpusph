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
    #include "/home/calderon/SPH/aquagpusph.cmake/resources/Scripts/types/2D.h"
    #include "/home/calderon/SPH/aquagpusph.cmake/resources/Scripts/KernelFunctions/Wendland2D.hcl"
#else
    #include "/home/calderon/SPH/aquagpusph.cmake/resources/Scripts/types/3D.h"
    #include "/home/calderon/SPH/aquagpusph.cmake/resources/Scripts/KernelFunctions/Wendland3D.hcl"
#endif

__kernel void predictor(__global int* imove,
                        __global float* T,
                        __global float* dTdt,
		        __global float* T_in,
		        __global float* dTdt_in,
                        unsigned int N,
                        float dt)
{
    // find position in global arrays
    unsigned int i = get_global_id(0);
    if(i >= N)
        return;

    float DT = dt;
    if(imove[i] <= 0)
    // Same boundary condition as in basic predictor 
        DT = 0.f;
	return;

    dTdt_in[i] = dTdt[i];
    T_in[i] = T[i] + DT * dTdt[i];
    // Guess that this means temperature does not evolve for boundaries
}


__kernel void corrector(__global int* imove,
                        __global float* T,
                        __global float* dTdt,
                        __global float* T_in,
                        __global float* dTdt_in,
                        unsigned int N,
                        float dt)
{
    // find position in global arrays
    unsigned int i = get_global_id(0);
    if(i >= N)
        return;
    if(imove[i] != 1)
        return;

    T[i] += 0.5f * dt * (dTdt[i] - dTdt_in[i]);
    T_in[i] = T[i];
    dTdt_in[i] = dTdt[i];
}

__kernel void   sort(const __global float *T_in, __global float *T,
                     const __global float *dTdt_in, __global float *dTdt,
                     const __global unit *id_sorted,
                     unsigned int N)
{
    uint i = get_global_id(0);
    if(i >= N)
        return;

    const uint i_out = id_sorted[i];

    T[i_out] = T_in[i];
    dTdt[i_out] = dTdt_in[i];
}


__kernel void  lap(const __global int* imove,
                   const __global vec* r,
                   const __global float* rho,
                   const __global float* m,
                   const __global float* T,
                   __global float* lap_T,
                   const float kc,
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
    if(imove[i] != 1){
        return;
    }

    const vec_xyz r_i = r[i].XYZ;
    const float rho_i = rho[i];
    const float T_i = T[i];

    // Initialize the output
    #ifndef LOCAL_MEM_SIZE
        #define _LAPT_ lap_T[i]
    #else
        #define _LAPT_ lap_T_l[it]
        __local float lap_T_l[LOCAL_MEM_SIZE];
        _LAPT_ = 0.f;
    #endif

    BEGIN_LOOP_OVER_NEIGHS(){
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
            const float m_j = m[j];
            const float T_j = T[j];
            const float f_ij = kernelF(q) * CONF * m[j];

            _LAPT_ += 2.f * kc * (T_j - T_i) * f_ij / (rho_i * rho_j);
        }
    }END_LOOP_OVER_NEIGHS()

    #ifdef LOCAL_MEM_SIZE
        lap_T[i] = _LAPT_;
    #endif
}

__kernel void  pdivu(const __global int* imove,
                     const __global vec* r,
                     const __global float* rho,
                     const __global float* m,
                     const __global float* p,
                     const __global float* div_u,
                     __global float* pdiv_u,
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
    if(imove[i] != 1){
        return;
    }

    const vec_xyz r_i = r[i].XYZ;
    const float rho_i = rho[i];
    const float p_i = p[i];
    const float m_i = m[i];
    const float div_u_i = div_u[i];

    // Initialize the output
    #ifndef LOCAL_MEM_SIZE
        #define _PDIVU_ pdiv_u[i]
    #else
        #define _PDIVU_ pdiv_u_l[it]
        __local float pdiv_u_l[LOCAL_MEM_SIZE];
        _PDIVU_ = 0.f;
    #endif

    BEGIN_LOOP_OVER_NEIGHS(){
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
            const float m_j = m[j];
            const float p_j = p[j];
            const float f_ij = kernelF(q) * CONF * m[j];

            _PDIVU_ += p_j * div_u_i * f_ij / rho_j;
        }
    }END_LOOP_OVER_NEIGHS()

    #ifdef LOCAL_MEM_SIZE
        pdiv_u[i] = _PDIVU_;
    #endif
}


// Here, rate of variation of dTdt will be computed

__kernel void rates(const __global uint* iset,
                    const __global int* imove,
                    const __global float* rho,
                    const __global float* m,
                    const __global float* lap_T,
                    const __global float* pdiv_u,
                    __global float* dTdt,
                    const float cp,
                    unsigned int N,
                    vec g)
{
    unsigned int i = get_global_id(0);
    if(i >= N)
        return;
    if(imove[i] != 1)
        return;

    // Temperature equation
      dTdt[i] = (lap_T[i] + pdiv_u[i]) / cp ;

}