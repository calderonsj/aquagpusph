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

/** Quasi-second order time integration predictor stage.
 */
__kernel void main(__global int* imove,
                   __global unsigned int* iset,
                   __global vec* pos,
                   __global vec* v,
                   __global vec* dvdt,
                   __global float* rho,
                   __global float* drhodt,
                   __global vec* pos_in,
                   __global vec* v_in,
                   __global vec* dvdt_in,
                   __global float* rho_in,
                   __global float* drhodt_in,
                   __global float* p_in,
                   __constant float* gamma,
                   __constant float* refd,
                   unsigned int N,
                   float dt,
                   float cs,
                   vec g,
                   float rho_min,
                   float rho_max)
{
	// find position in global arrays
	unsigned int i = get_global_id(0);
	if(i >= N)
		return;

	// Momentum equation is solved just for the fluid particles
	float DT = dt;
	if(imove[i] <= 0)
		DT = 0.f;
	v_in[i] = v[i] + DT * (dvdt[i] + g);
	pos_in[i] = pos[i] + DT * v[i] + 0.5f * DT * DT * (dvdt[i] + g);
    
	// Continuity equation must be solved for the fixed particles too
	if(imove[i] != -1){
        DT = dt;
    }
    rho_in[i] = rho[i] + DT * drhodt[i];
	if(rho_in[i] < rho_min) rho_in[i] = rho_min;
	if(rho_in[i] > rho_max) rho_in[i] = rho_max;

	// Batchelor 1967
	{
		const float ddenf = rho_in[i] / refd[iset[i]];
		const float prb = cs * cs * refd[iset[i]] / gamma[iset[i]];
		p_in[i] = prb * (pow(ddenf, gamma[iset[i]]) - 1.f);
	}

	dvdt_in[i] = VEC_ZERO;
	drhodt_in[i] = 0.f;
    
    /*
	// Predictor step for the fluid and walls
	v[i] = v_in[i] + DT * (dvdt_in[i] + g);
	pos[i] = pos_in[i] + DT * v_in[i] + 0.5f * DT * DT * (dvdt_in[i] + g);

	if(imove[i] != -1){
		// Continuity equation must be solved for the fixed particles too
        DT = dt;
    }
    rho[i] = rho_in[i] + DT * drhodt_in[i];
	if(rho[i] < rho_min) rho[i] = rho_min;
	if(rho[i] > rho_max) rho[i] = rho_max;

	// Batchelor 1967
	{
		const float ddenf = rho[i] / refd[iset[i]];
		const float prb = cs * cs * refd[iset[i]] / gamma[iset[i]];
		press[i] = prb * (pow(ddenf, gamma[iset[i]]) - 1.f);
	}
	// Output variables reinitialization
	dvdt[i] = VEC_ZERO;
	drhodt[i] = 0.f;
    */
}
