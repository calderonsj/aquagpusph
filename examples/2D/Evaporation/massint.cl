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
    #include "/home/calderon/SPH/aquagpusph.cmake/resources/Scripts/types/2D.h"
#else
    #include "/home/calderon/SPH/aquagpusph.cmake/resources/Scripts/types/3D.h"
#endif

__kernel void mint(__global int* imove,
                  __global vec* dshepard,
                  __global float* dshepard_i,
                  __global float* mint,
		          const __global float* m,
		          const __global float* rho,
                  unsigned int N)
{
    // find position in global arrays
    unsigned int i = get_global_id(0);

    if(i >= N)
        return;

    mint[i] = 0.0;

    if(imove[i] == 1){
        dshepard_i[i] = length(dshepard[i]);
        if (dshepard_i[i] > 0.6){
        mint[i] = m [i];
        }
	    return;
    }
}