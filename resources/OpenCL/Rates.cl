/*
 *  This file is part of AQUA-gpusph, a free CFD program based on SPH.
 *  Copyright (C) 2012  Jose Luis Cercos Pita <jl.cercos@upm.es>
 *
 *  AQUA-gpusph is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  AQUA-gpusph is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with AQUA-gpusph.  If not, see <http://www.gnu.org/licenses/>.
 */

// To use Gaussian kernel please compile AQUAgpusph with Gauss kernel option
#ifndef HAVE_3D
	#include "KernelFunctions/Wendland2D.hcl"
#else
	#include "KernelFunctions/Wendland3D.hcl"
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

#ifdef __NO_LOCAL_MEM__
	#define _SIGMA_ sigma[labp]
	#define _F_ f[labp]
	#define _DRDT_ drdt[labp]
	#define _DRDT_F_ drdt_F[labp]
	#define _SHEPARD_ shepard[labp]
	#define _GRADW_ gradW[labp]
#else
	#define _SIGMA_ lSigma[it]
	#define _F_ lF[it]
	#define _DRDT_ lDrdt[it]
	#define _DRDT_F_ lDrdt_F[it]
	#define _SHEPARD_ lShepard[it]
	#define _GRADW_ lGradW[it]
#endif

#ifndef M_PI
	#define M_PI 3.14159265359f
#endif
#ifndef iM_PI
	#define iM_PI 0.318309886f
#endif

#ifndef uint
	#define uint unsigned int
#endif

#ifdef _g
	#error '_g' is already defined.
#endif
#define _g __global

#ifdef _c
	#error '_c' is already defined.
#endif
#define _c __constant

#ifdef _l
	#error '_l' is already defined.
#endif
#define _l __local

/** Sort all input data by cells (the output data will remain unsorted).
 * @param iFluidin Unsorted fluid identifier.
 * @param iFluid Sorted fluid identifier.
 * @param posin Unsorted positions.
 * @param vin Unsorted velocity.
 * @param hpin Unsorted kernels height.
 * @param densin Unsorted density.
 * @param pressin Unsorted pressure.
 * @param pmassin Unsorted mass.
 * @param pos Sorted positions.
 * @param v Sorted velocity.
 * @param hp Sorted kernels height.
 * @param dens Sorted density.
 * @param press Sorted pressure.
 * @param pmass Sorted mass.
 * @param dPermut Sorted space -> unsorted space permutations.
 * @param iPermut Unsorted space -> sorted space permutations.
 * @param N Number of particles.
 */
__kernel void SortData( _g int* iFluidin, _g int* iFluid,
                        _g int* iMovein, _g int* iMove,
                        _g vec* posin, _g vec* vin, _g float* hpin,
                        _g float* densin, _g float* pressin, _g float* pmassin,
                        _g vec* pos, _g vec* v, _g float* hp,
                        _g float* dens, _g float* press, _g float* pmass,
                        _g uint *dPermut, _g uint *iPermut,
                        uint N)
{
	uint i = get_global_id(0);
	if(i >= N)
		return;

	// ---- | ------------------------ | ----
	// ---- V ---- Your code here ---- V ----

	// We assume i in the unsorted space, and labp in the sorted one 
	uint labp = iPermut[i];

	iFluid[labp] = iFluidin[i];
	iMove[labp]  = iMovein[i];
	v[labp]      = vin[i];
	hp[labp]     = hpin[i];
	dens[labp]   = densin[i];
	press[labp]  = pressin[i];
	pmass[labp]  = pmassin[i];
	pos[labp]    = posin[i];

	// ---- A ---- Your code here ---- A ----
	// ---- | ------------------------ | ----
}

/** Compute the rates of variation due to the fluid (fixed particles will be
 * included here). During this stage some other operations are performed as
 * well, like the values interpolation in the boundaries (for DeLeffe boundary
 * conditions), the sensors meassurement, or the Shepard factor computation.
 * @param iFluid Fluid identifier.
 * @param pos Position.
 * @param v Velocity.
 * @param dens Density.
 * @param hp Kernel height.
 * @param pmass Mass
 * @param press Pressure.
 * @param Visckin Kinetic viscosity (one per fluid)
 * @param Viscdyn Dynamic viscosity (one per fluid)
 * @param shepard Shepard term (0th correction).
 * @param gradW Shepard term gradient.
 * @param f Acceleration.
 * @param drdt Rate of change of the density.
 * @param sigma Viscosity time step term.
 * @param shepard Shepard factor.
 * @param gradW Shepard factor gradient.
 * @param lcell Cell where the particle is situated.
 * @param ihoc Head particle of cell chain.
 * @param validCell Mark cells that have at least one fluid particle.
 * @param dPermut Sorted space -> unsorted space permutations.
 * @param iPermut Unsorted space -> sorted space permutations.
 * @param n Number of particles.
 * @param N Number of particles & sensors.
 * @param lvec Number of cells at each direction.
 * @param grav Gravity acceleration.
 */
__kernel void Rates( _g int* iFluid, _g int* iMove,
                     _g vec* pos, _g vec* v, _g float* dens, _g float* hp,
                     _g float* pmass, _g float* press, _c float* Visckin,
                     _c float* Viscdyn, _g vec* f, _g float* drdt,
                     _g float* drdt_F, _g float* sigma, _g float* shepard,
                     _g vec* gradW,
                     // Link-list data
                     _g uint *lcell, _g uint *ihoc, _g short* validCell,
                     _g uint *dPermut, _g uint *iPermut,
                     // Specific sensors data
                     _g ushort *sensorMode,
                     // Simulation data
                     uint n, uint N, uivec lvec, vec grav
                     #ifdef __DELTA_SPH__
                         // Continuity equation diffusive term data
                         , _c float* refd, _c float* delta
                         , float dt, float cs
                     #endif
                     #ifdef __NO_LOCAL_MEM__
                         )
                     #else
                         // Local memory to accelerate writes
                         , _l float* lSigma, _l vec* lF, _l float* lDrdt
                         , _l float* lDrdt_F, _l float* lShepard
                         , _l vec* lGradW)
                     #endif
{
	// find position in global arrays
	uint i = get_global_id(0);			// Particle at sorted space
	uint it = get_local_id(0);			// Particle at local memory (temporal storage)
	if(i >= N)
		return;

	// ---- | ------------------------ | ----
	// ---- V ---- Your code here ---- V ----

	/* All the data has previously been sorted, so two spaces must be
	 * considereed:
	 * # Sorted space, where i is the index of the particle.
	 * # Unsorted space, where labp is the index of the particle.
	 * In the sorted space are stored the input data in order to read
	 * the variables in a convenient coalescing way, while in the unsorted
	 * space are stored the output data to can easily manage the time
	 * integration and output files writting.
	 */

	// Kernel variables
	float hav, dist, conw, conf, wab, fab;
	// Particle data
	int iIFluid, iIMove;
	uint j,labp, lc;
	vec iPos, iV;
	float iHp, iDens, iPress, iVisckin, iViscdyn;
	// Neighbours data
	uint cellCount, lcc;
	vec r,dv;
	float r1, vdr, pDens, pMass, prfac, viscg;

	j = i;              // Backup of the variable, in order to compare later
	labp = dPermut[i];  // Particle index at unsorted space
	lc = lcell[i];      // Cell of the particle
	if(!validCell[lc]){
		// We should not waste time computing boundaries without fluid.
		shepard[labp] = 0.f;
		gradW[labp]   = VEC_ZERO;
		return;
	}
	// Output initialization
	_SIGMA_   = 1000.f;
	_F_       = VEC_ZERO;
	_DRDT_    = 0.f;
	_DRDT_F_  = 0.f;
	_SHEPARD_ = 0.f;
	_GRADW_   = VEC_ZERO;

	// Particle data
	iHp = hp[i];
	iPos = pos[i];
	iV = v[i];
	iDens = dens[i];
	iPress = press[i];
	iIFluid = iFluid[i];
	iVisckin = Visckin[iIFluid];
	iViscdyn = Viscdyn[iIFluid];
	iIMove = iMove[i];
	#ifdef __DELTA_SPH__
		float rDens, iDelta, drfac, rdr;
		rDens  = refd[iIFluid];
		iDelta = delta[iIFluid];
	#endif
	//! 3th.- Loop over all neightbour particles
	{
		//! 3.a.- Home cell, starting by next particle
		i++;
		while( (i<N) && (lcell[i]==lc) ) {
			// Sensor specific computation
			if(!iIMove){
				#include "RatesSensors.hcl"
			}
			else{
				#if __BOUNDARY__==0
					// ElasticBounce boundary condition
					if(iIMove<0){
						#include "RatesBounds.hcl"
					}
					else{
						#include "Rates.hcl"
					}
				#elif __BOUNDARY__==1
					// Fix particles
					#include "Rates.hcl"
				#elif __BOUNDARY__==2
					// DeLeffe boundary condition
					if(iIMove<0){
						#include "RatesBounds.hcl"
					}
					else{
						#include "Rates.hcl"
					}
				#else
					#error Unknow boundary condition
				#endif
			}
			i++;
		}
		//! 3.b.- Neighbour cells
		for(cellCount=1;cellCount<NEIGH_CELLS;cellCount++) {
			// Loop over 8 neighbour cells, taking cell index (lcc)
			switch(cellCount) {
				// Cells at the same Z than main cell
				case 0: lcc = lc + 0; break;
				case 1: lcc = lc + 1; break;
				case 2: lcc = lc - 1; break;
				case 3: lcc = lc + lvec.x; break;
				case 4: lcc = lc + lvec.x + 1; break;
				case 5: lcc = lc + lvec.x - 1; break;
				case 6: lcc = lc - lvec.x; break;
				case 7: lcc = lc - lvec.x + 1; break;
				case 8: lcc = lc - lvec.x - 1; break;
				#ifdef HAVE_3D
					// Cells bellow main cell
					case 9 : lcc = lc + 0          - lvec.x*lvec.y; break;
					case 10: lcc = lc + 1          - lvec.x*lvec.y; break;
					case 11: lcc = lc - 1          - lvec.x*lvec.y; break;
					case 12: lcc = lc + lvec.x     - lvec.x*lvec.y; break;
					case 13: lcc = lc + lvec.x + 1 - lvec.x*lvec.y; break;
					case 14: lcc = lc + lvec.x - 1 - lvec.x*lvec.y; break;
					case 15: lcc = lc - lvec.x     - lvec.x*lvec.y; break;
					case 16: lcc = lc - lvec.x + 1 - lvec.x*lvec.y; break;
					case 17: lcc = lc - lvec.x - 1 - lvec.x*lvec.y; break;
					// Cells over main cell
					case 18: lcc = lc + 0          + lvec.x*lvec.y; break;
					case 19: lcc = lc + 1          + lvec.x*lvec.y; break;
					case 20: lcc = lc - 1          + lvec.x*lvec.y; break;
					case 21: lcc = lc + lvec.x     + lvec.x*lvec.y; break;
					case 22: lcc = lc + lvec.x + 1 + lvec.x*lvec.y; break;
					case 23: lcc = lc + lvec.x - 1 + lvec.x*lvec.y; break;
					case 24: lcc = lc - lvec.x     + lvec.x*lvec.y; break;
					case 25: lcc = lc - lvec.x + 1 + lvec.x*lvec.y; break;
					case 26: lcc = lc - lvec.x - 1 + lvec.x*lvec.y; break;
				#endif
			}
			// Sub-loop over particles into neighbour cells
			i = ihoc[lcc];
			while( (i<N) && (lcell[i]==lcc) ) {
				if(!iIMove){
					#include "RatesSensors.hcl"
				}
				else{
					#if __BOUNDARY__==0
						// ElasticBounce boundary condition
						if(iIMove<0){
							#include "RatesBounds.hcl"
						}
						else{
							#include "Rates.hcl"
						}
					#elif __BOUNDARY__==1
						// Fix particles
						#include "Rates.hcl"
					#elif __BOUNDARY__==2
						// DeLeffe boundary condition
						if(iIMove<0){
							#include "RatesBounds.hcl"
						}
						else{
							#include "Rates.hcl"
						}
					#else
						#error Unknow boundary condition
					#endif
		                }
				i++;
			}			
		}
		//! 3.c.- Home cell, starting from head of chain.
		i = ihoc[lc];
		while(i < j) {
			if(!iIMove){
				#include "RatesSensors.hcl"
			}
			else{
				#if __BOUNDARY__==0
					// ElasticBounce boundary condition
					if(iIMove<0){
						#include "RatesBounds.hcl"
					}
					else{
						#include "Rates.hcl"
					}
				#elif __BOUNDARY__==1
					// Fix particles
					#include "Rates.hcl"
				#elif __BOUNDARY__==2
					// DeLeffe boundary condition
					if(iIMove<0){
						#include "RatesBounds.hcl"
					}
					else{
						#include "Rates.hcl"
					}
				#else
					#error Unknow boundary condition
				#endif
			}
			i++;
		}
	}
	//! 4th.- Append own particle as part of shepard term
	// Sensors not included
	if(iIMove){
		#if __BOUNDARY__==0 || __BOUNDARY__==2
			// Contour not included
			if(iIMove>0){
				#ifndef HAVE_3D
					conw = 1.f/(iHp*iHp);				// Different for 1d and 3d
				#else
					conw = 1.f/(iHp*iHp*iHp);			// Different for 1d and 3d
				#endif
				wab = kernelW(0.f)*conw*pmass[j];
				_SHEPARD_ += wab/iDens;
			}
		#else
			#ifndef HAVE_3D
				conw = 1.f/(iHp*iHp);				// Different for 1d and 3d
			#else
				conw = 1.f/(iHp*iHp*iHp);			// Different for 1d and 3d
			#endif
			wab = kernelW(0.f)*conw*pmass[j];
			_SHEPARD_ += wab/iDens;
		#endif
	}
	//! 5th.- Write output into global memory (at unsorted space)
	#ifndef __NO_LOCAL_MEM__
		sigma[labp]   =  _SIGMA_;
		f[labp]       =  _F_;
		drdt[labp]    =  _DRDT_ + _DRDT_F_;
		drdt_F[labp]  =  _DRDT_F_;
		shepard[labp] =  _SHEPARD_;
		gradW[labp]   =  _GRADW_;
	#else
		drdt[labp]   +=  _DRDT_F_;
	#endif

	// ---- A ---- Your code here ---- A ----
	// ---- | ------------------------ | ----

}