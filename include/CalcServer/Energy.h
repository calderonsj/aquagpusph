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

#ifndef ENERGY_H_INCLUDED
#define ENERGY_H_INCLUDED

#include <CalcServer/Kernel.h>
#include <CalcServer/Reduction.h>

namespace Aqua{ namespace CalcServer{

/** @class Energy Energy.h CalcServer/Energy.h
 * @brief Computes the fluid energy components:
 *   -# Potential energy: \f$ E_{pot} = - \sum_i m_i
     \mathbf{g} \cdot \mathbf{r}_i \f$.
 *   -# Kinetic energy: \f$ E_{kin} = \sum_i \frac{1}{2} m_i
     \vert \mathbf{u}_i \vert^2 \f$.
 *   -# Internal energy: \f$ U = \int_0^t \sum_i \frac{p_i}{\rho_i^2}
     \left(
        \frac{\mathrm{d} \rho_i}{\mathrm{d} t}
        - \left. \frac{\mathrm{d} \rho_i}{\mathrm{d} t} \right\vert_F
     \right) m_i \mathrm{d}t \f$.
 *   -# Enthalpy: \f$ H = \int_0^t \sum_i \frac{p_i}{\rho_i^2}
     \frac{\mathrm{d} \rho_i}{\mathrm{d} t} m_i \mathrm{d}t \f$.
 *   -# Entropy: \f$ TS = U - H \f$.
 *   -# Total energy: \f$ E = U + E_{kin} \f$.
 * @remarks Since some energy components must be integrated in time, a low
 * energy file output/update frequency may imply too big time steps for the
 * numerical integration process whith poor results.
 * @note The energy computation is accelerated with OpenCL, so its computation
 * should not significantly affect to the preformance.
 */
class Energy : public Aqua::CalcServer::Kernel
{
public:
	/** Constructor.
	 */
	Energy();

	/** Destructor.
	 */
	~Energy();

	/** Get the resultant total energy.
	 * @return Total energy: \f$ E = U + E_{kin} \f$.
	 * @warning The viscous dissipation is not implemented yet.
	 */
	float energy(){return _E.x + _E.w;}

	/** Get the internal energy.
	 * @return Internal energy: \f$ U = \int_0^t \sum_i \frac{p_i}{\rho_i^2}
	     \left(
            \frac{\mathrm{d} \rho_i}{\mathrm{d} t}
            - \left. \frac{\mathrm{d} \rho_i}{\mathrm{d} t} \right\vert_F
	     \right) m_i \mathrm{d}t \f$.
	 * @warning The viscous dissipation is not implemented yet.
	 */
	float internalEnergy(){return _E.x;}

	/** Get the enthalpy.
	 * @return Enthalpy: \f$ H = \int_0^t \sum_i \frac{p_i}{\rho_i^2}
         \frac{\mathrm{d} \rho_i}{\mathrm{d} t} m_i \mathrm{d}t \f$.
	 * @warning The viscous dissipation is not implemented yet.
	 */
	float enthalpy(){return _E.y;}

	/** Get the entropy.
	 * @return Entropy: \f$ TS = U - H \f$.
	 * @warning The viscous dissipation is not implemented yet.
	 */
	float entropy(){return _E.x - _E.y;}

	/** Get the potential energy.
	 * @return Potential energy: \f$ E_{pot} = - \sum_i m_i
	     \mathbf{g} \cdot \mathbf{r}_i \f$.
	 */
	float potentialEnergy(){return _E.z;}

	/** Get the total kinetic energy.
	 * @return Kinetic energy: \f$ E_{kin} = \sum_i \frac{1}{2} m_i
	     \vert \mathbf{u}_i \vert^2 \f$.
	 */
	float kineticEnergy(){return _E.w;}

	/** Compute the energy.
	 * @return false if all gone right, true otherwise.
	 */
	bool execute();

private:
	/** Setup the energy OpenCL stuff.
	 * @return false if all gone right, true otherwise.
	 */
	bool setupEnergy();

	/** Setup the reduction tool
	 * @return false if all gone right, true otherwise.
	 */
	bool setupReduction();

	/// Server allocated energy array.
	cl_mem _device_energy;
	/// Host allocated energy rate of change.
	vec4 _dEdt;
	/// Host allocated energy.
	vec4 _E;
	/// Last time when the energy was computed
	float _time;
	/// Kernel path
	char *_path;
	/// OpenCL program
	cl_program _program;
	/// OpenCL kernel
	cl_kernel _kernel;
	/// Global work size
	size_t _global_work_size;
	/// Local work size
	size_t _local_work_size;
    /// Energy values reduction tool
    Reduction *_reduction;
};

}}  // namespace

#endif // ENERGY_H_INCLUDED