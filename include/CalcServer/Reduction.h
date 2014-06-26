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

#ifndef REDUCTION_H_INCLUDED
#define REDUCTION_H_INCLUDED

#include <deque>

#include <CalcServer.h>
#include <CalcServer/Kernel.h>

namespace Aqua{ namespace CalcServer{

/** @class Reduction Reduction.h CalcServer/Reduction.h
 * @brief Array reduction tool. It could every prefix scan operation that you
 * want.
 */
class Reduction : public Aqua::CalcServer::Tool
{
public:
	/** Constructor.
	 * @param name Tool name.
	 * @param input_name Variable to be reduced name.
	 * @param output_name Variable where the reduced value will be stored.
	 * @param operation The reduction operation. For instance:
	 *   - "a += b;"
	 *   - "a.x = (a.x < b.x) ? a.x : b.x; a.y = (a.y < b.y) ? a.y : b.y;"
	 * @param null_val The value considered as the null one, i.e. INFINITY for
	 * float min value reduction, or (vec2)(0.f,0.f) for a 2D vec sum reduction.
	 * @note Some helpers are available for null_val:
	 *   - VEC_ZERO: Zeroes vector.
	 *   - VEC_ONE: Ones vector, in 3D cases the last component will be zero.
	 *   - VEC_ALL_ONE: Equal to VEC_ONE, but in 3D cases the last component will be one as well.
	 *   - VEC_INFINITY: INFINITY components vector, in 3D cases the last component will be zero.
	 *   - VEC_ALL_INFINITY: Equal to VEC_INFINITY, but in 3D cases the last component will be INFINITY as well.
	 *   - VEC_NEG_INFINITY: -VEC_INFINITY
	 *   - VEC_ALL_NEG_INFINITY: -VEC_ALL_INFINITY.
	 */
	Reduction(const char *name,
              const char *input_name,
              const char *output_name,
              const char* operation,
              const char* null_val);

	/** Destructor.
	 */
	~Reduction();

    /** Initialize the tool.
     * @return false if all gone right, true otherwise.
     */
    bool setup();

	/** Compute the reduction.
	 * @return false if all gone right, true otherwise.
	 */
	bool execute();

    /** Number of steps needed to compute the reduction
     * @return Number of steps needed.
     */
    unsigned int nSteps(){return _global_work_sizes.size();}

private:
    /** Get the input and output variables
     * @return false if all gone right, true otherwise
     */
    bool variables();

	/** Setup the OpenCL stuff
	 * @return false if all gone right, true otherwise.
	 */
	bool setupOpenCL();

    /** Compile the source code and generate the corresponding kernel
     * @param source Source code to be compiled.
     * @param local_work_size Desired local work size.
     * @return Kernel instance, NULL if error happened.
     */
    cl_kernel compile(const char* source, size_t local_work_size);

    /** Update the input looking for changed value.
     * @return false if all gone right, true otherwise.
     */
    bool setVariables();

	/// Input variable name
	char* _input_name;
	/// Output variable name
	char* _output_name;
	/// Operation to be computed
	char* _operation;
	/// Considered null val
	char* _null_val;

    /// Input variable
    InputOutput::ArrayVariable *_input_var;
    /// Output variable
    InputOutput::Variable *_output_var;

    /// Set input variable
    cl_mem *_input;

	/// OpenCL kernels
	std::deque<cl_kernel> _kernels;

    /// Global work sizes in each step
    std::deque<size_t> _global_work_sizes;
    /// Local work sizes in each step
    std::deque<size_t> _local_work_sizes;
    /// Number of work groups in each step
    std::deque<size_t> _number_groups;
    /// Number of input elements for each step
    std::deque<size_t> _n;

    /// Memory objects
    std::deque<cl_mem> _mems;


};

}}  // namespace

#endif // REDUCTION_H_INCLUDED
