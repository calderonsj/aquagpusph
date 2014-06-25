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

#ifndef CALCSERVER_H_INCLUDED
#define CALCSERVER_H_INCLUDED

#include <CL/cl.h>

#include <deque>

#include <sphPrerequisites.h>
#include <Variable.h>
#include <Singleton.h>
#include <CalcServer/Tool.h>

#ifndef _ITEMS
	#define _ITEMS  128
#endif

#ifndef _GROUPS
	#define _GROUPS 32
#endif


namespace Aqua{
/** @namespace CalcServer Calculation server name space.
 */
namespace CalcServer{

/** @class CalcServer CalcServer.h CalcServer.h
 * @brief Entity that perform the main work of the simulation. Therefore this
 * class have an internal loop which iterates while no data must be sent to
 * the host for an output file writing process.
 * @remarks Some output files are managed internally by this class, like the
 * log file, the energy file, pressure sensors file, or bounds file.
 */
class CalcServer : public Aqua::Singleton<Aqua::CalcServer::CalcServer>
{
public:
	/** Constructor
	 */
	CalcServer();

	/** Destructor
	 */
	~CalcServer();

	/** Internal server loop. Calculation server will be iterating while the
	 * next output file event is reached (or the simulation is finished).
	 * @return false if all gone right, true otherwise.
	 */
	bool update();

	/** Setup the calculation server with the data recopilated by the host
	 * during the initialization process.
	 * @return false if the caluclation server has been succesfully setup,
     * true otherwise
	 */
	bool setup();

    /** Get the variables manager
     * @return Variables manager
     */
    InputOutput::Variables* variables() const {return _vars;}

    /** Get the active context
     * @return OpenCL context
     */
    cl_context context() const{return _context;}

    /** Get the platform
     * @return OpenCL platform
     */
    cl_platform_id platform() const{return _platform;}

    /** Get the device
     * @return OpenCL device
     */
    cl_device_id device() const{return _device;}

    /** Get the command queue
     * @return OpenCL command queue
     */
    cl_command_queue command_queue() const{return _command_queue;}
private:
	/** Setup the OpenCL stuff.
	 * @return false if the OpenCL environment has been succesfully built,
	 * true otherwise
	 */
	bool setupOpenCL();
	/** Prints all the available platforms and devices returned by OpenCL.
	 * @return false if the OpenCL environment can be succesfully built,
	 * true otherwise
	 */
	bool queryOpenCL();
	/** Get a platform from the available ones.
	 * @return false if a platform could be obtained, true otherwise
	 */
	bool setupPlatform();
	/** Get the available devices in the selected platform.
	 * @return false if the devices have been succesfully obtained, true
	 * otherwise
	 */
	bool setupDevices();

	/// Number of available platforms
	cl_uint _num_platforms;
	/// Array of platforms
	cl_platform_id *_platforms;
	/// Number of devices
	cl_uint _num_devices;
	/// Array of devices
	cl_device_id *_devices;
	/// OpenCL context
	cl_context _context;
	/// OpenCL command queue
	cl_command_queue *_command_queues;
	/// Selected platform
	cl_platform_id _platform;
	/// Selected device
	cl_device_id _device;
	/// Selected command queue
	cl_command_queue _command_queue;

    /// User registered variables
    InputOutput::Variables *_vars;

    /// User registered tools
    std::deque<Tool*> _tools;
};

}}  // namespace

#endif // CALCSERVER_H_INCLUDED
