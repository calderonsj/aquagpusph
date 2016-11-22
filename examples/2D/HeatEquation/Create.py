#! /usr/bin/env python
#########################################################################
#                                                                       #
#            #    ##   #  #   #                           #             #
#           # #  #  #  #  #  # #                          #             #
#          ##### #  #  #  # #####  ##  ###  #  #  ## ###  ###           #
#          #   # #  #  #  # #   # #  # #  # #  # #   #  # #  #          #
#          #   # #  #  #  # #   # #  # #  # #  #   # #  # #  #          #
#          #   #  ## #  ##  #   #  ### ###   ### ##  ###  #  #          #
#                                    # #             #                  #
#                                  ##  #             #                  #
#                                                                       #
#########################################################################
#
#  This file is part of AQUA-gpusph, a free CFD program based on SPH.
#  Copyright (C) 2012  Jose Luis Cercos Pita <jl.cercos@upm.es>
#
#  AQUA-gpusph is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  AQUA-gpusph is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with AQUA-gpusph.  If not, see <http://www.gnu.org/licenses/>.
#
#########################################################################

import os.path as path
import math
import numpy as np

# Input data
# ==========
hfac = 2.0
g = 0.0
cs = 2.0
gamma = 1.0
refd = 998.0
# Tank dimensions
H = 1.0
L = 1.0
# Number of particles in x direction
nx = 100
# Dimensions and number of particles readjustment
# ===============================================

dr = L / nx
ny = int(round(H / dr))
n = nx * ny
print ("dr = {}".format(dr))
# Solid boundary elements
Nx = nx
Ny = ny + 2
print ("Nx = {}".format(Nx))
print ("Ny = {}".format(Ny))
print ("N = {}".format(Nx * Ny))
# Correct tank dimensions
L = Nx * dr
H = Ny * dr

DeLeffeDistFactor = 1
Nx = DeLeffeDistFactor * Nx
N = 2 * Nx
# Particles generation
# ====================
prb = cs * cs * refd / gamma

print("Opening output file...")
output = open("Fluid.dat", "w")
string = """#############################################################
#                                                           #
#    #    ##   #  #   #                           #         #
#   # #  #  #  #  #  # #                          #         #
#  ##### #  #  #  # #####  ##  ###  #  #  ## ###  ###       #
#  #   # #  #  #  # #   # #  # #  # #  # #   #  # #  #      #
#  #   # #  #  #  # #   # #  # #  # #  #   # #  # #  #      #
#  #   #  ## #  ##  #   #  ### ###   ### ##  ###  #  #      #
#                            # #             #              #
#                          ##  #             #              #
#                                                           #
#############################################################
"""
output.write(string)
print(string)

string = """
    Writing fluid particles...
"""
print(string)

Percentage = -1
for i in range(0, n):
    if Percentage != (i * 100) / n:
        Percentage = (i * 100) / n
        if not Percentage % 10:
            string = '    {}%'.format(Percentage)
            print(string)
    j = i
    idx = j % nx
    idy = j / nx
    imove = 1
    pos = ((0.5 + idx) * dr,
           (0.5 + idy) * dr)
    press = refd * g * (H - pos[1])
    dens = refd + press / cs**2 
    mass = dens * dr**2.0
    T = 0.0    #Set the fluid particles temperature
    dTdt = 0.0    
    if pos[1] > 0.47 and pos[1] < 0.53:
        T = 1.0
        dTdt = 1.0    
    string = ("{} {}, " * 4 + "{}, {}, {}, {}, {}, {}\n").format(
        pos[0], pos[1],
        0.0, 0.0,
        0.0, 0.0,
        0.0, 0.0,
        dens,
        0.0,
        mass,
        T,
        dTdt, # dT/dt
        imove)
    output.write(string)
print('    100%')

#string = """
#    Writing the boundary elements...
#"""
#print(string)
#Percentage = -1
#for i in range(0, N):
#    if Percentage != (i * 100) / N:
#        Percentage = (i * 100) / N
#        if not Percentage % 10:
#            string = '    {}%'.format(Percentage)
#            print(string)
    # Bottom
#    if (i < Nx):
#        j = i
#        idx = j 
#        idy = 0.0
#        normal = [0.0, -1.0]
    # Roof    
#    elif(i < 2 * Nx):
#        j = i - Nx
#        idx = j
#        idy = Ny - 2
#        normal = [0.0, 1.0]
#    pos =  ((0.5 + idx) * dr / DeLeffeDistFactor,
#           idy * dr / DeLeffeDistFactor)
#    dens = refd
#    press = prb*pow(dens / refd, gamma - 1.0)    
#    imove = -3
#    mass = dr / DeLeffeDistFactor
#    T = 0.0
#    dTdt = 0.0
#    string = ("{} {}, " * 4 + "{}, {}, {}, {}, {}, {}\n").format(
#        pos[0], pos[1],
#        0.0, -1.0,
#        0.0, 0.0,
#        0.0, 0.0,
#        dens,
#        0.0,
#        mass,
#        T,
#        dTdt, # dT/dt
#        imove)
#    output.write(string)
#print('    100%')
output.close()
