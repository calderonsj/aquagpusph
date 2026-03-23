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

import os
import sys
script_folder = os.path.dirname(os.path.realpath(__file__))
sys.path.append(os.path.join(script_folder, "../../"))
import aqua_example_utils as utils
import math


g = 9.81
# Tank dimensions
d = 1.0
D = 2.0 * d
L = 5.367 * d
# Fluid
H = d
B = 2. * d
# Discretization
ny = 100
dr = H / ny
nx = int(round(B / dr))
H = ny * dr
Nx = int(round(L / dr))
L = Nx * dr
Ny = int(round(D / dr))
D = Ny * dr
# SPH
hfac = 2.0
h = hfac * dr
courant = 0.2
Ma = 0.1
cs = (g * H)**0.5 / Ma
refd = 998.0
alpha = 0.0
visc_dyn = 8.9e-4
visc_dyn = max(alpha * refd * hfac * dr * cs / 8.0, visc_dyn)
delta = 0.0
# Time
T = 7.15 / (g / H)**0.5
FPS = 100

LLes = 2. * hfac * dr

Uref = Ma * cs

e0 = 0.0
def particle(f, r):
    global e0
    m = dr**2 * refd
    e0 += m * g * r[1]
    # Hydrostatic
    p = refd * g * (H - r[1])
    rho = refd + p / cs**2
    f.write("{}, {}, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, {}, 0.0, {}, 1\n".format(
        r[0], r[1], rho, m))

def boundary(f, r, n):
    press = max(0.0, refd * g * (H - r[1]))
    rho = refd + press / cs**2
    f.write("{}, {}, {}, {}, 0.0, 0.0, 0.0, 0.0, {}, 0.0, {}, -3\n".format(
        r[0], r[1], n[0], n[1], rho, dr))


def sensor(f, r):
    f.write("{}, {}, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, {}, 0.0, {}, 0\n".format(
        r[0], r[1], refd, dr))


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
n = nx * ny + 2 * (Nx + Ny)
with open("Fluid.dat", "w") as f:
    f.write(string)
    # Fluid particles
    x = 0.5 * dr
    while x < B:
        y = 0.5 * dr
        while y < H:
            particle(f, (x, y))
            y += dr
        x += dr
    # Bottom wall
    x = 0.5 * dr
    while x < L:
        boundary(f, (x, 0), (0, -1))
        x += dr
    # Top wall
    x = 0.5 * dr
    while x < L:
        boundary(f, (x, D), (0, 1))
        x += dr
    # Left wall
    y = 0.5 * dr
    while y < D:
        boundary(f, (0.0, y), (-1, 0))
        y += dr
    # Right wall
    y = 0.5 * dr
    while y < D:
        boundary(f, (L, y), (1, 0))
        y += dr

with open("Sensors.dat", "w") as f:
    sensor(f, (L + 0.5 * dr, 0.003))
    sensor(f, (L + 0.5 * dr, 0.015))
    sensor(f, (L + 0.5 * dr, 0.030))
    sensor(f, (L + 0.5 * dr, 0.080))

domain_min = (-1.5 * L, -0.5 * D)
domain_min = str(domain_min).replace('(', '').replace(')', '')
domain_max = (1.5 * L, 1.5 * D)
domain_max = str(domain_max).replace('(', '').replace(')', '')

data = {'DR':str(dr), 'HFAC':str(hfac), 'CS':str(cs), 'COURANT':str(courant),
        'G':str(g), 'L':str(L), 'D':str(D), 'B':str(H), 'H':str(H),
	'DOMAIN_MIN':domain_min, 'DOMAIN_MAX':domain_max,
        'REFD':str(refd), 'T':str(T), 'FPS':str(FPS), 'N':str(n),
        'VISC_DYN':str(visc_dyn), 'E0':str(e0), 'N_SENSORS':str(4),
	'DELTA':str(delta), 'UREF':str(Uref), 'LLES':str(LLes)}
utils.configure(data, os.path.join(script_folder, "templates"))
