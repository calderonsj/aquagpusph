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

import numpy as np
import os.path as path
import aquagpusph as aqua


def main():
    # We should compute mgas y hg
    Ltot = aqua.get("Ltot")
    Htot = aqua.get("Htot")
    Vgas = Ltot * Htot - aqua.get("Vol")
    rhog = aqua.get("rhog")    
    mgas = Vgas / rhog
    hg = Vgas / Ltot
    aqua.set("mgas", mgas)
    aqua.set("hg", hg)
    # Compute the rate of variation of the gas temperature
    kg = aqua.get("kgas")
#    Tsolid = aqua.get("Tsolid")
    q = aqua.get("q")
    Tgas = aqua.get("Tgas")
    cvgas = aqua.get("cvgas")
#    dTdt = kg / (mgas * hg * cvgas) * (Tsolid - Tgas) #Para dif temperatura
    dTdt = q * Ltot / (cvgas * mgas)
    # Integrate the temperature
    dt = aqua.get("dt")
    Tgas += dTdt * dt
    aqua.set("Tgas", Tgas)
    # Compute the vapour pressure
    Ru = aqua.get("Ru")
    M = aqua.get("molargas")
    Rgas = Ru / M
    pgas = rhog * Rgas * Tgas
    aqua.set("pgas", pgas)
    # Compute the rate of variation of the saturation temperature
    kl = aqua.get("kliq")
    cvliq = aqua.get("cvliq")
    Tsat = aqua.get("Tsat")
    Mint = aqua.get("Mint")
    hint = aqua.get("hint")    
    dTsatdt = kl / (mgas * hint * cvliq) * (Tgas - Tsat)    
    # Compute the saturation temperature at the interface
    Tsat += dTsatdt *dt
    aqua.set("Tsat", Tsat)
    # Compute the saturation pressure (Clasius-Clapeyron)
    A = aqua.get("actant")
    B = aqua.get("bctant")    
    pSat = np.exp(( -A/Tsat ) + B)
    aqua.set("pSat", pSat)
    # Compute the mass transfer
    alphac = aqua.get("alphac")
    alphav = aqua.get("alphav")
    mCoeff = np.sqrt( M / (2 * np.pi * Ru ) )
    dmgasdt = Ltot * mCoeff * ( - alphac * (pgas / np.sqrt(Tgas) ) + alphav * (pSat / np.sqrt(Tsat)) )
#    dmgasdt = Ltot * mCoeff * (1.0/np.sqrt(Tsat)) * ( - alphac * pgas + alphav * pSat )  
    # Ltot is to be replaced by the actual length of the free surface
    # dmgasdt is given in kg/s       
    dmgas = dt*dmgasdt
    aqua.set("dmgas",dmgas)
    mgas += dt*dmgasdt
    aqua.set("mgas",mgas)
    Mint -= dt*dmgasdt
    aqua.set("Mint",Mint)
    # Update other variables
    rhog = mgas / Vgas
    aqua.set("rhog",rhog)
    return True