#    Copyright 2015 Adriana Supady 
#    adriana.supady@gmail.com
#
#    This file is part of fafoom.
#
#   Fafoom is free software: you can redistribute it and/or modify
#   it under the terms of the GNU Lesser General Public License as published by
#   the Free Software Foundation, either version 3 of the License, or
#   (at your option) any later version.
#
#   Fafoom is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU Lesser General Public License for more details.
#
#  You should have received a copy of the GNU Lesser General Public License
#   along with fafoom.  If not, see <http://www.gnu.org/licenses/>.
''' This module measures and set the torsions values'''
import sys
import openbabel
from operator import itemgetter

from utilities import sdf_to_sdf_readstring, sdf_to_xyz_writestring, xyz_to_xyz_readstring

def angle_measure(xyz_string, cistrans, torsion):
    ''' angle_measure(xyz_string, cistrans, torsion) -> vector1, vector2'''
    molecule = openbabel.OBMol()
    xyz_to_xyz_readstring(molecule, xyz_string)
    measured_torsion = []
    measured_cistrans = []
    if len(cistrans ) != 0: 
        for x in range(len(cistrans)):
            tmp = float(molecule.GetTorsion(itemgetter(0)(cistrans[x]),itemgetter(1)(cistrans[x]),itemgetter(2)(cistrans[x]),itemgetter(3)(cistrans[x])))
            measured_cistrans.append(float('{0:.2f}'.format(tmp)))
    if len(torsion) !=0:
        for y in range(len(torsion)):
            tmp = float(molecule.GetTorsion(itemgetter(0)(torsion[y]),itemgetter(1)(torsion[y]),itemgetter(2)(torsion[y]),itemgetter(3)(torsion[y])))
            measured_torsion.append(float('{0:.2f}'.format(tmp)))
    if len(cistrans) == 0 and len(torsion) == 0:
        sys.exit("There are no rotatable bonds and no cis/trans bonds")

    return(measured_cistrans, measured_torsion)           

    
    
def angle_set(sdf_string, cistrans, torsion, values_cistrans_n, values_tor_n):
    '''angle_set(sdf_string, cistrans, torsion, values_cistrans_n, values_tor_n) -> string'''
    molecule = openbabel.OBMol()
    sdf_to_sdf_readstring(molecule, sdf_string)
    deg_to_rad = 0.0174532925
    if len(cistrans) != 0 and len(torsion) != 0 : # creates the candidate molecule with previously chosen values 
        for x in range(len(cistrans)):
            for y in range(len(torsion)):
                molecule.SetTorsion(itemgetter(0)(torsion[y]),itemgetter(1)(torsion[y]),itemgetter(2)(torsion[y]),itemgetter(3)(torsion[y]),itemgetter(y)(values_tor_n)*deg_to_rad)
                molecule.SetTorsion(itemgetter(0)(cistrans[x]),itemgetter(1)(cistrans[x]),itemgetter(2)(cistrans[x]),itemgetter(3)(cistrans[x]),itemgetter(x)(values_cistrans_n)*deg_to_rad)
    elif len(cistrans) == 0 and len(torsion) !=0: # only torsions
        for y in range(len(torsion)):
            molecule.SetTorsion(itemgetter(0)(torsion[y]),itemgetter(1)(torsion[y]),itemgetter(2)(torsion[y]),itemgetter(3)(torsion[y]),itemgetter(y)(values_tor_n)*deg_to_rad)
    elif len(cistrans) != 0 and len(torsion) == 0: # only cis/trans
        for x in range(len(cistrans)):    
            molecule.SetTorsion(itemgetter(0)(cistrans[x]),itemgetter(1)(cistrans[x]),itemgetter(2)(cistrans[x]),itemgetter(3)(cistrans[x]),itemgetter(x)(values_cistrans_n)*deg_to_rad)
    else:                 # program closes if there are no degrees of freedom
        sys.exit("There are no rotatable bonds and no cis/trans bonds")
    
    return sdf_to_xyz_writestring(molecule)
