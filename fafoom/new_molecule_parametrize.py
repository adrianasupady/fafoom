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

''' This module creates the molecules from the smile code'''

import numpy as np
import StringIO
import pybel

from copy import copy
from operator import itemgetter

from utilities import print_output
from check_geometry import is_geometry_valid

def parametrize(smile, custom, smart_torsion, smart_cistrans, smart_custom):
    '''parametrize(smile, custom, smart_torsion, smart_cistrans, smart_custom) -> int, int, vector, vector, vector 
    
    This function outputs the torsion and cistrans positions as weel as the number of atoms and number of bonds
    '''
    mol = pybel.readstring("smi", smile)
    torsion = sorted(list(smart_torsion.findall(mol)), key=itemgetter(1))  
    cistrans = sorted(list(smart_cistrans.findall(mol)), key=itemgetter(1))
    custom_torsion = []

    def ig(x):
        return itemgetter(x)

    def cleaner(list_to_clean):
        for_remove = []
        for x in reversed(range(len(list_to_clean))):
            for y in reversed(range(x)):
                if (ig(1)(list_to_clean[x]) == ig(1)(list_to_clean[y]) and ig(2)(list_to_clean[x]) == ig(2)(list_to_clean[y])) or (ig(1)(list_to_clean[x]) == ig(2)(list_to_clean[y]) and ig(2)(list_to_clean[x]) == ig(1)(list_to_clean[y])):
                 for_remove.append(y)
        clean_list = [v for i, v in enumerate(list_to_clean) if i not in set(for_remove)]
        return clean_list

    if custom:
        custom = sorted(list(smart_custom.findall(mol)), key=itemgetter(1))
        to_del_bef_custom = []

        for x in reversed(range(len(torsion))):
            for y in reversed(range(len(custom))):
                if (itemgetter(1)(torsion[x]) == itemgetter(1)(custom[y]) and itemgetter(2)(torsion[x]) == itemgetter(2)(custom[y])) or (itemgetter(1)(torsion[x]) == itemgetter(2)(custom[y]) and itemgetter(2)(torsion[x]) == itemgetter(1)(custom[y])):
                    to_del_bef_custom.append(x)

        custom_torsion = copy(torsion)
        custom_torsion = [v for i, v in enumerate(custom_torsion) if i not in set(to_del_bef_custom)]
        custom_torsion = cleaner(custom_torsion)

    cistrans = cleaner(cistrans)
    torsion = cleaner(torsion)
    mol.make3D()    # creates a 3D representation applying force field!
    atoms = mol.OBMol.NumAtoms()
    bonds = mol.OBMol.NumBonds()
    
    return (atoms, bonds, torsion, cistrans, custom_torsion)

def template_sdf(smile, atoms, bonds, distance_cutoff_1, distance_cutoff_2):
    '''template_sdf(smile, atoms, bonds, distance_cutoff_1, distance_cutoff_2) -> array, string
    
    This function generates a template sdf file together with a matrix of connections that contains the information about atoms involved in bonds.
    '''
    cnt = 0
    sdf_check = True    
    while sdf_check:
        mol = pybel.readstring("smi", smile)
        mol.make3D()    # creates a 3D representation applying force field!
        mol.write("sdf", "mol.sdf", overwrite=True) # creates the sdf of the molecule of interest
        xyz_string = mol.write(format='xyz')
        ct_string = mol.write(format='ct')
        sdf_string = mol.write(format='sdf')
        lines = StringIO.StringIO(ct_string).readlines()
        con = np.array([[i for i in line.strip().split()] for line in lines[(atoms + 2):]])  # reads only the lines with indicies of connected atoms 
        for i in xrange(len(con)):
            if len(con[i][0]) == 5: # only if the indicies 'glued' in the sdf file
                a = list(con[i][0])
                b = con[i][1]
                c = con[i][2]
                con[i][0]=str(a[0])+str(a[1])
                con[i][1]=str(a[2])+str(a[3])+str(a[4])
                con[i][2]=b
                con[i].append(c)
                
            elif len(con[i][0]) == 6 : # only if the indicies 'glued' in the sdf file
                a = list(con[i][0])
                b = con[i][1]
                c = con[i][2]
                con[i][0]=str(a[0])+str(a[1])+str(a[2])
                con[i][1]=str(a[3])+str(a[4])+str(a[5])
                con[i][2]=b
                con[i].append(c)
    
        con1 = []
        con2 = []
        for i in xrange(len(con)):
            con1.append(con[i][0])
            con2.append(con[i][1])
        
        matrix_of_connections = np.zeros((atoms, atoms))
        x = 0
        while x < bonds:        # the atoms that are connected get 1 instead of 0
            k = int(con1[x]) - 1
            v = int(con2[x]) - 1
            matrix_of_connections[k][v] = 1
            matrix_of_connections[v][k] = 1
            x = x + 1
        
        check = is_geometry_valid(xyz_string, atoms, bonds, matrix_of_connections, distance_cutoff_1, distance_cutoff_2)
        if check == True: sdf_check = False
        else: cnt+=1
    return (matrix_of_connections, sdf_string)

