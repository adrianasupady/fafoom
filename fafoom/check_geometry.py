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

''' This module is responsible for checking the validity of the 3D geometry'''
import numpy as np
import StringIO
## Returns 'True' if the geometry is acceptable.

def is_geometry_valid(xyz_string, atoms, bonds, matrix_of_connections, distance_cutoff_1, distance_cutoff_2):
    '''is_geometry_valid(xyz_string, atoms, bonds, matrix_of_connections, distance_cutoff_1, distance_cutoff_2) -> bool
    
    This function checks, if there are no sterical clashes in the 3D geometry and if bonded pairs of atoms are not too far from each other.
    '''
    lines = StringIO.StringIO(xyz_string)
    array = np.loadtxt(lines, comments='#', delimiter=None, converters=None, skiprows=2, usecols=xrange(1,4), unpack=False )
    dist = np.zeros( (len(range(atoms)),len(range(atoms))))

    def distance(x,y):
        '''Calculates distance between two points in 3D '''
        return np.sqrt((x[0]-y[0])**2+(x[1]-y[1])**2+(x[2]-y[2])**2)
    
    for x in range(atoms):
        for y in range(atoms):
            a = np.array(array[float(x)])
            b = np.array(array[float(y)])
            dist[x][y] = distance(a,b)
    
    matrix = matrix_of_connections
    tocheck = dist + matrix # this ensures that distances between bonded atoms are "safe"  
    x,y = 0,0
    condition_clashes = True
    short_dist = []
    
    while True:    
        for x in range(atoms):
            for y in range(atoms):
                dis = tocheck[x][y]
                if (dis < distance_cutoff_1 and dis > 0): 
                    condition_clashes = False
                    check = False
                    return check
                else:
                    continue
        if condition_clashes == True:
            tocheck2 = np.multiply(matrix, dist)
            x,y = 0,0
            condition_bonded = True
            for x in range(atoms):
                for y in range(atoms):
                    dis2 = tocheck2[x,y]
                    if (dis2 > distance_cutoff_2): 
                        condition_bonded=False
                        break
                    else:
                        continue
                    
            if (condition_bonded == True) :
                check = True          
            else:
                check = False
        else:
            check = False
        return check