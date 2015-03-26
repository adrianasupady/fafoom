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
''' Collection of diverse help/convert functions '''
import os
import numpy as np
import math
import shutil

# Flow-handling


def backup(filename, obj):
    """ Write the representation of an object (or objects) to a file."""
    with open(filename, 'w') as outf:
        if hasattr(obj, "__len__"):
            for i in range(len(obj)):
                outf.write("%s\n" % repr(obj[i]))
        else:
            outf.write("%s\n" % repr(obj))


def boolean(string):
    """Recover the boolean value from a string and return it."""
    if string in ["False", "false", "FALSE"]:
        return False
    if string in ["True", "true", "TRUE"]:
        return True
    raise ValueError("Cannot be converted to a boolean type")


def print_output(text):
    """Write text to the 'output.txt'. Create if not existing."""
    if os.path.isfile("output.txt"):
        f = open("output.txt", "a")
        f.write(str(text)+'\n')
        f.close()
    else:
        f = open("output.txt", "w")
        f.write(str(text)+'\n')
        f.close()


def remover_file(instance):
    """Remove a file (if it exists)."""
    try:
        os.remove(instance)
    except OSError:
        pass


def remover_dir(instance):
    """Remove a directory (if it exists)."""
    try:
        shutil.rmtree(instance)
    except OSError:
        pass


def file2string(input_file):
    """Read file to a string and return it."""
    with open(input_file, 'r') as f:
        string = f.read()
    f.close()
    return string


def string2file(string, filename):
    """Write string to a file"""
    with open(filename, 'w') as target:
        target.write(string)
    target.close()

# Help vector/matrix functions


def get_vec(vec1, vec2):
    """Calculate difference between vectors of angles.
    Args:
        vec1 (list)
        vec2 (list)
    Returns:
        numpy array
    Raises:

        ValueError: if the length of the lists differ

    Warning1: the vectors contain periodic values, i.e -185 -> 175
    Warning2: symmetry is not included here, but can be easily added if the
    index of 'symmetric' torsion is known
    """
    if len(vec1) != len(vec2):
        raise ValueError("No length match between the lists")
    diff_vec = np.zeros((1, len(vec1)))
    for i in range(1, len(vec1)+1):
        tor_diff = abs(vec1[i-1]-vec2[i-1])
        diff_vec[0][i-1] = min(abs(tor_diff), abs(360-tor_diff))/180.0
    return diff_vec[0]


def tor_rmsd(p, vec):
    """Calculate the modified p norm.The difference from standard norm is the
    fact that the addends are divided by the length of the vector."""
    summe = 0
    for i in range(0, len(vec)):
        summe += math.pow(abs(vec[i]), p)
    return math.pow(summe/len(vec), (1.0/p))


def find_two_in_list(list_sum, list_to_search):
    """A list is mapped to a segment of a line which length is equal to 1.
    The lengths of the segments are proportional to the corresponding list
    values. Next, two random numbers between 0 and 1 are generated and the
    segments containing these random numbers are returned."""
    rn1 = list_sum*np.random.rand()
    found1 = False
    index = 1
    while not found1:
        if rn1 < list_to_search[:index].sum(axis=0):
            found1 = index
        else:
            index += 1
    equal = True
    while equal:
        rn2 = list_sum*np.random.rand()
        found2 = False
        index = 1
        while not found2:
            if rn2 < list_to_search[:index].sum(axis=0):
                found2 = index
            else:
                index += 1
        if found2 != found1:
            equal = False
    return found1, found2


def check_geo_sdf(sdf_string, cutoff1, cutoff2):
    """Check geometry from a sdf_string for clashes.

    Args:
        sdf_string (str)
        distance_cutoff_1 (float): min distance between non-bonded atoms [A]
        distance_cutoff_2 (float): max distance between bonded atoms [A]
    Returns:
        True for clash-free geometries and False for invalid geometries
    Raises:
        ValueError: if distance cutoffs are non-positive
    """
    if cutoff1 <= 0 or cutoff2 <= 0:
        raise ValueError("Distance cutoff needs to be a positive float")

    def distance(x, y):
        """"Calculate distance between two points in 3D."""
        return np.sqrt((x[0]-y[0])**2+(x[1]-y[1])**2+(x[2]-y[2])**2)

    atoms = int(sdf_string.split('\n')[3].split()[0])
    bonds = int(sdf_string.split('\n')[3].split()[1])
    coordinates = np.zeros((atoms, 3))
    bonds_pairs = np.zeros((bonds, 2))
    for i in range(4, atoms+4):
        coordinates[i-4][0] = sdf_string.split('\n')[i].split()[0]
        coordinates[i-4][1] = sdf_string.split('\n')[i].split()[1]
        coordinates[i-4][2] = sdf_string.split('\n')[i].split()[2]
    for i in range(atoms+4, atoms+bonds+4):
        bonds_pairs[i-atoms-4][0] = sdf_string.split('\n')[i].split()[0]
        bonds_pairs[i-atoms-4][1] = sdf_string.split('\n')[i].split()[1]

    dist = np.zeros((atoms, atoms))
    for x in range(atoms):
        for y in xrange(x, atoms):
            a = np.array(coordinates[float(x)])
            b = np.array(coordinates[float(y)])
            dist[x][y] = distance(a, b)
            dist[y][x] = dist[x][y]
    list_of_bonds = []
    for i in range(bonds):
        list_of_bonds.append([bonds_pairs[i][0], bonds_pairs[i][1]])

    def check_distance():
        for x in range(atoms):
            for y in xrange(x+1, atoms):
                if [x+1, y+1] not in list_of_bonds and \
                   [y+1, x+1] not in list_of_bonds and dist[x][y] < cutoff1:
                    check = False
                    return check
                elif ([x+1, y+1] in list_of_bonds and dist[x][y] > cutoff2) or\
                     ([y+1, x+1] in list_of_bonds and dist[x][y] > cutoff2):
                    check = False
                    return check
        return True

    return check_distance()

# Format conversions


def sdf2aims(sdf_string):
    """Convert a sdf string to a aims string."""
    atoms = int(sdf_string.split('\n')[3].split()[0])
    coord = []
    for i in range(4, 4+atoms):
        x = float(sdf_string.split('\n')[i].split()[0])
        y = float(sdf_string.split('\n')[i].split()[1])
        z = float(sdf_string.split('\n')[i].split()[2])
        name = sdf_string.split('\n')[i].split()[3]
        coord.append('%s%10.4f%10.4f%10.4f%2s' % ('atom', x, y, z, name))
        if i != 3+atoms:
            coord.append('\n')
    aims_string = ''.join(coord)
    return aims_string


def sdf2xyz(sdf_string):
    """Convert a sdf string to a xyz string."""
    atoms = int(sdf_string.split('\n')[3].split()[0])
    coord = [str(atoms)+('\n')]
    for i in range(4, 4+atoms):
        x = float(sdf_string.split('\n')[i].split()[0])
        y = float(sdf_string.split('\n')[i].split()[1])
        z = float(sdf_string.split('\n')[i].split()[2])
        name = sdf_string.split('\n')[i].split()[3]
        coord.append('\n%2s%10.4f%10.4f%10.4f' % (name, x, y, z))
    xyz_string = ''.join(coord)
    return xyz_string


def aims2sdf(aims_string, sdf_template_string):
    """Convert a aims string to a sdf string. Template for the sdf string
    required."""
    atoms = len(aims_string.splitlines())
    sdf_form = sdf_template_string.splitlines()
    c = []
    cnt = 0
    for i in range(len(sdf_form)):
        if i > 3 and i < 4+atoms:
            line = sdf_form[i].split()
            line[0] = aims_string.split()[5*cnt+1]
            line[1] = aims_string.split()[5*cnt+2]
            line[2] = aims_string.split()[5*cnt+3]
            cnt += 1
            c.append('%10.4f%10.4f%10.4f%2s' % (float(line[0]),
                                                float(line[1]),
                                                float(line[2]), line[3]))
            for j in xrange(4, len(line)):
                if j == 4:
                    c.append('%4d' % int(line[j]))
                elif j == len(line)-1:
                    c.append('%3d\n' % int(line[j]))
                else:
                    c.append('%3d' % int(line[j]))
        else:
            if i != len(sdf_form)-1:
                c.append(''.join(sdf_form[i])+'\n')
            else:
                c.append(''.join(sdf_form[i]))
    sdf_string = ''.join(c)
    return sdf_string


def xyz2sdf(xyz_string, sdf_template_string):
    """Convert a xyz string to a sdf string. Template for the sdf string
    required."""
    arr = xyz_string.splitlines()
    atoms = int(arr[0].split()[0])
    xyz_string_cut = '\n'.join(arr[2:])
    sdf_form = sdf_template_string.splitlines()
    c = []
    cnt = 0
    for i in range(len(sdf_form)):
        if i > 3 and i < 4+atoms:
            line = sdf_form[i].split()
            line[0] = xyz_string_cut.split()[4*cnt+1]
            line[1] = xyz_string_cut.split()[4*cnt+2]
            line[2] = xyz_string_cut.split()[4*cnt+3]
            cnt += 1
            c.append('%10.4f%10.4f%10.4f%2s' % (float(line[0]),
                                                float(line[1]),
                                                float(line[2]), line[3]))
            for j in xrange(4, len(line)):
                if j == 4:
                    c.append('%4d' % int(line[j]))
                elif j == len(line)-1:
                    c.append('%3d\n' % int(line[j]))
                else:
                    c.append('%3d' % int(line[j]))
        else:
            if i != len(sdf_form)-1:
                c.append(''.join(sdf_form[i])+'\n')
            else:
                c.append(''.join(sdf_form[i]))
    sdf_string = ''.join(c)
    return sdf_string


def mirror_sdf(sdf_string):
    """Mirror the geometry from a sdf string. Return a new sdf string."""
    atoms = int(sdf_string.split('\n')[3].split()[0])
    sdf_form = sdf_string.splitlines()
    c = []
    cnt = 0
    for i in range(len(sdf_form)):
        if i > 3 and i < 4+atoms:
            line = sdf_form[i].split()
            line[0] = -1.0*float(line[0])
            line[1] = -1.0*float(line[1])
            line[2] = -1.0*float(line[2])
            cnt += 1
            c.append('%10.4f%10.4f%10.4f%2s' % (float(line[0]),
                                                float(line[1]),
                                                float(line[2]), line[3]))
            for j in xrange(4, len(line)):
                if j == 4:
                    c.append('%4d' % int(line[j]))
                elif j == len(line)-1:
                    c.append('%3d\n' % int(line[j]))
                else:
                    c.append('%3d' % int(line[j]))
        else:
            if i != len(sdf_form)-1:
                c.append(''.join(sdf_form[i])+'\n')
            else:
                c.append(''.join(sdf_form[i]))
    mirror_sdf_string = ''.join(c)
    return mirror_sdf_string
