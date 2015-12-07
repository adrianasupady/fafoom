#    Copyright 2015 Adriana Supady
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
import ConfigParser

from rdkit import Chem
from rdkit.Chem import AllChem

from operator import itemgetter

# Flow-handling


def backup(filename, obj):
    """ Write the representation of an object (or objects) to a file."""
    with open(filename, 'w') as outf:
        if hasattr(obj, "__len__"):
            for i in range(len(obj)):
                outf.write("%s\n" % repr(obj[i]))
        else:
            outf.write("%s\n" % repr(obj))
    outf.close()


def boolean(string):
    """Recover the boolean value from a string and return it."""
    if string in ["False", "false", "FALSE"]:
        return False
    if string in ["True", "true", "TRUE"]:
        return True
    raise ValueError("Cannot be converted to a boolean type")


def number(s):
    """Convert to integer of float if needed"""
    try:
        return int(s)
    except ValueError:
        return float(s)


def print_output(text):
    """Write text to the 'output.txt'. Create it if needed."""
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
    """Read a file to a string and return it."""
    with open(input_file, 'r') as f:
        string = f.read()
    f.close()
    return string


def string2file(string, filename):
    """Write a string to a file"""
    with open(filename, 'w') as target:
        target.write(string)
    target.close()


def set_default(params, dict_default):
    """Set defaults for missing keys and add the key:value pairs to the
    dict."""
    for key in dict_default:
        if key not in params:
            print_output("Setting a default value for "+str(key)+": " +
                         str(dict_default[key]))
            params[str(key)] = dict_default[key]
    return params


def file2dict(filename, sections):
    """Parse a file and create a dictionary"""
    config = ConfigParser.RawConfigParser()
    config.read(filename)
    new_dict = {}
    for section in sections:
        if config.has_section(section):
            for key, value in config.items(section):
                new_dict[str(key)] = eval(value)
    return new_dict

# Help vector/matrix functions


def ig(x):
    return itemgetter(x)


def cleaner(list_to_clean):
    """ Remove duplicate torsion definion from a list of atom ind. tuples."""
    for_remove = []
    for x in reversed(range(len(list_to_clean))):
        for y in reversed(range(x)):
            ix1, ix2 = ig(1)(list_to_clean[x]), ig(2)(list_to_clean[x])
            iy1, iy2 = ig(1)(list_to_clean[y]), ig(2)(list_to_clean[y])
            if (ix1 == iy1 and ix2 == iy2) or (ix1 == iy2 and ix2 == iy1):
                for_remove.append(y)
    clean_list = [v for i, v in enumerate(list_to_clean)
                  if i not in set(for_remove)]
    return clean_list


def get_vec(vec1, vec2):
    """Calculate difference between vectors of angles [in rad!].
    Args:
        vec1 (list) in deg
        vec2 (list) in deg
    Returns:
        numpy array
    Raises:

        ValueError: if the length of the lists differ

    Warning1: the vectors contain periodic values, i.e -185 -> 175
    Warning2: symmetry is not included here, but can be easily added if the
    index of the 'symmetric' torsion is known
    """
    if len(vec1) != len(vec2):
        raise ValueError("No length match between the lists")
    diff_vec = np.zeros(len(vec1))
    for i in range(0, len(vec1)):
        tor_diff = abs(vec1[i]-vec2[i])
        diff_vec[i] = min(abs(tor_diff), abs(360-tor_diff))/180.0
    return diff_vec


def tor_rmsd(p, vec):
    """Calculate the modified p norm.The difference from standard norm is the
    fact that the addends are divided by the length of the vector."""
    summe = 0
    for i in range(0, len(vec)):
        summe += math.pow(abs(vec[i]), p)
    return math.pow(summe/len(vec), (1.0/p))


def get_cartesian_rms(sdf_string1, sdf_string2):
    """Return the optimal RMS after aligning two structures."""
    ref = Chem.MolFromMolBlock(sdf_string1, removeHs=False)
    probe = Chem.MolFromMolBlock(sdf_string2, removeHs=False)
    rms = AllChem.GetBestRMS(ref, probe)
    return rms


def lowest_cartesian(string1, string2, **linked_strings):
    """Select lowest Cartesian RMS for two structures (for nonchiral and
    previously optimized structures)."""
    values = []
    get_cartesian_rms(string1, string2)
    values.append(get_cartesian_rms(string1, string2))
    if linked_strings:
        for string in linked_strings:
            values.append(get_cartesian_rms(string1, string))

    return min(values)


def find_one_in_list(sum_array, list_to_search):
    """Generate a random number and return the corresponding index from a
    list. See the description of the method find_two_in_list."""
    nparray_to_search = np.array(list_to_search)
    rn = sum_array*np.random.rand()
    found = False
    index = 0
    while not found:
        if rn <= nparray_to_search[:index+1].sum(axis=0):
            found = True
        else:
            index += 1
    return index


def find_two_in_list(list_sum, nparray_to_search):
    """A numpy array is mapped to a segment of a line which length is equal to
    1. The lengths of the segments are proportional to the corresponding numpy
    array values. Next, two random numbers between 0 and 1 are generated and
    the segments containing these random numbers are returned."""
    rn1 = list_sum*np.random.rand()
    found1 = False
    index1 = 0
    while not found1:
        if rn1 < nparray_to_search[:index1+1].sum(axis=0):
            found1 = True
        else:
            index1 += 1
    equal = True
    while equal:
        rn2 = list_sum*np.random.rand()
        found2 = False
        index2 = 0
        while not found2:
            if rn2 < nparray_to_search[:index2+1].sum(axis=0):
                found2 = True
            else:
                index2 += 1
        if index2 != index1:
            equal = False
    return index1, index2


def find_closest(numb, list_of_values, periodic=False):
    """For a given number, return the closest value(s) from a given list"""
    all_dist = []
    for value in list_of_values:
        if periodic:
            all_dist.append(min(abs(numb-value), (360-abs(numb-value))))
        else:
            all_dist.append(abs(numb-value))
    m = min(all_dist)
    closest_ind = [i for i, j in enumerate(all_dist) if j == m]
    closest = []
    for ind in closest_ind:
        closest.append(list_of_values[ind])
    return closest


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

    atoms, bonds = get_ind_from_sdfline(sdf_string.split('\n')[3])
    coordinates = np.zeros((atoms, 3))
    bonds_pairs = np.zeros((bonds, 2))
    for i in range(4, atoms+4):
        coordinates[i-4][0] = sdf_string.split('\n')[i].split()[0]
        coordinates[i-4][1] = sdf_string.split('\n')[i].split()[1]
        coordinates[i-4][2] = sdf_string.split('\n')[i].split()[2]
    for i in range(atoms+4, atoms+bonds+4):
        i1, i2 = get_ind_from_sdfline(sdf_string.split('\n')[i])
        bonds_pairs[i-atoms-4][0] = i1
        bonds_pairs[i-atoms-4][1] = i2
    dist = np.zeros((atoms, atoms))
    for x in range(atoms):
        for y in xrange(x, atoms):
            a = np.array(coordinates[int(x)])
            b = np.array(coordinates[int(y)])
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


def get_ind_from_sdfline(sdf_line):
    """Extract the indicies from the sdf string (for molecules with more than
    99 atoms)"""
    l = len(sdf_line.split()[0])
    if l < 4:
        ind1 = int(sdf_line.split()[0])
        ind2 = int(sdf_line.split()[1])
    else:
        list_ind = list(sdf_line.split()[0])
        if len(list_ind) == 5:
            ind1 = int(list_ind[0]+list_ind[1])
            ind2 = int(list_ind[2]+list_ind[3]+list_ind[4])
        if len(list_ind) == 6:
            ind1 = int(list_ind[0]+list_ind[1]+list_ind[2])
            ind2 = int(list_ind[3]+list_ind[4]+list_ind[5])

    return ind1, ind2

# Format conversions


def sdf2aims(sdf_string):
    """Convert a sdf string to a aims string."""
    atoms = get_ind_from_sdfline(sdf_string.split('\n')[3])[0]
    coord = []
    for i in range(4, 4+atoms):
        x = float(sdf_string.split('\n')[i].split()[0])
        y = float(sdf_string.split('\n')[i].split()[1])
        z = float(sdf_string.split('\n')[i].split()[2])
        name = sdf_string.split('\n')[i].split()[3]
        coord.append('%s%10.4f%10.4f%10.4f%4s' % ('atom', x, y, z, name))
        coord.append('\n')
    aims_string = ''.join(coord)
    return aims_string


def sdf2xyz(sdf_string):
    """Convert a sdf string to a xyz string."""
    atoms = get_ind_from_sdfline(sdf_string.split('\n')[3])[0]
    coord = [str(atoms)+('\n')]
    for i in range(4, 4+atoms):
        x = float(sdf_string.split('\n')[i].split()[0])
        y = float(sdf_string.split('\n')[i].split()[1])
        z = float(sdf_string.split('\n')[i].split()[2])
        name = sdf_string.split('\n')[i].split()[3]
        coord.append('\n%2s%10.4f%10.4f%10.4f' % (name, x, y, z))
    coord.append('\n')
    xyz_string = ''.join(coord)
    return xyz_string


def aims2sdf(aims_string, sdf_template_string):
    """Convert a aims string to a sdf string. Template for the sdf string is
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
            c.append('%10.4f%10.4f%10.4f%s%-2s' % (float(line[0]),
                                                   float(line[1]),
                                                   float(line[2]), str(' '),
                                                   line[3]))
            for j in xrange(4, len(line)):
                if j == 4:
                    c.append('%3d' % int(line[j]))
                elif j == len(line)-1:
                    c.append('%3d\n' % int(line[j]))
                else:
                    c.append('%3d' % int(line[j]))
        else:
            c.append(''.join(sdf_form[i])+'\n')

    sdf_string = ''.join(c)
    return sdf_string


def xyz2sdf(xyz_string, sdf_template_string):
    """Convert a xyz string to a sdf string. Template for the sdf string is
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
            c.append('%10.4f%10.4f%10.4f%s%-2s' % (float(line[0]),
                                                   float(line[1]),
                                                   float(line[2]), str(' '),
                                                   line[3]))
            for j in xrange(4, len(line)):
                if j == 4:
                    c.append('%3d' % int(line[j]))
                elif j == len(line)-1:
                    c.append('%3d\n' % int(line[j]))
                else:
                    c.append('%3d' % int(line[j]))
        else:
            c.append(''.join(sdf_form[i])+'\n')

    sdf_string = ''.join(c)
    return sdf_string


def mirror_sdf(sdf_string):
    """Mirror the geometry from a sdf string. Return a new sdf string."""
    atoms = get_ind_from_sdfline(sdf_string.split('\n')[3])[0]
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
            c.append('%10.4f%10.4f%10.4f%s%-2s' % (float(line[0]),
                                                   float(line[1]),
                                                   float(line[2]), str(' '),
                                                   line[3]))
            for j in xrange(4, len(line)):
                if j == 4:
                    c.append('%3d' % int(line[j]))
                elif j == len(line)-1:
                    c.append('%3d\n' % int(line[j]))
                else:
                    c.append('%3d' % int(line[j]))
        else:
            c.append(''.join(sdf_form[i])+'\n')
    mirror_sdf_string = ''.join(c)
    return mirror_sdf_string


def sdf2coord_list(sdf_string):
    """Extract the information about the atoms(coordinates and symbol)
    from the sdf string to a list"""

    atoms = get_ind_from_sdfline(sdf_string.split('\n')[3])[0]
    tmp_arr = sdf_string.splitlines()
    coord_list = []
    for i in range(atoms):
        elem = tmp_arr[4+i].split()
        coord_list.append([elem[3], float(elem[0]), float(elem[1]),
                           float(elem[2])])
    return coord_list


def coord_list2sdf(coord_list, sdf_template_string):
    """Recreate a sdf_string from a coordinate list using a template."""

    atoms = len(coord_list)
    sdf_form = sdf_template_string.splitlines()
    c = []
    cnt = 0
    for i in range(len(sdf_form)):
        if i > 3 and i < 4+atoms:
            line = sdf_form[i].split()
            line[0] = float(coord_list[i-4][1])
            line[1] = float(coord_list[i-4][2])
            line[2] = float(coord_list[i-4][3])
            cnt += 1
            c.append('%10.4f%10.4f%10.4f%s%-2s' % (line[0], line[1], line[2],
                                                   str(' '), line[3]))
            for j in xrange(4, len(line)):
                if j == 4:
                    c.append('%3d' % int(line[j]))
                elif j == len(line)-1:
                    c.append('%3d\n' % int(line[j]))
                else:
                    c.append('%3d' % int(line[j]))
        else:
            c.append(''.join(sdf_form[i])+'\n')

    sdf_string = ''.join(c)
    return sdf_string


def center_mass(list_coord, ignoreH=False):
    """Calculate the center of mass from 3D coordinates and the atomic masses

    Args:
        list_coord - list with lists of coordinates with atom symbols
    Args(optional):
        ignoreH - default=False, if set True - the hydrogens will be ignored
    Returns:
        numpy ndarray of 3D coordinates of the center of mass

    """
    coord = np.zeros([len(list_coord), 4])
    i = 0
    for elem in list_coord:
        atom_type = elem[0]
        coord[i][1:4] = [float(x) for x in elem[1:]]
        if ignoreH and atom_type == "H":
            coord[i][0] = 0
        else:
            coord[i][0] = get_atomic_weight(atom_type)
        i += 1
    center_of_mass = np.average(coord[:, 1:4], axis=0, weights=coord[:, 0])

    return center_of_mass


def translate(list_coord, t_vec):
    """ Perform translation of the 3D coordinates with a vector."""

    t_list_coord = [[x[0], float(x[1])+t_vec[0], float(x[2])+t_vec[1],
                     float(x[3])+t_vec[2]] for x in list_coord]
    return t_list_coord


def rotate_point(point_coord, quaternion):
    """ Rotate a 3D point with a quaternion."""
    if np.linalg.norm(quaternion) != 1.0:
        quaternion = quaternion/np.linalg.norm(quaternion)
    s, x, y, z = quaternion
    p1, p2, p3 = point_coord
    p1_rot = p1*(s*s-z*z-y*y+x*x)+p2*(2*x*y-2*s*z)+p3*(2*x*z+2*y*s)
    p2_rot = p1*(2*x*y+2*s*z)+p2*(s*s-z*z+y*y-x*x)+p3*(2*y*z-2*x*s)
    p3_rot = p1*(2*x*z-2*s*y)+p2*(2*y*z+2*s*x)+p3*(s*s+z*z-y*y-x*x)

    return [p1_rot, p2_rot, p3_rot]


def sum_of_products(list_coord1, list_coord2):
    """Build a matrix of sums of products of coordinates betweent two lists of
    coordinates."""
    sxx, syy, szz = 0, 0, 0
    sxy, sxz, syz = 0, 0, 0
    syx, szx, szy = 0, 0, 0
    for i in range(len(list_coord1)):
        sxx += list_coord1[i][1]*list_coord2[i][1]
        syy += list_coord1[i][2]*list_coord2[i][2]
        szz += list_coord1[i][3]*list_coord2[i][3]
        syx += list_coord1[i][2]*list_coord2[i][1]
        sxy += list_coord1[i][1]*list_coord2[i][2]
        szx += list_coord1[i][3]*list_coord2[i][1]
        sxz += list_coord1[i][1]*list_coord2[i][3]
        szy += list_coord1[i][3]*list_coord2[i][2]
        syz += list_coord1[i][2]*list_coord2[i][3]

    return [[sxx, sxy, sxz], [syx, syy, syz], [szx, szy, szz]]


def get_N_matrix(matrix):
    """Build the N matrix from the sum of products matrix."""

    bigN = np.zeros([4, 4])
    bigN[0][0] = matrix[0][0]+matrix[1][1]+matrix[2][2]
    bigN[1][1] = matrix[0][0]-matrix[1][1]-matrix[2][2]
    bigN[2][2] = -matrix[0][0]+matrix[1][1]-matrix[2][2]
    bigN[3][3] = -matrix[0][0]-matrix[1][1]+matrix[2][2]

    bigN[0][1] = matrix[1][2]-matrix[2][1]
    bigN[1][0] = bigN[0][1]
    bigN[0][2] = matrix[2][0]-matrix[0][2]
    bigN[2][0] = bigN[0][2]
    bigN[0][3] = matrix[0][1]-matrix[1][0]
    bigN[3][0] = bigN[0][3]

    bigN[1][2] = matrix[0][1]+matrix[1][0]
    bigN[2][1] = bigN[1][2]
    bigN[1][3] = matrix[2][0]+matrix[0][2]
    bigN[3][1] = bigN[1][3]
    bigN[2][3] = matrix[1][2]+matrix[2][1]
    bigN[3][2] = bigN[2][3]

    return bigN


def generate_random_quaternion():
    """Generate a random and valid quaternion."""
    while True:
        random_quat = [np.random.uniform(-1, 1) for x in range(4)]
        if np.linalg.norm(random_quat) < 1:
            random_quat = random_quat/np.linalg.norm(random_quat)
        break
    return random_quat


def get_quaternion(list_coord1, list_coord2):
    """Obtain the quaternion between two lists of coordinates."""
    matrix_of_sums = sum_of_products(list_coord1, list_coord2)
    N_matrix = get_N_matrix(matrix_of_sums)
    moments, vectors = np.linalg.eigh(N_matrix)
    q = vectors[:, np.argmax(moments)]
    q = q/np.linalg.norm(q)
    return q


def get_atomic_weight(atom_type):
    """Get the atomic weight from the atom symbol."""
    periodic_table = AllChem.GetPeriodicTable()
    return periodic_table.GetAtomicWeight(atom_type)


def get_inertia_tensor(list_coord):
    """Obtain the inertia tensor from a list of coordinates."""
    m_center = center_mass(list_coord)
    inertia_tensor = np.zeros([3, 3])
    ixx, iyy, izz = 0, 0, 0
    ixy, ixz, iyz = 0, 0, 0
    for item in list_coord:
        atom_mass = get_atomic_weight(item[0])
        ixx += atom_mass*((item[2]-m_center[1])**2+(item[3]-m_center[2])**2)
        iyy += atom_mass*((item[1]-m_center[0])**2+(item[3]-m_center[2])**2)
        izz += atom_mass*((item[1]-m_center[0])**2+(item[2]-m_center[1])**2)

        ixy += atom_mass*((item[2]-m_center[1])*(item[1]-m_center[0]))
        ixz += atom_mass*((item[3]-m_center[2])*(item[1]-m_center[0]))
        iyz += atom_mass*((item[3]-m_center[2])*(item[2]-m_center[1]))
    inertia_tensor = [[ixx, -ixy, -ixz], [-ixy, iyy, -iyz], [-ixz, -iyz, izz]]
    return inertia_tensor


def solve_inertia_tensor(inertia_tensor):
    """Solve the tensor to obtain the principle moments and the axes"""
    moments, tensor_axes = np.linalg.eig(inertia_tensor)
    ind = np.argsort(moments)
    return moments, tensor_axes[:, ind]


def rotate_matrix(list_coord, rot):
    """Rotate the coordinates with a rot_matrix"""
    rot_list = []
    for i in xrange(0, len(list_coord)):
        t1 = list_coord[i][1]*rot[0][0]+list_coord[i][2]*rot[0][1]+list_coord[i][3]*rot[0][2]
        t2 = list_coord[i][1]*rot[1][0]+list_coord[i][2]*rot[1][1]+list_coord[i][3]*rot[1][2]
        t3 = list_coord[i][1]*rot[2][0]+list_coord[i][2]*rot[2][1]+list_coord[i][3]*rot[2][2]
        rot_list.append([list_coord[i][0], t1, t2, t3])
    return rot_list


def rotate_quaternion(list_coord, quaternion):
    """Rotate the coordinates with a quaternion"""
    rot_coord = []
    for item in list_coord:
        rot_pos = []
        rot_pos.append(item[0])
        rot_pos.extend((rotate_point(item[1:], quaternion)))
        rot_coord.append(rot_pos)
    return rot_coord


def make_canon_unique(list_coord):
    """Bring the canonincal from to unique"""
    unique_coord = []
    x1, y1, z1 = list_coord[0][1:]
    for item in list_coord:
        unique_pos = []
        unique_pos.append(item[0])
        unique_pos.extend([np.sign(x1)*item[1], np.sign(y1)*item[2], np.sign(z1)*item[3]])
        unique_coord.append(unique_pos)
    return unique_coord


def get_canonical(list_coord):
    """Obtain the canonical form"""
    tens = get_inertia_tensor(list_coord)
    principle_moments, principle_axes = solve_inertia_tensor(tens)
    axes_transp = principle_axes.T
    canon_list = rotate_matrix(list_coord, axes_transp)
    unique_canon_list = make_canon_unique(canon_list)

    return unique_canon_list
