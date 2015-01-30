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
''' This module collects diverse help/convert functions '''
import os
import glob
import openbabel
import pickle
import pybel
import numpy as np
import math
import random
import shutil

def backup(filename, object):
    ''' This function write the representation of an object (or group) of objects to a file'''
    with open(filename, 'w') as outf:
        if hasattr(object,"__len__"):
            for i in range(len(object)):
                outf.write("%s\n"%repr(object[i]))
        else: outf.write("%s\n"%repr(object))


def boolean(string):
    ''' This function recovers the boolean value from a string'''
    if string in ["False", "false", "FALSE"]:
        return False
    if string in ["True", "true", "TRUE"]:
        return True
    raise Exception("Cannot be converted to a boolean type")
    

def get_vec(vec1,vec2):
    '''This function compares two vectors and calculate their difference. No symmetry handling included yet!'''
    diff_vec=np.zeros((1,len(vec1)))
    for i in range(1, len(vec1)+1):
        tor_diff=abs(vec1[i-1]-vec2[i-1])
        diff_vec[0][i-1]=min(abs(tor_diff),abs(360-tor_diff))/180.0
    return diff_vec[0]

def tor_rmsd(p,vec):
    ''' This function return the torsional rmsd values from a difference vector, p is the norm'''
    summe=0
    for i in range(0,len(vec)):
        summe+=math.pow(abs(vec[i]),p)
    return math.pow(summe/len(vec), (1.0/p))

def find_two_in_list(vector_sum, vector):
    ''' A list is mapped to a segment of a line which length is equal to 1. The lengths of the segments are proportional to the 
    corresponding list values. Next, two random numbers between 0 and 1 are generated and the segements containing these random numbers are returned. '''
    rn1 = vector_sum*np.random.rand()
    found1 = False
    index = 1
    while not found1:
        if rn1 < vector[:index].sum(axis=0):
            found1 = index
        else: index += 1
    equal = True
    while equal:
        rn2 = vector_sum*np.random.rand()
        found2 = False
        index = 1
        while not found2:
            if rn2 < vector[:index].sum(axis=0):
                found2 = index
            else: index += 1
        if found2 != found1:
            equal = False

    return found1, found2

def xyz_to_fhiaims_string(xyz_string):
    mymol = pybel.readstring("xyz", xyz_string)
    return mymol.write("fhiaims")

def xyz_to_fhiaims_file(xyz_string, filename):
    mymol = pybel.readstring("xyz", xyz_string)
    return mymol.write("fhiaims", filename, overwrite=True)

def fhiaims_to_xyz_string(fhiaims_string):
    mymol = pybel.readstring("fhiaims", fhiaims_string)
    return mymol.write("xyz")

def fhiaims_to_xyz_file(fhiaims_string, filename):
    mymol = pybel.readstring("fhiaims", fhiaims_string)
    return mymol.write("xyz", filename, overwrite=True)

def xyz_to_sdf_string(xyz_string):
    mymol = pybel.readstring("xyz", xyz_string)
    return mymol.write("sdf")

def xyz_to_sdf_file(xyz_string, filename):
    mymol = pybel.readstring("xyz", xyz_string)
    return mymol.write("sdf", filename, overwrite=True)

def sdf_to_xyz_string(sdf_string):
    mymol = pybel.readstring("sdf", sdf_string)
    return mymol.write("xyz")

def sdf_to_xyz_file(sdf_string, filename):
    mymol = pybel.readstring("sdf", sdf_string)
    return mymol.write("xyz", filename, overwrite=True)

def sdf_to_sdf_readstring(molecule, sdf_string):
    conv = openbabel.OBConversion() 
    conv.SetInAndOutFormats("sdf", "sdf")
    conv.ReadString(molecule, sdf_string)
    
def sdf_to_xyz_writestring(molecule):
    conv = openbabel.OBConversion() 
    conv.SetInAndOutFormats("sdf", "xyz")
    return conv.WriteString(molecule)

def xyz_to_xyz_readstring(molecule, xyz_string):
    conv = openbabel.OBConversion() 
    conv.SetInAndOutFormats("xyz", "xyz")
    conv.ReadString(molecule, xyz_string)

def mirrorer(filename, filename_mirr):
    with open(filename_mirr, 'w') as df:
        with open(filename, 'r') as f:
            for line in f:
                newline = line.split()
                if len(newline)>2:
                    newline[1]=str(-1.0*float(newline[1]))
                    newline[2]=str(-1.0*float(newline[2]))
                    newline[3]=str(-1.0*float(newline[3]))
                df.write(' '.join(newline)+'\n')          
    df.close()
    f.close() 
    
def print_output(text):
    if os.path.isfile("output.txt"):
        f = open("output.txt", "a")
        f.write(str(text)+'\n')
        f.close()
    else: 
        f = open("output.txt", "w")
        f.write(str(text)+'\n')
        f.close()

def remover_file(instance):
    try: os.remove(instance)
    except OSError: pass
    
def remover_dir(instance):
    try: shutil.rmtree(instance)
    except OSError: pass
