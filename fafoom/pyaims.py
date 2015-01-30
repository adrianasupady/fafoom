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

'''This module is a wrapper for FHI-aims.'''
import shutil
import os
import subprocess
import itertools

from utilities import remover_file, xyz_to_fhiaims_string


class AimsObject():
    '''This class allows for creating and handling objects FHI-aims objects'''
    def __init__(self, parameter_file):
        ''' __init__(parameter_file) -> None 
        
        parameter_file needs to contain following parameters:
        
        sourcedir (directory containing the control.in file) '''
        with open(parameter_file) as fin:
            parameter_dict = dict(line.strip().partition(' ')[::2] for line in fin)
        self.sourcedir = parameter_dict['sourcedir']
        
    def generate_input(self,xyz_string):
        '''generate_input(xyz_string) -> None
        
        This function reformats the XYZ-string to a FHI-aims-string, creates the geometry.in file 
        and copies the control.in file in to the working directory'''
        self.aims_string = xyz_to_fhiaims_string(xyz_string)
        file=open("geometry.in", "w")
        file.write(self.aims_string)
        file.close()
        name='control.in'
        src = os.path.join(self.sourcedir, name)
        shutil.copy(src, os.getcwd())
        
    def build_storage(self, dirname):
        ''' build_storage(dirname) -> None
        
        This function creates a directory for future storing of the FHI-aims input and output'''
        
        self.dirname = dirname
        os.mkdir(self.dirname)
        shutil.copy('geometry.in', self.dirname)
        shutil.copy('control.in', self.dirname)

    def run_aims(self,execution_string):
        ''' This function runs FHI-aims. The input files need to be present in the working directory.
        It uses the execution string to invoke the run.
        
        Warning: this function uses subprocessing to invoke the MPI-run. The subprocess's shell is set to TRUE.
        
        Immediately after launching, an output file ('result.out') will be created. 
        Once the FHI-aims run is finished, the output file is searched for the energy values.
        If the SCF cycle of the FHI-aims hasn't converged, the kill.dat file will be created
        in the working directory (i.e. the code will terminate after this iteration).        
        
        The optimized geometry is written to 'geometry.out' and subsequently read from this file and assigned to a new object attribute (aims_string_opt)'''
    
        aims = subprocess.Popen(execution_string, stdout=subprocess.PIPE, shell=True)
        out = subprocess.Popen(['cat'], stdin=aims.stdout, stdout=open('result.out', 'w'), shell=True)
        out.wait()
        
        searchfile = open("result.out", "r")
        for line in searchfile:
            if "Total energy of the DFT / Hartree-Fock s.c.f. calculation      :" in line:
                a = line.split(" ") 
                self.energy = float('{0:.4f}'.format(float(a[-2]))) #energy update from the DFT 
            if  "*** scf_solver: SCF cycle not converged." in line:
                killfile = open("kill.dat", "w")
                killfile.close()
        searchfile.close()
        searchfile = open("result.out", "r")        
        for i,line in enumerate(searchfile,1):
            if "Final atomic structure:" in line:
                line_number=int(i)                
        searchfile.close()

        atoms =len([word for word in self.aims_string.split() if word.startswith('atom')])
        with open('geometry.out', 'w') as file_geo:
            try:
                with open('result.out') as f:
                    for line in itertools.islice(f, line_number+1, line_number+1+atoms):
                        file_geo.write(line)
            except IOError or NameError or UnboundLocalError:        
                with open('geometry.in') as f:
                    for line in f:
                        file_geo.write(line)
   
        with open('geometry.out', 'r') as f:
            self.aims_string_opt = f.read()
        f.closed
        
    def get_energy(self):
        '''  get_energy() -> energy
        
        This function retrieves the energy value or raises an exception if it is not present. '''
        
        if self.energy == None: 
            raise Exception("The calculation wasn't performed yet. You need to run run_aims first.")
        else: return self.energy
    def get_aims_string_opt(self):
        ''' get_aims_string_opt() -> aims_string_opt
        
        This function retrieves the optimized string in the FHI-aims format and raises an exception if it is not present.'''
        if self.aims_string_opt == None:
            raise Exception("The calculation wasn't performed yet. You need to run run_aims first.")
        else: return self.aims_string_opt
        
    def clean_and_store(self):
        ''' This function cleans the working directory once the FHI-aims run is finished. It moves the files to dedicated directory.'''
        os.remove('geometry.in')
        os.remove('control.in')
        shutil.copy('result.out', self.dirname)
        os.remove('result.out')
        remover_file('geometry.in.next_step')
        shutil.copy('geometry.out',self.dirname)
        os.remove('geometry.out')