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

'''This module is a wrapper for OpenBabel force-field utility'''
import openbabel
import pybel

from utilities import xyz_to_xyz_readstring

class FFObject():
    '''This class allows for creating and handling force-field objects '''
    def __init__(self, parameter_file):
        ''' __init__(parameter_file) -> None 
        
        Required parameters in the parameter_file:
        force_field
        
        Optional parameters for the SteepestDescent method used for optimization:
        steps (default=3000)
        threshold (default=1e-10)
        '''
        
        with open(parameter_file) as fin:
            parameter_dict = dict(line.strip().partition(' ')[::2] for line in fin)
        self.force_field = parameter_dict['force_field']
        
        s =[('steps', 3000), ('threshold', 1.0e-10)] 
        for k, v in s:
            if k not in parameter_dict: parameter_dict[k]=v
                
        self.steps = parameter_dict['steps']
        self.threshold = parameter_dict['threshold']
    
    def run_ff(self, xyz_string):
        ''' run_ff(xyz_string) -> None
        
        This function performs a force field optimization. '''
        molecule = openbabel.OBMol()
        xyz_to_xyz_readstring(molecule,xyz_string)
        ff = openbabel.OBForceField.FindForceField(self.force_field)
        ff.Setup(molecule)
        ff.SteepestDescent(int(self.steps), float(self.threshold))
        ff.UpdateCoordinates(molecule)
        self.energy= float('{0:.4f}'.format(ff.Energy(molecule)))
        pybelmol = pybel.Molecule(molecule)
        self.xyz_string_optimized = pybelmol.write(format = 'xyz')

    def get_energy(self):
        '''  get_energy() -> energy
        
        This function retrieves the energy value or raises an exception if it is not present. '''
        
        if self.energy == None: 
            raise Exception("The calculation wasn't performed yet. You need to run ff first.")
        else: return self.energy
    def get_xyz_string_opt(self):
        ''' get_aims_string_opt() -> aims_string_opt
        
        This function retrieves the optimized string in the FHI-aims format and raises an exception if it is not present.'''
        if self.xyz_string_optimized == None:
            raise Exception("The calculation wasn't performed yet. You need to run ff first.")
        else: return self.xyz_string_optimized
    
