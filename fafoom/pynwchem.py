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

'''Wrapper for NWCHEM'''

import glob
import os
import subprocess

from utilities import sdf2xyz

prefix_name = 'geo'


class NWChemObject():
    '''Create and handle NWCHEM objects '''
    def __init__(self, parameter_file):
        """Get the parameters from file.
        Args(required):
            parameter_file (str): name of the parameter file
        Raises:
            KeyError: if the functional or basis set is not defined
        """
        with open(parameter_file) as fin:
            parameter_dict = dict(line.strip().partition(' ')[::2] for line in fin)
        fin.close()
        for key in ['functional', 'basis_set']:
            if key not in parameter_dict:
                raise KeyError("The functional or basis_set is not defined.")
        self.functional = parameter_dict['functional']
        self.basis_set = parameter_dict['basis_set']

    def generate_input(self, sdf_string):
        """Create input files for NWChem.
        Args:
            sdf_string (str)
        """
        xyz_string = sdf2xyz(sdf_string)
        coord = xyz_string.split('\n')
        string1 = 'title "Molecule DFT geometry optimization"\ngeometry\n'
        string2 = 'end \nbasis\n * library '+str(self.basis_set)
        string3 = 'end\ndft\n  iterations 50\n  xc '+str(self.functional)
        string4 = 'end\ndriver\n  xyz '+prefix_name
        string5 = '  loose\nend\ntask dft optimize'
        with open('nwchem_molecule.nw', 'w') as f:
            f.write(string1)
            f.write('\n'.join(coord[2:]))
            f.write('\n'+string2)
            f.write('\n'+string3)
            f.write('\n'+string4)
            f.write('\n'+string5)
        f.close()

    def run_nwchem(self, execution_string):
        """Run NWChem and write output to 'result.out'. The optimized
        geometry is written to 'geo.xyz'.

        Warning: this function uses subprocessing to invoke the run.
        The subprocess's shell is set to TRUE.
        Args:
            execution_string (str): e.g. mpirun -n 4 nwchem
        Raises:
            OSError: if nwchem_molecule.nw not present in the working directory
        """
        if os.path.exists('nwchem_molecule.nw') == False:
                raise OSError("Required input file not present.")
        nwchem = subprocess.Popen(
            execution_string+str(" nwchem_molecule.nw"),
            stdout=subprocess.PIPE, shell=True)
        out = subprocess.Popen(
            ['cat'], stdin=nwchem.stdout,
            stdout=open('result.out', 'w'), shell=True)
        out.wait()
        searchfile = open("result.out", "r")
        for line in searchfile:
            if "Total DFT energy" in line:
                a = line.split("=")[-1].split(" ")[-1].split('\n')
                energy_tmp = float('{0:.4f}'.format(float(a[0])))
        searchfile.close()
        self.energy = energy_tmp
        all_images = glob.glob(prefix_name+'-*')
        # sort based on the numerical value of the index of the image file
        all_images.sort(key=lambda x: int(x.split('-')[-1].split('.')[0]))
        with open(all_images[-1], 'r') as f:
            self.xyz_string_opt = f.read()
        f.close()

    def get_energy(self):
        """Get the energy of the molecule.

        Returns:
            energy (float)
        Raises:
            AttributeError: if energy hasn't been calculated yet
        """
        if not hasattr(self, 'energy'):
            raise AttributeError("The calculation wasn't performed yet.")
        else:
            return self.energy

    def get_xyz_string_opt(self):
        """Get the optimized xyz string.

        Returns:
            optimized xyz string (str)
        Raises:
            AttributeError: if the optimization hasn't been performed yet
        """
        if not hasattr(self, 'xyz_string_opt'):
            raise AttributeError("The calculation wasn't performed yet.")
        else:
            return self.xyz_string_opt

    def clean(self):
        """Clean the working direction after the NWChem calculation has been
        completed.
        """
        for f in glob.glob("nwchem_molecule.*"):
            os.remove(f)

        for f in glob.glob(prefix_name+'-*'):
            os.remove(f)
