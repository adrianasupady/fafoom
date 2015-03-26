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
''' Handle the molecule and its 3D structures.'''

from rdkit import Chem
from rdkit.Chem import AllChem
import numpy as np
import os
import glob
from copy import deepcopy


from new_molecule_parametrize import parametrize, template_sdf
from angle import angle_measure, angle_set
from genetic_operations import crossover, mutation_tor, mutation_cistrans
from pyaims import AimsObject
from pyff import FFObject
from pynwchem import NWChemObject
from utilities import boolean, print_output, check_geo_sdf, get_vec, tor_rmsd,\
 string2file, aims2sdf, xyz2sdf, mirror_sdf
from random import randint, choice


class MoleculeDescription:
    """Create the molecule."""
    newline = "NEWLINE"

    def __init__(self, parameter_file=None, **kwargs):
        """Initialize the molecule. Get the parameters from file (if present) or
        read keyword args. The keyword args overwrite the file values."""
        if parameter_file is not None:
            with open(parameter_file) as fin:
                parameter_dict = dict(line.strip().partition(' ')[::2] for line in fin)
            self.smile = parameter_dict['smile']
            self.custom = parameter_dict['custom']
            self.rmsd_cutoff_uniq = float(parameter_dict['rmsd_cutoff_uniq'])
            self.distance_cutoff_1 = float(parameter_dict['distance_cutoff_1'])
            self.distance_cutoff_2 = float(parameter_dict['distance_cutoff_2'])
            self.chiral = parameter_dict['chiral']
            self.rmsd_type = parameter_dict['rmsd_type']
            self.smart_torsion = parameter_dict['smart_torsion']
            self.smart_cistrans = parameter_dict['smart_cistrans']
            self.smart_custom = parameter_dict['smart_custom']
            fin.close()
        for key in kwargs.keys():
            if hasattr(self, str(key)):
                print_output("Overwriting the values for keyword "+str(key))
                print_output("Old value: "+str(getattr(self, str(
                                    key)))+", new value: "+str(kwargs[key]))
            if key not in ["template_sdf_string"]:
                setattr(self, str(key), kwargs[key])
            else:
                setattr(self, str(key), kwargs[key].replace(
                                        MoleculeDescription.newline, "\n"))

    def __repr__(self):
        """Create an unambiguous object representation. The resulting string
        is an one-liner with the newline parameter replacing the original
        '\n' sign in the template sdf_string attribute."""
        repr_list = []
        for att_name in self.__dict__.keys():
            if att_name in ["smile",  "rmsd_type", "smart_torsion", "custom",
                            "chiral", "smart_cistrans", "smart_custom"]:
                repr_list.append('%s="%s"' %
                                 (att_name, getattr(self, att_name)))
            elif att_name in ["rmsd_cutoff_uniq", "distance_cutoff_1",
                              "distance_cutoff_2", "atoms", "bonds",
                              "torsion", "cistrans", "custom_torsion"]:
                repr_list.append('%s=%s' %
                                 (att_name, repr(getattr(self, att_name))))
            elif att_name in ["template_sdf_string"]:
                repr_list.append("%s='%s'" % (
                    att_name, getattr(
                        self, att_name).replace("\n",
                                                MoleculeDescription.newline,)))
            else:
                print_output("Unknown type of attribute "+str(att_name))
        return "%s(%s)" % (self.__class__.__name__, ', '.join(repr_list))

    def __eq__(self, other):
        """Compare all attribute values of two objects. Returns True if all
        values are identical."""
        for att_name in self.__dict__.keys():
            if getattr(self, att_name) == getattr(other, att_name):
                continue
            else:
                return False
        return True

    def get_mol_parameters(self):
        """Assign new attributes (atoms, bonds, torsion, cistrans,
        custom_torsion) to the object."""
        (self.atoms, self.bonds,
         self.torsion, self.cistrans,
         self.custom_torsion) = parametrize(
         self.smile, self.custom, self.smart_torsion,
         self.smart_cistrans, self.smart_custom)
        if self.custom:
            self.torsion = self.custom_torsion

    def create_template_sdf(self):
        """Assign new attribute (template_sdf_string) to the object."""
        self.template_sdf_string = template_sdf(self.smile,
                                                self.distance_cutoff_1,
                                                self.distance_cutoff_2)


class Structure:
    """Create 3D structures."""
    index = 0
    newline = "NEWLINE"

    def __init__(self, arg=None, **kwargs):
        """Initialize the 3D structure: (1) from MoleculeDescription class object
        or from (2) from previously created object of the Structure class. Any
        present and valid keyword args overwrite the old values.
        Warning: the may be more attributes in the (2) case."""
        if isinstance(arg, MoleculeDescription):
            self.mol_info = arg
            Structure.index += 1
            self.index = Structure.index
        elif isinstance(arg, Structure):
            self.mol_info = arg.mol_info
            Structure.index += 1
            self.index = Structure.index
            for att_name in arg.__dict__.keys():
                if att_name != "mol_info" and att_name != "index":
                    setattr(self, str(att_name),
                            deepcopy(getattr(arg, str(att_name))))
        else:
            print_output("Initialization can't be performed. Check the input")
        for key in kwargs.keys():
            if key != "index":
                if hasattr(self, str(key)):
                    print_output("Overwriting the values for keyword "+str(key))
                    print_output("Old value: "+str(getattr(self, str(
                                        key)))+", new value: "+str(kwargs[key]))
                if key not in ["sdf_string", "initial_sdf_string"]:
                    setattr(self, str(key), kwargs[key])
                else:
                    setattr(self, str(key),
                            kwargs[key].replace(Structure.newline, "\n"))

    def __repr__(self):
        """Create an unambiguous object representation. The resulting string
        is an one-liner with the newline parameter replacing the original
        '\n' sign in the sdf_string and initial_sdf_string attribute."""
        repr_list = []
        for att_name in self.__dict__.keys():
            if att_name in ["energy", "index", "values_cistrans", "values_tor",
                            "initial_values_cistrans", "initial_values_tor"]:
                repr_list.append('%s=%s' % (
                    att_name, repr(getattr(self, att_name))))
            elif att_name in ["sdf_string", "initial_sdf_string"]:
                repr_list.append("%s='%s'" % (
                    att_name, getattr(
                        self, att_name).replace("\n",
                                                Structure.newline)))
            elif att_name in ["mol_info"]:
                pass
            else:
                print_output("Unknown type of attribute "+str(att_name))
        return "%s(mol, %s)" % (self.__class__.__name__, ', '.join(repr_list))

    def __str__(self):
        """Return the object index."""
        return "%s %d" % (self.__class__.__name__, self.index)

    def __float__(self):
        """Return the object energy."""
        return float(self.energy)

    def generate_random_structure(self, cistrans1=0, cistrans2=180):
        """Generate a random structure and assign new attributes to the object
        (sdf_string, values_cistrans, values_tor)

        Args(optional):
            cistrans1 (int): first option for the cistrans bond value
            cistrans2 (int): second option for the cistrans bond value
        If the arguments are not specified, the cistrans bond will 0 or
        180 degrees.
        """
        self.values_cistrans, self.values_tor = [], []
        for x in range(len(self.mol_info.cistrans)):
            self.values_cistrans.append(float(choice([cistrans1, cistrans2])))
        for y in range(len(self.mol_info.torsion)):
            self.values_tor.append(randint(0.0, 359.0)-179.0)
        self.sdf_string = angle_set(self.mol_info.template_sdf_string,
                                    self.mol_info.cistrans,
                                    self.mol_info.torsion,
                                    self.values_cistrans, self.values_tor)
        self.values_cistrans, self.values_tor = angle_measure(
                                            self.sdf_string,
                                            self.mol_info.cistrans,
                                            self.mol_info.torsion)

    def generate_structure_from_values(self, values_cis, values_tor):
        """Generate structure from lists of valuse and assign new attributes
        to the object (sdf_string, values_cistrans, values_tor).

        Args:
            values_cis (list): values for the cistrans bonds, put [] if no
            cistrans bonds are present
            values_tor (list): values for the torsions, put [] if no torsions
            are present
        Raises:
            ValueError: if the length of the corresponding lists differ
        """
        if len(values_cis) != len(self.mol_info.cistrans):
            raise ValueError("No match between the number of values to assign \
and the number of cistrans bonds.")
        if len(values_tor) != len(self.mol_info.torsion):
            raise ValueError("No match between the number of values to assign \
and the number of torsions.")
        self.values_cistrans, self.values_tor = values_cis, values_tor
        self.sdf_string = angle_set(self.mol_info.template_sdf_string,
                                    self.mol_info.cistrans,
                                    self.mol_info.torsion,
                                    self.values_cistrans, self.values_tor)
        self.values_cistrans, self.values_tor = angle_measure(
                                            self.sdf_string,
                                            self.mol_info.cistrans,
                                            self.mol_info.torsion)

    def is_geometry_valid(self):
        """Return True if the geometry is valid."""
        check = check_geo_sdf(self.sdf_string, self.mol_info.distance_cutoff_1,
                              self.mol_info.distance_cutoff_2)
        return check

    def __eq__(self, other):
        """Calculate the rmsd for an object pair.

        Returns:
            True, if the objects are 'similar'
        Raises:
            ValueError: if the rmsd type differs between the objects
            ValueErrof: if the rmsd type is unknown
        """
        if self.mol_info.rmsd_type != other.mol_info.rmsd_type:
            raise ValueError("The type of the rmsd differ for the objects")
        if self.mol_info.rmsd_type not in ['cartesian', 'torsional']:
            raise ValueError("Unknown type of rmsd.")

        if self.mol_info.rmsd_type == 'cartesian':
            obj1, obj2 = self, other
            if hasattr(self, "initial_sdf_string"):
                obj1, obj2 = obj2, obj1
            if hasattr(obj1, "initial_sdf_string"):
                raise Exception("Both structures are already relaxed")
            ref = Chem.MolFromMolBlock(obj1.sdf_string, removeHs=False)
            probe = Chem.MolFromMolBlock(obj2.sdf_string, removeHs=False)
            bestrms = AllChem.GetBestRMS(ref, probe)

            if hasattr(obj2, "initial_sdf_string"):
                probe_ini = Chem.MolFromMolBlock(
                        obj2.initial_sdf_string, removeHs=False)
                bestrms_ini = AllChem.GetBestRMS(ref, probe_ini)
                if bestrms_ini < bestrms:
                    bestrms = bestrms_ini
            if not boolean(obj1.mol_info.chiral):
                obj2_mirror_string = mirror_sdf(obj2.sdf_string)
                probe_mirror = Chem.MolFromMolBlock(
                        obj2_mirror_string, removeHs=False)
                bestrms_mirror = AllChem.GetBestRMS(ref, probe_mirror)
                if hasattr(obj2, "initial_sdf_string"):
                    obj2_mirror_string_ini = mirror_sdf(
                                        obj2.initial_sdf_string)
                    probe_mirror_ini = Chem.MolFromMolBlock(
                                obj2_mirror_string_ini, removeHs=False)
                    bestrms_mirror_ini = AllChem.GetBestRMS(
                                    ref, probe_mirror_ini)
                    if bestrms_mirror_ini < bestrms_mirror:
                        bestrms_mirror = bestrms_mirror_ini
                if bestrms_mirror < bestrms:
                    bestrms = bestrms_mirror
            if bestrms > obj1.mol_info.rmsd_cutoff_uniq:
                return False
            else:
                return True
        if self.mol_info.rmsd_type == 'torsional':
            obj1, obj2 = self, other
            if hasattr(self, "initial_sdf_string"):
                obj1, obj2 = obj2, obj1
            if hasattr(obj1, "initial_sdf_string"):
                raise Exception("Both structures are already relaxed.")
            rmsd1 = tor_rmsd(2, get_vec(
                np.concatenate((
                        obj1.values_cistrans, obj1.values_tor), axis=0),
                np.concatenate((
                        obj2.values_cistrans, obj2.values_tor), axis=0)))
            rmsd = rmsd1
            if hasattr(obj2, "initial_sdf_string"):
                rmsd2 = tor_rmsd(2, get_vec(
                    np.concatenate((
                        obj1.values_cistrans, obj1.values_tor), axis=0),
                    np.concatenate((
                        obj2.initial_values_cistrans,
                        obj2.initial_values_tor), axis=0)))
                if rmsd2 < rmsd1:
                    rmsd = rmsd2
            if not boolean(obj1.mol_info.chiral):
                rmsd_mirr1 = tor_rmsd(2, get_vec(
                        np.concatenate((
                            obj1.values_cistrans, obj1.values_tor), axis=0),
                        -1*np.concatenate((
                            obj2.values_cistrans, obj2.values_tor), axis=0)))
                rmsd_mirr2 = tor_rmsd(2, get_vec(
                        np.concatenate((
                            obj1.values_cistrans, obj1.values_tor), axis=0),
                        -1*np.concatenate((
                            obj2.initial_values_cistrans,
                            obj2.initial_values_tor), axis=0)))
                rmsd_mirr = min(rmsd_mirr1, rmsd_mirr2)
                if rmsd_mirr < rmsd:
                    rmsd = rmsd_mirr
            if rmsd > obj1.mol_info.rmsd_cutoff_uniq:
                return False
            else:
                return True

    def __cmp__(self, other):
        """Compare two object basing on their energy values."""
        return cmp(self.energy, other.energy)

    def send_to_blacklist(self, address, array):
        """Append the structure to dedicated array. Write it to file and store.

        Args:
           address: name of directory to store the structure (file); if the
           directory is absent it will be created
           array: the array to append to
        Raise:
            NameError: if the array not defined
        """
        try:
            os.mkdir(address)
        except OSError:
            pass
        array.append(self)

        cnt_black = len(glob.glob(address+"/*.sdf"))
        string2file(self.sdf_string,
                    str(address)+"/black_"+str(cnt_black)+".sdf")
        if not self.mol_info.chiral:
            mirror_sdf_string = mirror_sdf(self.sdf_string)
            cnt_black = len(glob.glob(address+"/*.sdf"))
            string2file(mirror_sdf_string,
                        str(address)+"/black_"+str(cnt_black)+".sdf")

    def perform_aims(self, parameter_file, execution_string, dirname):
        """Generate the FHI-aims input, run FHI-aims, store the output, assign
        new attributes (initial_sdf_string, initial_values_cistrans,
        initial_values_tor) and update attribute values (sdf_string,
        values_cistrans, values_tor)."""

        aims_object = AimsObject(parameter_file)
        aims_object.generate_input(self.sdf_string)
        aims_object.build_storage(dirname)
        aims_object.run_aims(execution_string)
        aims_object.clean_and_store()
        self.energy = aims_object.get_energy()
        self.initial_sdf_string = self.sdf_string
        self.initial_values_cistrans = self.values_cistrans
        self.initial_values_tor = self.values_tor
        self.sdf_string = aims2sdf(
                    aims_object.get_aims_string_opt(),
                    self.mol_info.template_sdf_string)
        self.values_cistrans, self.values_tor = angle_measure(
                                            self.sdf_string,
                                            self.mol_info.cistrans,
                                            self.mol_info.torsion)

    def perform_nwchem(self, parameter_file, execution_string):
        """Generate the NWChem input, run NWChem, assign new attributes
        (initial_sdf_string, initial_values_cistrans, initial_values_tor) and
        update attribute values (sdf_string, values_cistrans, values_tor)."""
        nwchem_object = NWChemObject(parameter_file)
        nwchem_object.clean()
        nwchem_object.generate_input(self.sdf_string)
        nwchem_object.run_nwchem(execution_string)
        nwchem_object.clean()
        self.energy = nwchem_object.get_energy()
        self.initial_sdf_string = self.sdf_string
        self.initial_values_cistrans = self.values_cistrans
        self.initial_values_tor = self.values_tor
        self.sdf_string = xyz2sdf(
                    nwchem_object.get_xyz_string_opt(),
                    self.mol_info.template_sdf_string)
        self.values_cistrans, self.values_tor = angle_measure(
                                            self.sdf_string,
                                            self.mol_info.cistrans,
                                            self.mol_info.torsion)

    def perform_ff(self, parameter_file):
        """Generate the force-field input, run force=field calculation, assign
        new attributes (initial_sdf_string, initial_values_cistrans,
        initial_values_tor) and update attribute values (sdf_string,
        values_cistrans, values_tor)."""
        ff_object = FFObject(parameter_file)
        ff_object.run_ff(self.sdf_string)
        self.energy = ff_object.get_energy()
        self.initial_sdf_string = self.sdf_string
        self.initial_values_cistrans = self.values_cistrans
        self.initial_values_tor = self.values_tor
        self.sdf_string = ff_object.get_sdf_string_opt()
        self.values_cistrans, self.values_tor = angle_measure(
                                            self.sdf_string,
                                            self.mol_info.cistrans,
                                            self.mol_info.torsion)

    def crossover(self, other):
        """Perform the crossover."""
        start1, start2 = crossover(self.values_cistrans,
                                   other.values_cistrans)
        end1, end2 = crossover(self.values_tor, other.values_tor)
        child1 = Structure(self.mol_info)
        child1.generate_structure_from_values(start1, end1)
        child2 = Structure(self.mol_info)
        child2.generate_structure_from_values(start2, end2)
        return child1, child2

    def mutation_tor(self, max_mutations_torsions):
        """Mutate the torsion values list."""
        dt = len(self.mol_info.torsion)
        if dt > 0:
            mut_list = mutation_tor(self.values_tor,
                                    max_mutations_torsions)
            self.generate_structure_from_values(self.values_cistrans, mut_list)

    def mutation_cistrans(self, max_mutations_cistrans):
        """Mutate the cistrans bonds list."""
        dc = len(self.mol_info.cistrans)
        if dc > 0:
            mut_list = mutation_cistrans(self.values_cistrans,
                                         max_mutations_cistrans)
            self.generate_structure_from_values(mut_list, self.values_tor)
