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
''' Handle the molecule and its 3D structures.'''

from copy import deepcopy

from get_parameters import (
    create_dof_object,
    get_atoms_and_bonds,
    get_positions,
    template_sdf
)
from genetic_operations import crossover
from pyaims import AimsObject
from pyff import FFObject
from pynwchem import NWChemObject

from utilities import (
    aims2sdf,
    check_geo_sdf,
    file2dict,
    lowest_cartesian,
    mirror_sdf,
    print_output,
    set_default,
    xyz2sdf

)
import random
#from random import randint, choice


class MoleculeDescription:
    """Create the molecule."""
    newline = "NEWLINE"

    def __init__(self, parameter_file=None, **kwargs):
        """Initialize the molecule. Get the parameters from file (if present)
        or read keyword args. The keyword args overwrite the file values."""
        params = {}
        if parameter_file is not None:
            params = file2dict(parameter_file, ['Molecule'])

        else:
            for key in kwargs.keys():
                if key not in ["template_sdf_string"]:
                    params[str(key)] = kwargs[key]
                else:
                    params[str(key)] = kwargs[key].replace(
                        MoleculeDescription.newline, "\n")

        dict_default = {'rmsd_type': "cartesian", 'distance_cutoff_1': 1.3,
                        'distance_cutoff_2': 2.15, 'rmsd_cutoff_uniq': 0.2,
                        'chiral': True, 'optimize_torsion': True,
                        'smart_torsion':
                        "[*]~[!$(*#*)&!D1]-&!@[!$(*#*)&!D1]~[*]"}

        params = set_default(params, dict_default)

        for key in params:
            if not hasattr(self, str(key)):
                setattr(self, str(key), params[key])

    def __repr__(self):
        """Create an unambiguous object representation. The resulting string
        is an one-liner with the newline parameter replacing the original
        '\n' sign in the template sdf_string attribute."""
        repr_list = []
        for att_name in self.__dict__.keys():

            if type(self.__dict__[att_name]) in [str] and \
               att_name != "template_sdf_string":
                repr_list.append('%s="%s"' %
                                 (att_name, getattr(self, att_name)))

            elif type(self.__dict__[att_name]) in [int, float, bool, list]:

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
        for att_name in other.__dict__.keys():
            if getattr(other, att_name) == getattr(self, att_name):
                continue
            else:
                return False
        return True

    def get_parameters(self):
        """Assign permanent attributes (number of atoms, number of bonds and
        degrees of freedom related attributes) to the object."""
        self.atoms, self.bonds = get_atoms_and_bonds(self.smile)
        self_copy = deepcopy(self)
        dof_names = []
        for attr, value in self_copy.__dict__.iteritems():
            if str(attr).split('_')[0] == "optimize" and value:
                type_of_dof = str(attr).split('_')[1]
                linked_attr = {}
                for attr, value in self_copy.__dict__.iteritems():
                    if type_of_dof in str(attr).split('_'):
                        linked_attr[attr] = value
                pos = get_positions(type_of_dof, self.smile, **linked_attr)
                if len(pos) != 0:
                    setattr(self, type_of_dof, pos)
                    dof_names.append(type_of_dof)
                else:
                    print_output("The degree to optimize: "+str(type_of_dof) +
                                 " hasn't been found.")
        setattr(self, "dof_names", dof_names)

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
        """Initialize the 3D structure: (1) from MoleculeDescription class
        object or from (2) from previously created object of the Structure
        class. Any present and valid keyword args overwrite the old values.
        Warning: there may be more attributes in the (2) case."""
        if isinstance(arg, MoleculeDescription):

            self.mol_info = arg
            Structure.index += 1
            self.index = Structure.index
            dof = []
            for i in self.mol_info.dof_names:
                new_obj = create_dof_object(str(i), getattr(self.mol_info, i))
                dof.append(new_obj)
            setattr(self, "dof", dof)

        elif isinstance(arg, Structure):

            self.mol_info = arg.mol_info
            Structure.index += 1
            self.index = Structure.index
            for att_name in arg.__dict__.keys():

                if att_name not in ["mol_info", "index"]:
                    setattr(self, str(att_name),
                            deepcopy(getattr(arg, str(att_name))))

        else:
            print_output("Initialization can't be performed. Check the input")

        for key in kwargs.keys():
            if key != "index":
                if hasattr(self, str(key)):
                    print_output("Overwriting the value for keyword "+str(key))
                    print_output("Old value: "+str(getattr(self, str(key))) +
                                 ", new value: "+str(kwargs[key]))
                if key in ["sdf_string", "initial_sdf_string"]:
                    setattr(self, str(key),
                            kwargs[key].replace(Structure.newline, "\n"))

                elif key.split('_')[0] in self.mol_info.dof_names:
                    for dof in self.dof:
                        if key.split('_')[0] == dof.type:
                            if key.split('_')[1] == 'initial':
                                setattr(dof, 'initial_values', kwargs[key])
                            if key.split('_')[1] == 'values':
                                setattr(dof, 'values', kwargs[key])
                else:
                    setattr(self, str(key), kwargs[key])

    def __repr__(self):
        """Create an unambiguous object representation. The resulting string
        is an one-liner with the newline parameter replacing the original
        '\n' sign in the sdf_string and initial_sdf_string attribute."""
        repr_list = []
        for att_name in self.__dict__.keys():

            if att_name in ["sdf_string", "initial_sdf_string"]:
                repr_list.append("%s='%s'" % (
                    att_name, getattr(
                        self, att_name).replace("\n",
                                                Structure.newline)))
            else:
                if type(self.__dict__[att_name]) in [str]:
                    repr_list.append('%s=%s' % (
                                     att_name, repr(getattr(self, att_name))))
                elif type(self.__dict__[att_name]) in [int, float, bool]:
                    repr_list.append('%s=%s' % (
                                     att_name, repr(getattr(self, att_name))))
                elif att_name == 'dof':
                    for dof in self.dof:
                        repr_list.append('%s_%s=%s' % (
                                         dof.type, "values", repr(dof.values)))
                        try:
                            repr_list.append('%s_%s=%s' % (
                                             dof.type, "initial_values",
                                             repr(dof.initial_values)))
                        except:
                            pass
                elif att_name == 'mol_info':
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

    def generate_structure(self, values={}):
        """Generate a 3D structures. If no values are passed, a random
        structure will be generated (weights, associated with the degrees of
        freedom, will be taken into account)."""

        new_string = deepcopy(self.mol_info.template_sdf_string)
        for dof in self.dof:
            if dof.type in values.keys():
                new_string = dof.apply_on_string(new_string, values[dof.type])
            else:
                if hasattr(self.mol_info, "weights_"+str(dof.type)):
                    weights = getattr(self.mol_info, "weights_"+str(dof.type))
                    dof.get_weighted_values(weights)
                else:
                    dof.get_random_values()
                new_string = dof.apply_on_string(new_string)
        self.sdf_string = new_string
        for dof in self.dof:
            dof.update_values(self.sdf_string)

    def is_geometry_valid(self):
        """Return True if the geometry is valid."""
        check = check_geo_sdf(self.sdf_string, self.mol_info.distance_cutoff_1,
                              self.mol_info.distance_cutoff_2)
        return check

    def __eq__(self, other):
        """Deicide, if the objects are equal based on the rms values.

        Returns:
            True, if the objects are 'similar'
        Raises:
            ValueError: if the rmsd type differs between the objects
            ValueErrof: if the rmsd type is unknown (supported options are
            'cartesian' and 'internal_coord')
        """

        if self.mol_info.rmsd_type != other.mol_info.rmsd_type:
            raise ValueError("The type of the rmsd differ for the objects")
        if self.mol_info.rmsd_type not in ['cartesian', 'internal_coord']:
            raise ValueError("Unknown type of rmsd.")

        obj1, obj2 = self, other
        if hasattr(self, "initial_sdf_string"):
            obj1, obj2 = obj2, obj1
        if hasattr(obj1, "initial_sdf_string"):
            raise Exception("Both structures are already relaxed.")

        if obj1.mol_info.rmsd_type == 'cartesian':
            linked_strings = {}

            if hasattr(obj2, "initial_sdf_string"):
                n_str = str(obj2.initial_sdf_string)
                linked_strings[n_str] = obj2.initial_sdf_string

            if not obj1.mol_info.chiral:
                n_str = str(mirror_sdf(obj2.sdf_string))
                linked_strings[n_str] = mirror_sdf(obj2.sdf_string)
                if hasattr(obj2, "initial_sdf_string"):
                    n_str = str(mirror_sdf(obj2.initial_sdf_string))
                    linked_strings[n_str] = mirror_sdf(obj2.initial_sdf_string)

            bestrms = lowest_cartesian(obj1.sdf_string, obj2.sdf_string,
                                       **linked_strings)

            if bestrms > obj1.mol_info.rmsd_cutoff_uniq:
                return False
            else:
                return True

        if obj1.mol_info.rmsd_type == 'internal_coord':
            all_bool = []
            for dof1, dof2 in zip(obj1.dof, obj2.dof):
                all_bool.append(dof1.is_equal(dof2,
                                              obj1.mol_info.rmsd_cutoff_uniq,
                                              obj1.mol_info.chiral))

            if False in all_bool:
                return False
            else:
                return True

    def __cmp__(self, other):
        """Compare two object basing on their energy values."""
        return cmp(self.energy, other.energy)

    def send_to_blacklist(self, array):
        """Append the structure to dedicated array.

        Args:
           array: the array to append to
        Raise:
            NameError: if the array not defined
        """

        array.append(self)

    def perform_aims(self, sourcedir, execution_string, dirname):
        """Generate the FHI-aims input, run FHI-aims, store the output, assign
        new attributes and update attribute values."""

        aims_object = AimsObject(sourcedir)
        aims_object.generate_input(self.sdf_string)
        aims_object.build_storage(dirname)
        success = aims_object.run_aims(execution_string)
        if success:
            aims_object.clean_and_store()
            self.energy = aims_object.get_energy()
            self.initial_sdf_string = self.sdf_string
            self.sdf_string = aims2sdf(aims_object.get_aims_string_opt(),
                                       self.mol_info.template_sdf_string)

            for dof in self.dof:
                setattr(dof, "initial_values", dof.values)
                dof.update_values(self.sdf_string)
        else:
            print_output("The FHI-aims relaxation failed")

    def perform_nwchem(self, functional, basis_set, execution_string):
        """Generate the NWChem input, run NWChem, assign new attributes and
        update attribute values."""
        nwchem_object = NWChemObject(functional, basis_set)
        nwchem_object.clean()
        nwchem_object.generate_input(self.sdf_string)
        nwchem_object.run_nwchem(execution_string)
        nwchem_object.clean()
        self.energy = nwchem_object.get_energy()
        self.initial_sdf_string = self.sdf_string
        self.sdf_string = xyz2sdf(nwchem_object.get_xyz_string_opt(),
                                  self.mol_info.template_sdf_string)

        for dof in self.dof:
            setattr(dof, "initial_values", dof.values)
            dof.update_values(self.sdf_string)

    def perform_ff(self, force_field, **kwargs):
        """Generate the force-field input, run force=field calculation, assign
        new attributes and update attribute values."""
        ff_object = FFObject(force_field, **kwargs)
        ff_object.run_ff(self.sdf_string)
        self.energy = ff_object.get_energy()
        self.initial_sdf_string = self.sdf_string
        self.sdf_string = ff_object.get_sdf_string_opt()

        for dof in self.dof:
            setattr(dof, "initial_values", dof.values)
            dof.update_values(self.sdf_string)

    def crossover(self, other):
        """Perform the crossover."""

        child1 = Structure(self.mol_info)
        child2 = Structure(self.mol_info)

        for dof_par1, dof_par2, dof_child1, dof_child2 in zip(self.dof,
                                                              other.dof,
                                                              child1.dof,
                                                              child2.dof):
            if dof_par1.type == dof_par2.type:
                a, b = crossover(getattr(dof_par1, "values"),
                                 getattr(dof_par2, "values"))
                setattr(dof_child1, "values", a)
                setattr(dof_child2, "values", b)

        for child in child1, child2:
            new_string = deepcopy(child.mol_info.template_sdf_string)
            for dof in child.dof:
                new_string = dof.apply_on_string(new_string, dof.values)
            child.sdf_string = new_string
            for dof in child.dof:
                dof.update_values(child.sdf_string)

        return child1, child2

    def mutate(self, **kwargs):

        def call_mut(dof, max_mutations=None, weights=None):
            print_output("Performing mutation for: "+str(dof.type))
            if max_mutations is not None:
                if hasattr(self.mol_info, "weights_"+str(dof.type)):
                    weights = getattr(self.mol_info, "weights_"+str(dof.type))
                    dof.mutate_values(max_mutations, weights)
                else:
                    dof.mutate_values(max_mutations=max_mutations)
            else:
                if hasattr(self.mol_info, "weights_"+str(dof.type)):
                    weights = getattr(self.mol_info, "weights_"+str(dof.type))
                    dof.mutate_values(weights=weights)
                else:
                    dof.mutate_values()

        for dof in self.dof:
            if 'prob_for_mut_'+str(dof.type) in kwargs:
                if random.random() < kwargs['prob_for_mut_'+str(dof.type)]:
                    if 'max_mutations_'+str(dof.type) in kwargs:
                        call_mut(dof, kwargs['max_mutations_'+str(dof.type)])
                    else:
                        call_mut(dof)
            else:
                if 'max_mutations_'+str(dof.type) in kwargs:
                    call_mut(dof, kwargs['max_mutations_'+str(dof.type)])
                else:
                    call_mut(dof)

        new_string = deepcopy(self.sdf_string)
        for dof in self.dof:
            new_string = dof.apply_on_string(new_string, dof.values)
        self.sdf_string = new_string
        for dof in self.dof:
            dof.update_values(self.sdf_string)
