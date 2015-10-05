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
''' Collection of diverse run controlling utilites '''
import glob
import sys

from utilities import print_output, remover_file, remover_dir, backup


def simple_or_restart():
    """ Select the type of run. If the all backup files are present the run
    will be restarted. Otherwise, the directory will be cleaned and a new run
    will be started."""

    for_restart = ["backup_population.dat", "backup_mol.dat",
                   "backup_min_energy.dat", "backup_iteration.dat",
                   "backup_blacklist.dat"]
    opt = "restart"
    for filename in for_restart:
        if glob.glob(filename):
            pass
        else:
            opt = "simple"
    if opt == "simple":
        print_output("Cleaning up the directory")
        for d in glob.glob("initial_*"):
            remover_dir(d)
        for d in glob.glob("generation_*_child*"):
            remover_dir(d)
        remover_dir("blacklist")
        for f in ["mol.sdf", "control.in", "geometry.in", "output.txt",
                  "result.out", "kill.dat"]:
            remover_file(f)
        for f in for_restart:
            remover_file(f)
    if opt == "restart":
        remover_file("kill.dat")
    return opt


def str_info(struct):
    """ Prints the information about the structure to the output file"""
    print_output(struct)
    for dof in struct.dof:
        print_output("Values of "+str(dof.type)+": " + str(dof.values))


def relax_info(struct):
    """ Prints the information about the structure to the output file after
    the local optimization."""
    print_output(struct)
    for dof in struct.dof:
        print_output("Values of "+str(dof.type)+" before: " +
                     str(dof.initial_values)+" and after local opt: " +
                     str(dof.values))


def check_for_kill():
    """ Check if the kill.dat file is present in the directory or in the
    subdirectories. If discoveed the run will be aborted."""
    if len(glob.glob("*/kill.dat")) == 0 and len(glob.glob("kill.dat")) == 0:
        pass
    else:
        print_output("Kill.dat file discovered. The code terminates")
        sys.exit(0)


def detect_energy_function(params):
    """ Detect the energy function that will be used for local optimization."""
    if 'energy_function' not in params:
        print_output("The energy function hasn't been defined."
                     " The code terminates")
        sys.exit(0)
    else:
        if params['energy_function'] in ['aims', 'FHI-aims', 'FHIaims']:
            print_output("Local optimization will be performed with FHI-aims.")
            energy_function = "aims"
        elif params['energy_function'] in ['nwchem', 'NWChem']:
            print_output("Local optimization will be performed with NWChem.")
            energy_function = "nwchem"
        elif params['energy_function'] in ['ORCA', 'Orca', 'orca']:
            print_output("Local optimization will be performed with ORCA.")
            energy_function = "orca"
        elif params['energy_function'] in ['ff', 'force_field', 'RDKit',
                                           'rdkit']:
            print_output("Local optimization will be performed with RDKit.")
            energy_function = "ff"
        else:
            print_output("Unknown type of energy function."
                         " The code terminates.")
            sys.exit(0)
        return energy_function


def optimize(structure, energy_function, params, name=None):
    """Perform local optimization."""
    if energy_function == "aims":
        structure.perform_aims(params['sourcedir'], params['aims_call'], name)
    elif energy_function == "nwchem":
        structure.perform_nwchem(params['functional'], params['basis_set'],
                                 params['nwchem_call'])
    elif energy_function == "orca":
        linked_params = {}
        for key in ["chargemult", "nprocs", "optsteps"]:
            if key in params:
                linked_params[str(key)] = params[str(key)]
        structure.perform_orca(params['commandline'],
                               params['memory'],
                               params['orca_call'], **linked_params)
    elif energy_function == "ff":
        linked_params = {}
        for key in ["steps", "force_tol", "energy_tol"]:
            if key in params:
                linked_params[str(key)] = params[str(key)]
        structure.perform_ff(params['force_field'], **linked_params)


def perform_backup(mol, population, blacklist, iteration, min_energy):
    """Write object representation to files for a future restart."""
    backup("backup_mol.dat", mol)
    backup("backup_population.dat", population)
    backup("backup_blacklist.dat", blacklist)
    backup("backup_iteration.dat", iteration)
    backup("backup_min_energy.dat", min_energy)


def find_linked_params(mol, params):
    """ Check for parameters assosciated with the degrees of freedom."""
    linked_params = {}
    for dof_name in mol.dof_names:
        if 'prob_for_mut_'+str(dof_name) in params:
            linked_params['prob_for_mut_'+str(dof_name)] = \
                params['prob_for_mut_'+str(dof_name)]
        if 'max_mutations_'+str(dof_name) in params:
            linked_params['max_mutations_'+str(dof_name)] = \
                params['max_mutations_'+str(dof_name)]
    return linked_params


def check_for_convergence(iteration, params, min_energy):
    """Check the run for convergence"""
    if iteration >= params['iter_limit_conv']-1:
        print_output("Checking for convergence")
        d = abs(min_energy[iteration+1]-min_energy[iteration + 1 -
                                                   params['iter_limit_conv']])
        if 'energy_wanted' in params:
            if min_energy[-1] < params['energy_wanted'] or \
               d < params['energy_diff_conv']:
                print_output("Converged")
                killfile = open("kill.dat", "w")
                killfile.close()
                sys.exit(0)
            else:
                print_output("Not converged yet")
        else:
            if d < params['energy_diff_conv']:
                print_output("Converged")
                killfile = open("kill.dat", "w")
                killfile.close()
                sys.exit(0)
            else:
                print_output("Not converged yet")
    if iteration == params['max_iter']-1:
        print_output("Max. number of iterations reached. The code terminates")
        killfile = open("kill.dat", "w")
        killfile.close()
        sys.exit(0)
    else:
        print_output("Next iteration will be perfomed")
