import numpy as np
import sys

from fafoom import MoleculeDescription, Structure, selection, print_output,\
    remover_dir, set_default, file2dict
import fafoom.run_utilities as run_util

# Decide for restart or a simple run.
opt = run_util.simple_or_restart()
p_file = sys.argv[1]
# Build a dictionary from two section of the parameter file.
params = file2dict(p_file, ['GA settings', 'Run settings'])

dict_default = {'energy_var': 0.001, 'selection': "roulette_wheel",
                'fitness_sum_limit': 1.2, 'popsize': 10,
                'prob_for_crossing': 1.0, 'max_iter': 30,
                'iter_limit_conv': 20, 'energy_diff_conv': 0.001}
# Set defaults for parameters not defined in the parameter file.
params = set_default(params, dict_default)
energy_function = run_util.detect_energy_function(params)

cnt_max = 200
population, blacklist = [], []
min_energy = []

if opt == "simple":
    mol = MoleculeDescription(p_file)
    # Assign the permanent attributes to the molecule.
    mol.get_parameters()
    mol.create_template_sdf()
    # Check for potential degree of freedom related parameters.
    linked_params = run_util.find_linked_params(mol, params)
    print_output("Number of atoms: "+str(mol.atoms))
    print_output("Number of bonds: "+str(mol.bonds))
    for dof in mol.dof_names:
        print_output("Number of identified "+str(dof)+": " +
                     str(len(getattr(mol, dof))))
        print_output("Identified "+str(dof)+": "+str(getattr(mol, dof)))

    print_output("___Initialization___")
    cnt = 0
    # Generate sensible and unique 3d structures.
    while len(population) < params['popsize'] and cnt < cnt_max:
        print_output("New trial")
        str3d = Structure(mol)
        str3d.generate_structure()
        if not str3d.is_geometry_valid():
            print_output("The geometry of "+str(str3d)+" is invalid.")
            cnt += 1
            continue
        if str3d not in blacklist:
            name = "initial_%d" % (len(population))
            # Perform the local optimization
            run_util.optimize(str3d, energy_function, params, name)
            run_util.check_for_kill()
            str3d.send_to_blacklist(blacklist)
            population.append(str3d)
            print_output(str(str3d)+", energy: "+str(float(str3d)) +
                         ", was added to the population")
            run_util.relax_info(str3d)
            cnt += 1
        else:
            print_output("Geomerty of "+str(str3d)+" is fine, but already "
                         "known.")
            cnt += 1
    if cnt == cnt_max:
        print_output("The allowed number of trials for building the "
                     "population has been exceeded. The code terminates.")
        sys.exit(0)
    print_output("___Initialization completed___")
    population.sort()
    print_output("Initial population after sorting: ")
    for i in range(len(population)):
        print_output(str(population[i])+" "+str(float(population[i])))
    min_energy.append(population[0].energy)
    print_output("Blacklist: " + ', '.join([str(v) for v in blacklist]))
    iteration = 0


if opt == "restart":
    # Reconstruct the molecule, population, blacklist and the state of the run.
    print_output(" \n ___Restart will be performed___")
    with open("backup_mol.dat", 'r') as inf:
        mol = eval(inf.readline())
    inf.close()
    with open("backup_population.dat", 'r') as inf:
        for line in inf:
            population.append(eval(line))
    inf.close()
    with open("backup_blacklist.dat", 'r') as inf:
        for line in inf:
            blacklist.append(eval(line))
    inf.close()
    with open("backup_min_energy.dat", 'r') as inf:
        for line in inf:
            min_energy.append(eval(line))
    inf.close()
    with open("backup_iteration.dat", 'r') as inf:
        iteration_tmp = eval(inf.readline())
    inf.close()
    linked_params = run_util.find_linked_params(mol, params)
    population.sort()
    for i in range(len(population)):
        print_output(str(population[i])+" "+str(float(population[i])))
    print_output("Blacklist: " + ', '.join([str(v) for v in blacklist]))
    iteration = iteration_tmp+1
    print_output(" \n ___Reinitialization completed___")
    remover_dir('generation_'+str(iteration)+'_child1')
    remover_dir('generation_'+str(iteration)+'_child2')


def mutate_and_relax(candidate, name, iteration, cnt_max, **kwargs):
    print_output("__%s__" % name)
    found = False
    cnt = 0
    while found is False and cnt < cnt_max:
        candidate_backup = Structure(candidate)
        candidate.mutate(**kwargs)
        print_output("%s after mutation: " % name + str(candidate))
        run_util.str_info(candidate)
        if not candidate.is_geometry_valid():
            print_output(" The geometry of %s is invalid." % name)
            cnt += 1
            # Rebuild the structure
            candidate = candidate_backup
            continue

        if candidate not in blacklist:
            name = "generation_%d_%s" % (iteration, name)
            run_util.optimize(candidate, energy_function, params, name)
            run_util.check_for_kill()
            candidate.send_to_blacklist(blacklist)
            print_output(str(candidate)+":, energy: "+str(float(
                candidate))+", is temporary added to the population")
            run_util.relax_info(candidate)
            found = True
            population.append(candidate)
        else:
            print_output("Geomerty of "+str(candidate)+" is fine, but already "
                         "known.")
            cnt += 1
            candidate = candidate_backup
        if cnt == cnt_max:
            raise Exception("The allowed number of trials for generating"
                            " a unique child has been exceeded.")

while iteration < params['max_iter']:
    print_output(" \n ___Start of iteration " + str(iteration) + "___")
    (parent1, parent2, fitness) = selection(population, params['selection'],
                                            params['energy_var'],
                                            params['fitness_sum_limit'])
    param = np.random.rand()
    cnt = 0
    while param < params['prob_for_crossing'] and cnt < cnt_max:
        child1, child2 = Structure.crossover(parent1, parent2)
        if child1.is_geometry_valid() and child2.is_geometry_valid():
            print_output("Crossover outcome: "+str(child1)+(", ")+str(child2))
            break
        else:
            print_output("The geometries created via crossover are invalid.")
            cnt += 1
            continue
    else:
        child1, child2 = Structure(parent1), Structure(parent2)
        print_output("No crossover was performed. Children are copies of "
                     "parents: " + str(child1) + (": ") + str(child1) +
                     (", ") + str(child2) + (": ") + str(child2))
        # Delete inherited attributes.
        for child in child1, child2:
            attr_list = ["initial_sdf_string", "energy"]
            for attr in attr_list:
                delattr(child, attr)
            for dof in child.dof:
                delattr(dof, "initial_values")

    run_util.str_info(child1)
    run_util.str_info(child2)

    try:
        mutate_and_relax(child1, "child1", iteration, cnt_max, **linked_params)
    except Exception as exc:
        print_output(exc)
        sys.exit(0)
    try:
        mutate_and_relax(child2, "child2", iteration, cnt_max, **linked_params)
    except Exception as exc:
        print_output(exc)
        sys.exit(0)
    population.sort()
    print_output("Sorted population: " + ', '.join([
        str(v) for v in population]))
    del population[-1]
    del population[-1]
    print_output("Sorted population after removing two structures with highest"
                 " energy: " + ', '.join([str(v) for v in population]))
    min_energy.append(population[0].energy)
    print_output("Lowest energy of the population: %.3f" % min_energy[-1])
    print_output("Lowest energies in run: "+str(min_energy))
    run_util.perform_backup(mol, population, blacklist, iteration, min_energy)
    run_util.check_for_convergence(iteration, params, min_energy)
    run_util.check_for_kill()
    iteration += 1
