import numpy as np
import sys
import glob

from fafoom import MoleculeDescription, Structure, selection, print_output, backup

p_file = sys.argv[1]

with open(p_file) as f:
            p_dict = dict(line.strip().partition(' ')[::2] for line in f)
popsize = int(p_dict['popsize'])
cistrans1, cistrans2 = float(p_dict['cistrans1']), float(p_dict['cistrans2'])
black_dir = str(p_dict['black_dir'])
max_iter = int(p_dict['max_iter'])
energy_var = float(p_dict['energy_var'])
selection_type = str(p_dict['selection'])
fitness_sum_limit = str(p_dict['fitness_sum_limit'])
prob_for_crossing = float(p_dict['prob_for_crossing'])
prob_for_mut_cistrans = float(p_dict['prob_for_mut_cistrans'])
prob_for_mut_torsions = float(p_dict['prob_for_mut_torsions'])
max_mutations_cistrans = int(p_dict['max_mutations_cistrans'])
max_mutations_torsions = int(p_dict['max_mutations_torsions'])
nwchem_call = str(p_dict['nwchem_call'])
iter_limit_conv = int(p_dict['iter_limit_conv'])
energy_diff_conv = float(p_dict['energy_diff_conv'])
energy_wanted = float(p_dict['energy_wanted'])

cnt_max = 100

mol=MoleculeDescription(p_file)
mol.get_mol_parameters()
mol.create_template_sdf()

population, blacklist = [], []
min_energy = []
print_output("___Initialization___")
cnt =0 
while len(population) < popsize and cnt < cnt_max:
	print_output("New trial")
	str3d=Structure(mol)
	str3d.generate_random_structure(cistrans1,cistrans2)
	if not str3d.is_geometry_valid():
		print_output("The geometry of "+str(str3d)+" is invalid.")
		cnt +=1
		continue
	if str3d not in blacklist:
		name= "initial_%d" %(len(population))
		str3d.perform_nwchem(p_file, nwchem_call)
		str3d.send_to_blacklist(black_dir, blacklist)
		population.append(str3d)
		print_output(str(str3d)+", energy: "+str(float(str3d))+", was added to the population")
		cnt +=1
	else:
		print_output("Geomerty of "+str(str3d)+" is fine, but is already known.")
		cnt +=1
if cnt == cnt_max:
	print_output("The allowed number of trials for building the population has been exceeded. The code terminates.")
	sys.exit(0)
	
print_output("___Initialization completed___")		
population.sort()
print_output("Initial population after sorting: ")
for i in range(len(population)):
	print_output(str(population[i])+" "+str(float(population[i])))
min_energy.append(population[0].energy)
print_output("Blacklist: " +', '.join([str(v) for v in blacklist]))

iteration=0

def mutate_and_relax(candidate, name, prob_for_mut_torsions, prob_for_mut_cistrans, max_mutations_torsions, max_mutations_cistrans, iteration, cnt_max):
	print_output("__%s__" %name)
	found=False
	cnt = 0 
	while found == False and cnt < cnt_max:
		mut1, mut2 = np.random.rand(), np.random.rand()
		candidate_backup = Structure(candidate)	
		if mut1 < prob_for_mut_torsions:
			candidate.mutation_tor(max_mutations_torsions) 
			print_output("%s after mutation in torsions: " %name +str(candidate))
		if mut2 < prob_for_mut_cistrans:
			candidate.mutation_cistrans(max_mutations_cistrans) 
			print_output("%s after mutation in cistrans: " %name+str(candidate)) 
		if not candidate.is_geometry_valid(): 
			print_output(" The geometry of %s is invalid." %name)
			cnt +=1
			candidate = candidate_backup
			continue
		if candidate not in blacklist:  
			name= "generation_%d_%s" %(iteration, name)
			candidate.perform_nwchem(p_file, nwchem_call)
			candidate.send_to_blacklist(black_dir, blacklist) 
			print_output(str(candidate)+":, energy: "+str(float(candidate))+", is temporary added to the population")
			found = True
			population.append(candidate)
			
		else:
			print_output("Geomerty fine, but the structure already known")
			cnt +=1
			candidate = candidate_backup
		if cnt == cnt_max:
			raise Exception("The allowed number of trials for generating a unique child has been exceeded. The code terminates.")	


while (iteration < max_iter):

	print_output(" \n ___Start of iteration "+str(iteration)+"___")
	(parent1, parent2, fitness) = selection(population, selection_type, energy_var, fitness_sum_limit)
	param=np.random.rand()
	cnt = 0 
	while param < prob_for_crossing and cnt < cnt_max: 
		child1,child2 = Structure.crossover(parent1, parent2)
		if child1.is_geometry_valid() and child2.is_geometry_valid():
			print_output("Crossover outcome: "+str(child1)+(", ")+str(child2))
			break
		else: 
			print_output("The geometries created in the crossover are invalid.")
			cnt +=1
			continue
	else: 
		child1, child2 = Structure(parent1), Structure(parent2)
		print_output("No crossover was performed. Children are copies of parents: "+str(child1)+(": ")+str(child1)+(", ")+str(child2)+(": ")+str(child2))
		for child in child1, child2:
			for attr in "initial_sdf_string", "energy", "initial_values_cistrans", "initial_values_tor":
				delattr(child, attr)
	try:
		mutate_and_relax(child1, "child1", prob_for_mut_torsions, prob_for_mut_cistrans,max_mutations_torsions, max_mutations_cistrans, iteration,cnt_max)
	except Exception as exc:
		print_output(exc)
		sys.exit(0)
	try:
		mutate_and_relax(child2, "child2", prob_for_mut_torsions, prob_for_mut_cistrans,max_mutations_torsions, max_mutations_cistrans, iteration,cnt_max)
	except Exception as exc:
		print_output(exc)
		sys.exit(0)
	
	population.sort()
	print_output("Sorted population: " +', '.join([str(v) for v in population]))
	del population[-1]
	del population[-1]
	
	print_output("Sorted population after removing two structures with highest energy: " +', '.join([str(v) for v in population]))
	print_output("Lowest energy of the population: %.3f" %population[0].energy)
	min_energy.append(population[0].energy)
	print_output("Lowest energies in run: "+str(min_energy))
	
	def perform_backup():
		backup("backup_mol.dat", mol)
		backup("backup_population.dat", population)
		backup("backup_blacklist.dat", blacklist)
		backup("backup_iteration.dat", iteration)
		backup("backup_min_energy.dat", min_energy)	
	
		
	if iteration >= iter_limit_conv-1:
		print_output("Checking for convergence")
				
		if population[0].energy < energy_wanted or abs(min_energy[iteration+1]-min_energy[iteration+1-iter_limit_conv]) < energy_diff_conv:
			perform_backup()
			print_output("Converged")
			killfile = open("kill.dat", "w")
			killfile.close	
			sys.exit(0)
		else:
			perform_backup()
			print_output("Not converged yet")
	
	if iteration == max_iter-1:
		perform_backup()
		print_output("Maximal number of iterations reached. The code terminates")
		killfile = open("kill.dat", "w")
		killfile.close
		sys.exit(0)
	else:
		perform_backup() 
		print_output("Next iteration will be perfomed")
		
	if len(glob.glob("*/kill.dat")) == 0 and len(glob.glob("kill.dat")) == 0 : 
		iteration+=1
	else:
		print_output("Kill.dat file discovered. The code terminates")
		sys.exit(0)

