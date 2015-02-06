#    Copyright 2015 Adriana Supady 
#    adriana.supady@gmail.com
#
#    This file is part of fafoom.
#
#   Fafoom is free software: you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation, either version 3 of the License, or
#   (at your option) any later version.
#
#   Fafoom is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#   along with fafoom.  If not, see <http://www.gnu.org/licenses/>.

''' This module collects genetic operations'''

import numpy as np 
import random
import os
from copy import deepcopy

from utilities import find_two_in_list

def selection(pop_list, selection_type, energy_range, fitness_sum_limit):
    '''selection(pop_list, selection_type, energy_range, fitness_sum_limit) -> object, object, vector
    
    This function returns two selected objects from the population, together with the fitness values vector. No copies are created. 
    '''
    best_e=float(pop_list[0])
    worst_e=float(pop_list[-1])
    fitness = np.zeros(len(pop_list))
    
    for j in range(len(pop_list)):
        if ( worst_e-best_e > energy_range):   # only when energy differs among the candidates
            fitness[j] = (worst_e-float(pop_list[j]))/(worst_e-best_e)
        else: fitness[j] = 1.0 # equal fitness if diversity low
    
    fitness_sum=fitness.sum()
    
    if selection_type == "roulette_wheel":
        if fitness_sum > fitness_sum_limit:
            x,y = find_two_in_list(fitness_sum, fitness)       
            parent1 = pop_list[x-1]
            parent2 = pop_list[y-1]
        else: 
            parent1 = pop_list[0] # if the sum of all fitness values is under the limit, the best and a random object are selected
            parent2 = pop_list[int(np.ceil(np.random.rand()*(len(pop_list)-1)))]
    
    if selection_type == "roulette_wheel_reverse":
        if  fitness_sum > fitness_sum_limit: # in order to prevent numerical problems
            
            fitness_rev = np.zeros(len(fitness))
            for i in range(len(fitness)):
                fitness_rev[-(i+1)]=fitness[i] # swapping of fitness values
            x,y = find_two_in_list(fitness_sum, fitness_rev)
        
            parent1 = pop_list[x-1]
            parent2 = pop_list[y-1]
        else: 
            parent1 = pop_list[-1]
            parent2 = pop_list[int(np.ceil(np.random.rand()*(len(pop_list)-1))-1)]
    
    if selection_type == "random":
        
        parents_ind = random.sample(xrange(len(pop_list)), 2)
        parent1 = pop_list[parents_ind[0]]
        parent2 = pop_list[parents_ind[1]]
    
    return (parent1, parent2, fitness)        

    
def crossover_vec(vector1,vector2):
    '''crossover_vec(vector1, vector2) -> new_vector1, new_vector2
    
    This function cuts to vectors in a random position and exchanges the parts so that two new vectors can be created    
    '''  
    d=len(vector1)
    cross_point = int(np.ceil(np.random.rand()*d)-1)
    start1 = vector1[:cross_point]
    ind = 0
    end1 = np.zeros(d - cross_point)
    for n in range(len(end1)):
        end1[ind] = vector2[cross_point+n]
        ind +=1
    start2 = vector2[:cross_point]
    ind = 0
    end2 = np.zeros(d - cross_point)
    for n in range(len(end2)):
        end2[ind] = vector1[cross_point+n]
        ind +=1
    new_vector1 = np.append(start1,end1)
    new_vector2 = np.append(start2,end2)
    return new_vector1, new_vector2
        

    
def mutation_tor_vec(vector, max_mutations_torsions):
    '''mutation_tor_vec(vector, max_mutations_torsions) -> vector
    
    This function performs mutation in the vector containing the torsion values.
    The number of mutations is a random integer, but not higher than the max_mutations_torsions.
    For each mutation, a random position in the vector is selected. A random integer between(-180,179) is generated and assigned to the selected value.
    '''
    mut_numb = int(np.ceil(np.random.rand()*max_mutations_torsions))
    newtor1 = []
    for p in range(mut_numb):
        newtor1.append(np.ceil(np.random.rand()*359.0)-180.0)
        mutp2 = int(np.ceil(np.random.rand()*len(vector)))
        vector[mutp2-1] = newtor1[p]
    return vector
    
def mutation_cistrans_vec(vector, max_mutations_cistrans):
	'''mutation_cistrans_vec(vector, max_mutations_cistrans) -> vector
    
    This function performs mutation in the vector containing the cistrans values.
    The number of mutations is a random integer, but not higher than the max_mutations_cistrans.
    For each mutation, a random position in the vector is selected. Values higher than 90 (or lower than 90) are replaced with 0.
    Other values are replaced with 180.
    '''
	
	mut_numb = int(np.ceil(np.random.rand()*max_mutations_cistrans))
	for p in range(mut_numb):
		mutp1 = int(np.ceil(np.random.rand()*len(vector)))
		wert = vector[mutp1-1]
		if wert > 90.0 or wert < -90.0:
			vector[mutp1-1] = 0.0
		else:
			vector[mutp1-1] = 180.0
	return vector


