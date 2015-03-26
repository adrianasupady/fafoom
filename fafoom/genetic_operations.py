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

"""Collection of genetic operations."""

import numpy as np
import random

from utilities import find_two_in_list, print_output


def selection(pop_list, selection_type, energy_range, fitness_sum_limit):
    """Select two objects from a list.

    Args:
        pop_list (list): sorted list of population objects
        selection_type (string): method of selection
        energy_range (float): criterion for distinctive energy values
        fitness_sum_limit (float): criterion for distinctive fitness values
    Returns:
        selected objects and a list of all fitness values
    Raises:
        ValueError: if the input parameters are outside the range
    If selection_type is not properly defined, it will be set to random.
    """
    if len(pop_list) < 2:
        raise ValueError("The list needs to have at least two members.")
    if selection_type not in ["roulette_wheel",
                              "roulette_wheel_reverse", "random"]:
        print_output("Method not found. Selection method set to random")
        selection_type = "random"
    if energy_range < 0:
        raise ValueError("The energy range cannot be negative")
    if fitness_sum_limit < 1:
        raise ValueError("The fitness_sum_limit cannot be smaller than 1")
    best_e = float(pop_list[0])
    worst_e = float(pop_list[-1])
    fitness = np.zeros(len(pop_list))
    for j in range(len(pop_list)):
        if worst_e-best_e > energy_range:  # only when energy differs
            fitness[j] = (worst_e-float(pop_list[j]))/(worst_e-best_e)
        else:
            fitness[j] = 1.0  # equal fitness if diversity low
    fitness_sum = fitness.sum()
    if selection_type == "roulette_wheel":
        if fitness_sum > fitness_sum_limit:
            x, y = find_two_in_list(fitness_sum, fitness)
            parent1 = pop_list[x-1]
            parent2 = pop_list[y-1]
        else:  # if the sum is below the limit, best and a random are selected
            parent1 = pop_list[0]
            parent2 = pop_list[int(np.ceil(np.random.rand()*(len(pop_list)-1)))]
    if selection_type == "roulette_wheel_reverse":
        if fitness_sum > fitness_sum_limit:  # in order to prevent num problems
            fitness_rev = np.zeros(len(fitness))
            for i in range(len(fitness)):
                fitness_rev[-(i+1)] = fitness[i]  # swapping of fitness values
            x, y = find_two_in_list(fitness_sum, fitness_rev)
            parent1 = pop_list[x-1]
            parent2 = pop_list[y-1]
        else:
            parent1 = pop_list[-1]
            parent2 = pop_list[int(np.ceil(np.random.rand()*(len(pop_list)-1))-1)]
    if selection_type == "random":
        parents_ind = random.sample(xrange(len(pop_list)), 2)
        parent1 = pop_list[parents_ind[0]]
        parent2 = pop_list[parents_ind[1]]
    return parent1, parent2, fitness


def crossover(list1, list2):
    """Exchange parts of two lists.

    Args:
        list1 (list): list of values
        list2 (list): list of values
    Returns:
        two numpy arrays
    Raises:
        ValueError: if the length of the lists differ
    """
    if len(list1) != len(list2):
        raise ValueError("No length match between the lists")
    d = len(list1)
    if d > 0:
        cross_point = int(np.ceil(np.random.rand()*d)-1)
        start1 = list1[:cross_point]
        ind = 0
        end1 = np.zeros(d - cross_point)
        for n in range(len(end1)):
            end1[ind] = list2[cross_point+n]
            ind += 1
        start2 = list2[:cross_point]
        ind = 0
        end2 = np.zeros(d - cross_point)
        for n in range(len(end2)):
            end2[ind] = list1[cross_point+n]
            ind += 1
        new_list1 = np.append(start1, end1)
        new_list2 = np.append(start2, end2)
        return new_list1, new_list2
    else:
        return list1, list2


def mutation_tor(list_for_mut, max_mutations_torsions):
    """Mutate a list of torsion values.

    Args:
        list_for_mut (list): list of values
        max_mutations_torsions (int): maximal allowed number of mutations
    Returns:
        mutated list; integers from (-180,179) are allowed as new values
    Raises:
        TypeError: if the max_mutations_torsions is not an integer
        ValueError: if the max_mutations_torsions is negative
    """
    if not isinstance(max_mutations_torsions, int):
        raise TypeError("The maximal number of mutations needs to be an integer")
    if max_mutations_torsions < 0:
        raise ValueError("The maximal number of mutations cannot be negative")
    mut_numb = int(np.ceil(np.random.rand()*max_mutations_torsions))
    newtor1 = []
    for p in range(mut_numb):
        newtor1.append(np.ceil(np.random.rand()*359.0)-180.0)
        mutp2 = int(np.ceil(np.random.rand()*len(list_for_mut)))
        list_for_mut[mutp2-1] = newtor1[p]
    return list_for_mut


def mutation_cistrans(list_for_mut, max_mutations_cistrans):
    """Mutate a list of cistrans values.

    Args:
        list_for_mut (list): list of values
        max_mutations_cistrans (int): maximal allowed number of mutations
    Returns:
        mutated list; values from (-90,90> are changed to 180; otherwise to 0
    Raises:
        TypeError: if the max_mutations_cistrans is not an integer
        ValueError: if the max_mutations_cistrans is negative
    """
    if not isinstance(max_mutations_cistrans, int):
        raise TypeError("The maximal number of mutations needs to be an integer")
    if max_mutations_cistrans < 0:
        raise ValueError("The maximal number of mutations cannot be negative")
    mut_numb = int(np.ceil(np.random.rand()*max_mutations_cistrans))
    for p in range(mut_numb):
        mutp1 = int(np.ceil(np.random.rand()*len(list_for_mut)))
        wert = list_for_mut[mutp1-1]
        if wert > 90.0 or wert < -90.0:
            list_for_mut[mutp1-1] = 0.0
        else:
            list_for_mut[mutp1-1] = 180.0
    return list_for_mut
