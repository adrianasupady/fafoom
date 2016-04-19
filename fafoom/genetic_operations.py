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

"""Collection of genetic operations."""
from __future__ import division
import numpy as np
import random

from utilities import (
    find_closest,
    find_one_in_list,
    find_two_in_list,
    print_output
)


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
            parent1 = pop_list[x]
            parent2 = pop_list[y]
        else:  # if the sum is below the limit, best and a random are selected
            parent1 = pop_list[0]
            parent2 = pop_list[int(random.choice(range(1, len(pop_list))))]
    if selection_type == "roulette_wheel_reverse":
        if fitness_sum > fitness_sum_limit:  # in order to prevent num problems
            fitness_rev = np.zeros(len(fitness))
            for i in range(len(fitness)):
                fitness_rev[-(i+1)] = fitness[i]  # swapping of fitness values
            x, y = find_two_in_list(fitness_sum, fitness_rev)
            parent1 = pop_list[x]
            parent2 = pop_list[y]
        else:
            parent1 = pop_list[-1]
            parent2 = pop_list[int(random.choice(range(1, len(pop_list)))-1)]
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
        two lists (converted numpy arrays)
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
        new_list1 = list(np.append(start1, end1))
        new_list2 = list(np.append(start2, end2))
        return new_list1, new_list2
    else:
        return list1, list2


def mutation(list_for_mut, max_mutations, options, weights=None,
             periodic=False):
    """Mutate a list of values.

    Args(required):
        list_for_mut (list): list of values to be mutated
        max_mut: maximal number of changes to be made
        options: list of options for new values
    Args(optional):
        weights (list): weights for options
        periodic (bool)
    Returns:
        mutated list
    """

    if not isinstance(max_mutations, int):
        raise TypeError("The max. number of mutations needs to be an integer")
    if max_mutations < 0:
        raise ValueError("The max. number of mutations cannot be negative")
    mut_numb = random.randint(1, min(max_mutations, len(list_for_mut)))
    pos = random.sample(range(len(list_for_mut)), mut_numb)
    for p in pos:
        current_value = list_for_mut[p]
        banned = find_closest(current_value, options, periodic)
        cnt = 0
        while cnt < 100:
            if weights is None:
                new_value = random.sample(options, 1)[0]
            else:
                new_value = options[find_one_in_list(sum(weights), weights)]

            if len(banned) != len(options):
                if new_value not in banned:
                    list_for_mut[p] = new_value
                    break
                else:
                    cnt += 1
            else:
                # for very rare cases, i.e. if there are only two options
                # and the value is the mean of them!
                list_for_mut[p] = new_value
                break

    return list_for_mut
