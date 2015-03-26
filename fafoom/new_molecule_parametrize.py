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

"""Create the molecule from the smile code."""

from rdkit import Chem
from rdkit.Chem import AllChem
from copy import copy
from operator import itemgetter

from utilities import check_geo_sdf, boolean


def parametrize(smile, custom="false",
                smart_torsion="[C,N,O]~[!$(*#*)&!D1]-&!@[!$(*#*)&!D1]~[C,N,O]",
                smart_cistrans="C~[$(C=O)]-[$(NC)]~[C]", smart_custom=""):
    """Build the molecule and return parameters.

    Args(required):
        smile (str): one-line representation of the molecule
    Args(optional):
        custom (str): if true, the selection of torsions will be filtered
        smart_torsion (str): search pattern for rotatable bonds [4 atoms!]
        smart_cistrans (str): search pattern for cistrans bonds [4 atoms!]
        smart_custom (str): search pattern torsions to be ignored [4 atoms!]
    Returns:
        number of atoms, bonds, lists of tuples: defining torsions, defining
        cistrans bonds and, optionally, defining filtered torsions
    Raises:
        ValueError: if the molecule cannot be build form the smile
        ValueError: if custom true and smart_custom not defined
    """
    if boolean(custom) and smart_custom == "":
        raise ValueError("The custom_pattern is not defined")
    mol = Chem.MolFromSmiles(smile)
    if mol is None:
        raise ValueError("The smile is invalid")
    pattern_tor = Chem.MolFromSmarts(smart_torsion)
    pattern_cistrans = Chem.MolFromSmarts(smart_cistrans)
    torsion = list(mol.GetSubstructMatches(pattern_tor))
    cistrans = list(mol.GetSubstructMatches(pattern_cistrans))
    custom_torsion = []

    def ig(x):
        return itemgetter(x)

    def cleaner(list_to_clean):
        for_remove = []
        for x in reversed(range(len(list_to_clean))):
            for y in reversed(range(x)):
                ix1, ix2 = ig(1)(list_to_clean[x]), ig(2)(list_to_clean[x])
                iy1, iy2 = ig(1)(list_to_clean[y]), ig(2)(list_to_clean[y])
                if (ix1 == iy1 and ix2 == iy2) or (ix1 == iy2 and ix2 == iy1):
                    for_remove.append(y)
        clean_list = [v for i, v in enumerate(list_to_clean) if i not in set(for_remove)]
        return clean_list
    if boolean(custom):
        pattern_custom = Chem.MolFromSmarts(smart_custom)
        custom = list(mol.GetSubstructMatches(pattern_custom))
        to_del_bef_custom = []

        for x in reversed(range(len(torsion))):
            for y in reversed(range(len(custom))):
                ix1, ix2 = ig(1)(torsion[x]), ig(2)(torsion[x])
                iy1, iy2 = ig(1)(custom[y]), ig(2)(custom[y])
                if (ix1 == iy1 and ix2 == iy2) or (ix1 == iy2 and ix2 == iy1):
                    to_del_bef_custom.append(x)

        custom_torsion = copy(torsion)
        custom_torsion = [v for i, v in enumerate(custom_torsion) if i not in set(to_del_bef_custom)]
        custom_torsion = cleaner(custom_torsion)
    cistrans = cleaner(cistrans)
    torsion = cleaner(torsion)
    mol = Chem.AddHs(mol)
    atoms = mol.GetNumAtoms()
    bonds = mol.GetNumBonds()

    return (atoms, bonds, torsion, cistrans, custom_torsion)


def template_sdf(smile, distance_cutoff_1, distance_cutoff_2):
    """Create a template sdf string and writes it to file.

    Args:
        smile (str): one-line representation of the molecule
        distance_cutoff_1 (float): min distance between non-bonded atoms [A]
        distance_cutoff_2 (float): max distance between bonded atoms [A]
    Returns:
        sdf string
    """
    cnt = 0
    sdf_check = True
    while sdf_check:
        mol = Chem.MolFromSmiles(smile)
        mol = Chem.AddHs(mol)
        AllChem.EmbedMolecule(mol)
        AllChem.UFFOptimizeMolecule(mol)
        Chem.SDWriter('mol.sdf').write(mol)
        sdf_string = Chem.MolToMolBlock(mol)
        check = check_geo_sdf(sdf_string, distance_cutoff_1, distance_cutoff_2)
        if check:
            sdf_check = False
            Chem.SDWriter('mol.sdf').write(mol)
        else:
            cnt += 1
    return sdf_string
