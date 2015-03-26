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
"""Measure and set torsions values."""

from operator import itemgetter
from rdkit import Chem
from rdkit.Chem import rdMolTransforms


def ig(x):
    return itemgetter(x)


def angle_measure(sdf_string, cistrans, torsion):
    """Return the dihedral angles (in degrees).

    Args:
        sdf_string
        cistrans (list of tuples): positions of cistrans bonds
        torsion (list of tuples): positions of rotatable bonds
    Returns:
        list of values of cistrans bonds and list of values of rotatable
        bonds in deg
    Raises:
        ValueError: if both lists are empty

    """
    if len(cistrans) == 0 and len(torsion) == 0:
        raise ValueError("There are no rotatable bonds and no cis/trans bonds")
    mol = Chem.MolFromMolBlock(sdf_string, removeHs=False)
    measured_torsion = []
    measured_cistrans = []
    if len(cistrans) != 0:
        for x in range(len(cistrans)):
            tmp = float(rdMolTransforms.GetDihedralDeg(
                mol.GetConformer(),
                ig(0)(cistrans[x]), ig(1)(cistrans[x]),
                ig(2)(cistrans[x]), ig(3)(cistrans[x])))
            measured_cistrans.append(float('{0:.2f}'.format(tmp)))
    if len(torsion) != 0:
        for y in range(len(torsion)):
            tmp = float(rdMolTransforms.GetDihedralDeg(
                mol.GetConformer(),
                ig(0)(torsion[y]), ig(1)(torsion[y]),
                ig(2)(torsion[y]), ig(3)(torsion[y])))
            measured_torsion.append(float('{0:.2f}'.format(tmp)))

    return measured_cistrans, measured_torsion


def angle_set(sdf_string, cistrans, torsion, values_cistrans_n, values_tor_n):
    """Return a modifies sdf_string after setting torsions.

    Args:
        sdf_string
        cistrans (list of tuples): positions of cistrans bonds
        torsion (list of tuples): positions of rotatable bonds
        values_cistrans_n: cistrans bonds values to set (in degrees)
        values_tor_n: rotatable bonds values to set (in degrees)
    Returns:
        modified sdf_string
    Raises:
        ValueError: if both lists of positions are empty
        ValueError: if lengths of the corresponding lists do not match

    """
    if len(cistrans) == 0 and len(torsion) == 0:
        raise ValueError("There are no rotatable bonds and no cis/trans bonds")
    if len(cistrans) != len(
                    values_cistrans_n) or len(torsion) != len(values_tor_n):
        raise ValueError("No length match between the values and positions")
    mol = Chem.MolFromMolBlock(sdf_string, removeHs=False)
    if len(torsion) != 0:
        for x in range(len(torsion)):
            rdMolTransforms.SetDihedralDeg(
                mol.GetConformer(), ig(0)(torsion[x]),
                ig(1)(torsion[x]), ig(2)(torsion[x]),
                ig(3)(torsion[x]), ig(x)(values_tor_n))
    if len(cistrans) != 0:
        for x in range(len(cistrans)):
            rdMolTransforms.SetDihedralDeg(
                mol.GetConformer(), ig(0)(cistrans[x]),
                ig(1)(cistrans[x]), ig(2)(cistrans[x]),
                ig(3)(cistrans[x]), ig(x)(values_cistrans_n))

    return Chem.MolToMolBlock(mol)
