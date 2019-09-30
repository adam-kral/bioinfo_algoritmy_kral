#! /bin/env python3

from collections import OrderedDict
from unittest.mock import patch
from Bio.PDB import PDBParser

# in Python 3.8 there will be cached_property
from functools import lru_cache

from Bio.PDB.Structure import Structure as BioStructure
from Bio.PDB.Model import Model as BioModel
from Bio.PDB.Chain import Chain as BioChain
from Bio.PDB.Residue import Residue as BioResidue


class CountModelsMixin:
    @property
    @lru_cache(maxsize=1)
    def count_models(self):
        return sum(1 for _ in self.get_models())


class CountChainsMixin:
    @property
    @lru_cache(maxsize=1)
    def count_chains(self):
        return sum(1 for _ in self.get_chains())


class CountResiduesMixin:
    @property
    @lru_cache(maxsize=1)
    def count_residues(self):
        return sum(1 for _ in self.get_residues())


class CountAtomsMixin:
    @property
    @lru_cache(maxsize=1)
    def count_atoms(self):
        return sum(1 for _ in self.get_atoms())


class Structure(BioStructure, CountModelsMixin, CountChainsMixin, CountResiduesMixin, CountAtomsMixin):
    pass

class Model(BioModel, CountChainsMixin, CountResiduesMixin, CountAtomsMixin):
    pass

class Chain(BioChain, CountResiduesMixin, CountAtomsMixin):
    pass

class Residue(BioResidue, CountAtomsMixin):
    pass


# have to do it this way, Bio.PDB.StructureBuilder is not extensible at all, cannot provide my own subclasses to it, I would have to copy
# in a lot of its code

patch('Bio.PDB.StructureBuilder.Structure', Structure).__enter__()
patch('Bio.PDB.StructureBuilder.Model', Model).__enter__()
patch('Bio.PDB.StructureBuilder.Chain', Chain).__enter__()
patch('Bio.PDB.StructureBuilder.Residue', Residue).__enter__()


class DistanceAnalyzer:
    # entity can be a (single-model) structure, a model or a chain
    def __init__(self, entity):
        self.structure = entity
        self.atoms = {atom.get_serial_number(): atom for atom in entity.get_atoms()}

    def _atoms_within_d_to_atom(self, hetatm, d):
        for atom in structure.get_atoms():
            candidate_d = atom-hetatm
            if candidate_d <= d:
                yield candidate_d, atom

    # return list of tuples (distance, atom) sorted by distance to the hetero atom
    # complexity: Theta(#structure_atoms)
    # note: includes itself if `d` non-negative
    def atoms_within_d_to_atom(self, atom_serial, d):
        return sorted(list(self._atoms_within_d_to_atom(self.atoms[atom_serial], d)), key=lambda tup: tup[0])

    # return list of tuples (distance, residue) sorted by distance := d(closest residue's atom; hetero atom)
    # complexity: Theta(#structure_atoms)
    # note: includes itself if `d` non-negative
    def residues_with_d_to_atom(self, atom_serial, d):
        residues = OrderedDict()

        for dist, atom in self.atoms_within_d_to_atom(atom_serial, d):
            if atom.parent not in residues:
                residues[atom.parent] = dist

        return [(dist, residue) for residue, dist in residues.items()]


# entity can be a (single-model) structure, a model or a chain
# complexity: Theta(squared(#structure_atoms))
def structure_width(entity):
    atoms = list(entity.get_atoms())

    largest_d = 0

    for i in range(len(atoms)):
        for j in range(i+1, len(atoms)):
            d = atoms[i] - atoms[j]
            if d > largest_d:
                largest_d = d

    return largest_d


# entity can be a (single-model) structure, a model or a chain
def get_hetero_atoms(entity):
    for residue in entity.get_residues():
        # ignore non-hetero atoms and water ('W')
        if not residue.id[0].startswith('H_'):
            continue

        for atom in residue:
            yield atom


if __name__ == '__main__':
    import os
    structure = PDBParser(QUIET=True).get_structure('1tup', os.path.dirname(__file__) + os.path.sep + 'test_data/1tup.pdb')

    print(structure.count_residues)
    print(structure.count_atoms)

    hetero_atoms = list(get_hetero_atoms(structure))

    a = DistanceAnalyzer(structure)

    hetatm_serial = hetero_atoms[0].get_serial_number()

    print(a.atoms_within_d_to_atom(hetatm_serial, 10))
    print(a.residues_with_d_to_atom(hetatm_serial, 10))

    print(structure_width(structure))
    print(structure_width(structure[0]))

    for chain in structure[0]:
        print(chain)
        print(structure_width(chain))

    for model in structure:
        print(model)

        for chain in model:
            print('chain', chain)
            print('len of chain', len(chain))
            for residue in chain:
                if residue.id[0].startswith('H_'):
                    for atom in residue:
                        print(atom.get_name())
                        print(atom.get_serial_number())
                    print(residue.id)
