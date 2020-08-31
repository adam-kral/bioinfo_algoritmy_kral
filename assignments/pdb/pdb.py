#! /bin/env python3

from collections import OrderedDict
from unittest.mock import patch
from Bio.PDB import PDBParser, is_aa

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
    @property
    def aa_residues(self):
        for chain in self:
            yield from chain.aa_residues

class Chain(BioChain, CountResiduesMixin, CountAtomsMixin):
    @property
    def aa_residues(self):
        return filter(lambda r: r.is_aa, self.get_residues())

class Residue(BioResidue, CountAtomsMixin):
    def __repr__(self):
        # based on original
        """Return the residue full id."""
        resname = self.get_resname()
        hetflag, resseq, icode = self.get_id()

        chainstr = f'ch={self.get_parent().id} ' if self.get_parent() else ''

        return f'<Residue {resname} het={hetflag} {chainstr}resseq={resseq} icode={icode}>'

    @property
    def hetatm_flag(self):
        return self.id[0]

    @property
    def is_water(self):
        return self.hetatm_flag == 'W'

    @property
    def is_hetatm_excl_w(self):
        return self.hetatm_flag.starts_with('H_')

    @property
    def is_aa(self):
        return is_aa(self)
        # return not self.hetatm_flag.strip() and 'CA' in (a.get_name() for a in self.get_atoms())




# have to do it this way, Bio.PDB.StructureBuilder is not extensible at all, cannot provide my own subclasses to it, I would have to copy
# in a lot of its code
# patch is used as context manager, therefore calling __enter__
patch('Bio.PDB.StructureBuilder.Structure', Structure).__enter__()
patch('Bio.PDB.StructureBuilder.Model', Model).__enter__()
patch('Bio.PDB.StructureBuilder.Chain', Chain).__enter__()
patch('Bio.PDB.StructureBuilder.Residue', Residue).__enter__()


class DistanceAnalyzer:
    # entity can be a (single-model) structure, a model or a chain
    def __init__(self, entity):
        self.structure = entity

    def _atoms_within_d_to_atom(self, hetatm, d):
        for atom in self.structure.get_atoms():
            actual_d = atom-hetatm
            if actual_d <= d:
                yield actual_d, atom

    # return list of tuples (distance, atom) sorted by distance to the hetero atom
    # complexity: Theta(#structure_atoms) todo určitě ne!
    # note: includes itself if `d` non-negative
    def atoms_within_d_to_atom(self, atom_of_interest, d):
        return sorted(list(self._atoms_within_d_to_atom(atom_of_interest, d)), key=lambda tup: tup[0])

    # return list of tuples (distance, residue) sorted by distance := d(closest residue's atom; hetero atom)
    # complexity: Theta(#structure_atoms)
    # note: includes itself if `d` non-negative
    def residues_within_d_to_atom(self, atom_of_interest, d):
        residues = OrderedDict()

        for dist, atom in self.atoms_within_d_to_atom(atom_of_interest, d):
            # set lowest distance yet
            if atom.parent not in residues or residues[atom.parent] > dist:
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
        # # ignore non-hetero atoms and water ('W')
        # if not residue.id[0].startswith('H_'):
        #     continue
        # # alternatively residue.id[0] != '' (that's the hetero field, but could be also water as in their comments)

        if residue.id[0] == 'W' or residue.is_aa:
            continue

        for atom in residue:
            yield atom


if __name__ == '__main__':
    import os
    structure = PDBParser(QUIET=True).get_structure('1tup', os.path.dirname(__file__) + os.path.sep + 'test_data/1tup.pdb')

    print(structure.count_residues)
    print(structure.count_atoms)

    hetero_atoms = list(get_hetero_atoms(structure))
    #
    # a = DistanceAnalyzer(structure)
    #
    # hetatm = hetero_atoms[0]
    #
    # print(a.atoms_within_d_to_atom(hetatm, 10))
    # print(a.residues_within_d_to_atom(hetatm, 10))

    # print(structure_width(structure))
    # print(structure_width(structure[0]))

    for chain in structure[0]:
        print(chain)
        # print(structure_width(chain))

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

