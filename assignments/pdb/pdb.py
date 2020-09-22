from unittest.mock import patch
from Bio.PDB import PDBParser, is_aa

from Bio.PDB.Structure import Structure as BioStructure
from Bio.PDB.Model import Model as BioModel
from Bio.PDB.Chain import Chain as BioChain
from Bio.PDB.Residue import Residue as BioResidue


class CountModelsMixin:
    @property
    def count_models(self):
        return sum(1 for _ in self.get_models())


class CountChainsMixin:
    @property
    def count_chains(self):
        return sum(1 for _ in self.get_chains())


class CountResiduesMixin:
    @property
    def count_residues(self):
        return sum(1 for _ in self.get_residues())


class CountAtomsMixin:
    @property
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
        # customized BioPython's code to also show the chain name
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


def attach_custom_classes_to_pdb_structure_builder():
    """ have to do it this way, Bio.PDB.StructureBuilder is not extensible at all, cannot provide my own subclasses to it, I would have
    to copy in a lot of its code.
    From BioPython code comments: "The StructureBuilder class is used by the PDBParser classes to translate a file to a Structure object."
    """
    # Patch is used as a context manager, therefore calling __enter__
    patch('Bio.PDB.StructureBuilder.Structure', Structure).__enter__()
    patch('Bio.PDB.StructureBuilder.Model', Model).__enter__()
    patch('Bio.PDB.StructureBuilder.Chain', Chain).__enter__()
    patch('Bio.PDB.StructureBuilder.Residue', Residue).__enter__()


def _atoms_within_d_to_atom(entity, atom_of_interest, d):
    for atom in entity.get_atoms():
        actual_d = atom-atom_of_interest
        if actual_d <= d:
            yield atom, actual_d


def atoms_within_d_to_atom(entity, atom_of_interest, d):
    """ Returns list of tuples (atom, distance) sorted by distance := Euclidean_distance(atom_of_interest; entity's atom)

    :param entity: contains superset of the returned atoms
    :param atom_of_interest: distance of entity's atoms is measured to this atom
    complexity: O(n*logn); n := #entity atoms
    """
    return sorted(list(_atoms_within_d_to_atom(entity, atom_of_interest, d)), key=lambda tup: tup[1])


def residues_within_d_to_atom(entity, atom_of_interest, d):
    """ Returns list of tuples (residue, distance) sorted by distance := Euclidean_distance(atom_of_interest; closest residue's atom)
    complexity: O(n*logn); n := #entity atoms
    """
    residues = {}

    for atom, dist in atoms_within_d_to_atom(entity, atom_of_interest, d):
        # set lowest distance yet
        if atom.parent not in residues or residues[atom.parent] > dist:
            residues[atom.parent] = dist

    return sorted(residues.items(), key=lambda t: t[1])


def structure_width(entity):
    """ Returns maximum distance between entity's atoms

    :param entity: a (single-model) structure, a model or a chain
    complexity: Theta(squared(#structure_atoms))
    """
    atoms = list(entity.get_atoms())

    largest_d = 0

    for i in range(len(atoms)):
        for j in range(i+1, len(atoms)):
            d = atoms[i] - atoms[j]
            if d > largest_d:
                largest_d = d

    return largest_d


def get_hetero_atoms(entity):
    """ Returns all non-water and non-amino-acid atoms. (E.g. zinc HETATM or nucleotide's atoms)

    :param entity: a (single-model) structure, a model or a chain
    """
    for residue in entity.get_residues():
        # ignore non-hetero atoms and water ('W')

        # if not residue.id[0].startswith('H_'):
        #     continue
        # # alternatively residue.id[0] != '' (that's the hetero field, but could be also water as in their comments)

        if residue.id[0] == 'W' or residue.is_aa:
            continue

        for atom in residue:
            yield atom


if __name__ == '__main__':
    import os

    attach_custom_classes_to_pdb_structure_builder()
    structure = PDBParser(QUIET=True).get_structure('1tup', os.path.dirname(__file__) + os.path.sep + 'test_data/1tup.pdb')

    print(structure.count_residues)
    print(structure.count_atoms)

    hetero_atoms = list(get_hetero_atoms(structure))

    hetatm = hetero_atoms[0]

    print(atoms_within_d_to_atom(structure, hetatm, 10))
    print(residues_within_d_to_atom(structure, hetatm, 10))

    # print(structure_width(structure))

    for chain in structure[0]:
        print(chain)
        # print(structure_width(chain))

    for model in structure:
        print(model)

        for chain in model:
            print(chain)
            print('number of residues in chain:', len(chain))
            for residue in chain:
                if residue.id[0].startswith('H_'):
                    for atom in residue:
                        print(atom.get_name())
                        print(atom.get_serial_number())
                    print(residue.id)

