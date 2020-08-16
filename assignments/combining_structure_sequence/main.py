
# 1) load MSA (get some test data, maybe I already have)
# 2) load structure, get the model
# 3) determine active sites (heteroatoms/dna? within 2.5 A)
# 4) scatterplot (2xcontinuous) conservation score/active-nonactive  or box plot
import itertools

from Bio.Data.SCOPData import protein_letters_3to1
from Bio.SubsMat import MatrixInfo
import seaborn as sns


from ..clustal_parser import parse_clustal
from ..pdb.pdb import *


def aa_sequence_from_pdb_structure(structure):
    aas = filter(lambda r: r.is_aa, structure.get_residues())
    aas_in_order = [None] * len(aas)

    for aa in aas:
        _, seq_idx, _ = aa.id

        try:
            aas_in_order[seq_idx] = aa
        except IndexError:
            raise Exception(f'Inconsistent sequence indices of amino acids in structure, `{seq_idx}` too large')


class MismatchError(Exception):
    pass


def get_ungapped_aa_seq_to_msa_columns(msa_seq):
    # returns mapping from structure's (aa) sequence to column indices in the msa with corresponding aas

    ungapped_seq_index_to_msa_columns = []

    for i in range(len(msa_seq)):
        if msa_seq[i] != '-':  # is not gap
            ungapped_seq_index_to_msa_columns.append(i)

    return ungapped_seq_index_to_msa_columns


def assert_chain_aas_correspond_to_msa_seq_aas(chain_aas, msa_seq, ungapped_seq_index_to_msa_columns):
    for aa in chain_aas:
        chain_aa_code = protein_letters_3to1[aa.resname]
        aa_seq_index = aa.id[1]

        msa_aa_index = ungapped_seq_index_to_msa_columns[aa_seq_index - 1]  # -1 , because string indexing is zero-based, but first aa has
        # SSSEQ=1
        try:
            msa_aa_code = msa_seq[msa_aa_index]
        except IndexError:
            raise MismatchError(f'MSA sequence too short. Lacks amino acid at MSA column {msa_aa_index} (ungapped: {aa.id[1]})')

        if msa_aa_code != chain_aa_code:
            print(MismatchError(f'MSA amino acid `{msa_aa_code}` does not match `{chain_aa_code}`, the amino acid in the structure'
                                f' at MSA column {msa_aa_index} (ungapped: {aa.id[1]})'))


def conservation_rate_of_chain_aas(msa, chain_aas, aa_seq_to_msa_columns):
    for aa in chain_aas:
        yield msa.sum_of_pairs_column(aa_seq_to_msa_columns[aa.id[1]], MatrixInfo.blosum62)


if __name__ == '__main__':
    import os

    # Read and parse MSA.
    with open(os.path.dirname(__file__) + os.path.sep + 'test_data/clustal_assigned') as f:
        msa = parse_clustal(f)

    print(msa)

    model = PDBParser(QUIET=True).get_structure('1tup', os.path.dirname(__file__) + os.path.sep + 'test_data/1tup.pdb')[0]
    chain = model['A']

    # note that the whole structure, not the single chain!!!
    hetero_atoms = list(get_hetero_atoms(model))

    # now I'm interested just
    d_analyzer = DistanceAnalyzer(chain)

    active_site_residues = set()

    for hetatm in hetero_atoms:
        distances_residues_near = d_analyzer.residues_within_d_to_atom(hetatm, 2.5)
        active_site_residues.update(r for d, r in distances_residues_near if r.is_aa)

    non_active_site_residues = list(set(chain.aa_residues).difference(active_site_residues))
    active_site_residues = list(active_site_residues)  # making it list, next code depends on the order


    db_seq_id = 'UniRef90_P04637'
    msa_seq = next(filter(lambda seq_record: db_seq_id == seq_record.id, msa)).seq

    chain_aas = list(chain.aa_residues)

    # někde se to prostě liší -> warningy a to dané místo přeskočit

    ungapped_seq_index_to_msa_columns = get_ungapped_aa_seq_to_msa_columns(msa_seq)
    assert_chain_aas_correspond_to_msa_seq_aas(chain_aas, msa_seq, ungapped_seq_index_to_msa_columns)

    conservation_active = conservation_rate_of_chain_aas(msa, active_site_residues, ungapped_seq_index_to_msa_columns)
    conservation_non_active = conservation_rate_of_chain_aas(msa, non_active_site_residues, ungapped_seq_index_to_msa_columns)

    import matplotlib.pyplot as plt

    ax = sns.swarmplot(
        x=list(itertools.chain((False for _ in range(len(non_active_site_residues))), (True for _ in range(len(active_site_residues))))),
        y=list(itertools.chain(conservation_non_active, conservation_active)),
    )
    ax.set(xlabel='is active site', ylabel='conservation rate (SoP)')
    plt.show()


    # todo more test data (results not great...)
    # did test with 5 angstroms not much better (to incorporate also DNA-binding regions. The zinc binding were only seen when 2.5 A)


    # dbref -- reference to the entry having the sequence
    # 1) resname to one letter
    # 2) full_id[2][1]

    # map alignment (may have insertions! to sequence ids!!)

    # option to map sequence to structure (just ignore leading or trailing aas, but the sequence has to be contiguous)