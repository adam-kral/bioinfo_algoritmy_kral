# 1) load MSA (get some test data, maybe I already have)
# 2) load structure, get the model
# 3) determine active sites (heteroatoms/dna? within 2.5 A)
# 4) scatterplot (label x continuous) conservation score/active-nonactive  or box plot
import itertools

from Bio.Data.SCOPData import protein_letters_3to1
from Bio.SubsMat import MatrixInfo
import seaborn as sns


from ..clustal_parser import parse_clustal
from ..pdb.pdb import *


class MismatchError(Exception):
    pass


def get_resseq_to_msa_columns(msa_seq):
    """ Returns mapping from structure's (aa) sequence to column indices in the msa with corresponding aas.
    
    Skips gaps in MSA.
    disclaimer - structure should have atoms without insertion codes (sequential RESSEQ), to correspond to MSA seq
    """

    seq_index_to_msa_columns = [None]  # dummy element because python indexing is zero-based, but first aa has SSSEQ=1,
    # so None as a placeholder

    for i in range(len(msa_seq)):
        if msa_seq[i] != '-':  # is not gap
            seq_index_to_msa_columns.append(i)

    return seq_index_to_msa_columns


def assert_chain_aas_correspond_to_msa_seq_aas(chain_aas, msa_seq, resseq_to_msa_columns):
    for aa in chain_aas:
        chain_aa_code = protein_letters_3to1[aa.resname]
        aa_seq_index = aa.id[1]  # residue sequence number is, in biopython, an int (wonder how they handle '31A' situation in a pdb file)

        msa_aa_index = resseq_to_msa_columns[aa_seq_index]
        try:
            msa_aa_code = msa_seq[msa_aa_index]
        except IndexError:
            raise MismatchError(f'MSA sequence too short. Lacks amino acid at MSA column {msa_aa_index} (ungapped: {aa.id[1]})')

        if msa_aa_code != chain_aa_code:
            print(MismatchError(f'MSA amino acid `{msa_aa_code}` does not match `{chain_aa_code}`, the amino acid in the structure'
                                f' at MSA column {msa_aa_index} (ungapped: {aa.id[1]})'))


def conservation_rate_of_chain_aa(msa, aa, resseq_to_msa_columns):
    return msa.sum_of_pairs_column(resseq_to_msa_columns[aa.id[1]], MatrixInfo.blosum62)


if __name__ == '__main__':
    import os

    # Read and parse MSA.
    with open(os.path.dirname(__file__) + os.path.sep + 'test_data/clustal_assigned') as f:
        msa = parse_clustal(f)

    print(msa)

    # load structure, get the model
    attach_custom_classes_to_pdb_structure_builder()
    model = PDBParser(QUIET=True).get_structure('1tup', os.path.dirname(__file__) + os.path.sep + 'test_data/1tup.pdb')[0]
    chain = model['B']

    # note that the whole structure, not the single chain! In case we want to include nucleotides (another chain)
    hetero_atoms = list(get_hetero_atoms(model))

    active_site_residues = set()

    for hetatm in hetero_atoms:
        residues_and_distances = residues_within_d_to_atom(chain, hetatm, 2.5)  # 4 for dna binding  (donor to acceptor --
        # hydrogen non existent in X-ray structures)

        active_site_residues.update(r for r, d in residues_and_distances if r.is_aa)

    non_active_site_residues = list(set(chain.aa_residues).difference(active_site_residues))
    active_site_residues = list(active_site_residues)  # making it list, next code depends on the order

    msa_seq = next(filter(lambda seq_record: 'UniRef90_P04637' == seq_record.id, msa)).seq

    chain_aas = list(chain.aa_residues)

    resseq_to_msa_columns = get_resseq_to_msa_columns(msa_seq)
    assert_chain_aas_correspond_to_msa_seq_aas(chain_aas, msa_seq, resseq_to_msa_columns)

    conservation_active = [conservation_rate_of_chain_aa(msa, aa, resseq_to_msa_columns)
                           for aa in active_site_residues]

    conservation_non_active = [conservation_rate_of_chain_aa(msa, aa, resseq_to_msa_columns)
                               for aa in non_active_site_residues]

    import matplotlib.pyplot as plt

    print('avg conservation active: ', sum(conservation_active)/len(conservation_active))
    print('avg conservation non-active: ', sum(conservation_non_active)/len(conservation_non_active))

    ax = sns.swarmplot(
        x=list(itertools.chain(itertools.repeat(False, len(non_active_site_residues)), itertools.repeat(True, len(active_site_residues)))),
        y=list(itertools.chain(conservation_non_active, conservation_active)),
    )
    ax.set(xlabel='is active site', ylabel='conservation rate (SoP)')
    plt.show()
