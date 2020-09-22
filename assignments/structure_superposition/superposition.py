import itertools

from Bio.SubsMat import MatrixInfo
from Bio.PDB import PPBuilder
from Bio.PDB.PDBIO import PDBIO

import rmsd
import numpy as np

from assignments.structure_related_properties.lib import residues_get_min_neighbor_residues_count, residues_get_neighboring_water_count, \
    residues_get_min_neighbor_atoms_count
from ..alignment import AlignmentComputer

# cannot use Biopython pairwise2.align.globalcc(), because the match callback is called only with characters, not their positions
# -> will modify my EditDistance homework and use it

# looking at BioPython, I also could use Superimposer, but that can transform just the matched atoms (I match only c_alphas
# and then transform _all_ atoms in the chain)


class Blosum62AlignmentComputer(AlignmentComputer):
    # gap scores for alignment - https://www.drive5.com/usearch/manual/aln_params.html
    END_GAP_OPEN = -1
    END_GAP_EXT = -0.5

    GAP_OPEN = -17
    GAP_EXT = -1

    def match(self, i, j, paths):
        a = self.s1[i]
        b = self.s2[j]

        try:
            return MatrixInfo.blosum62[a, b]
        except KeyError:
            return MatrixInfo.blosum62[b, a]

    def _gap(self, i, j, paths, sequence_length, gap_position_in_sequence, gap_ext_flag):
        """  parameter `gap_position_in_sequence` is before which letter the gap is introduced (zero-based). E.g. 0 means leading gap
        (before 0th residue/letter)

        assuming gap extension penalties <= gap open penalties ( >= for `scores`)
        """

        optimal_step_to_previous_cell = paths[i, j]  # step to matrix[i,j]

        # leading or trailing gap
        if gap_position_in_sequence in (0, sequence_length):
            if optimal_step_to_previous_cell & gap_ext_flag:
                return self.END_GAP_EXT

            return self.END_GAP_OPEN

        # gap somewhere in the sequence
        if optimal_step_to_previous_cell & gap_ext_flag:
            return self.GAP_EXT

        return self.GAP_OPEN

    def gap_s1(self, i, j, paths):
        return self._gap(i, j, paths, len(self.s1), i, self.LEFT_FLAG)

    def gap_s2(self, i, j, paths):
        return self._gap(i, j, paths, len(self.s2), j, self.TOP_FLAG)


class StructuralFeatureMetric:
    """ Simple metric, initialized with values for each residue. Upon call with two amino acids returns ratio of residues that have
    values in between the two, if they were sorted. Thus the metric returns values in [0,1] interval. """

    def __init__(self, aa_feature_dict):
        sorted_features = sorted(aa_feature_dict.items(), key=lambda aa_and_value: aa_and_value[1])

        scale_max_distance_to_one = len(sorted_features) - 1

        self.distance_index = {}
        for i, (aa, value) in enumerate(sorted_features):
            self.distance_index[aa] = i / scale_max_distance_to_one

    def __call__(self, aa1, aa2):
        return abs(self.distance_index[aa1] - self.distance_index[aa2])


class StructuralFeaturesSequenceAlignmentComputer(Blosum62AlignmentComputer):
    def __init__(self, s1, s2, pp1, pp2, additive_scoring_metrics_with_weights, gap_weight):
        super().__init__(s1, s2)

        self.pp1 = pp1
        self.pp2 = pp2
        self.additive_scoring_metrics_with_weights = additive_scoring_metrics_with_weights
        self.gap_weight = gap_weight

    def match(self, i, j, paths):
        blosum62_score = super().match(i, j, paths)

        structural_score = 0
        for d_fn, weight in self.additive_scoring_metrics_with_weights:
            structural_score += weight * (1 - d_fn(self.pp1[i], self.pp2[j]))

        return blosum62_score + structural_score

    def _gap(self, i, j, paths, sequence_length, gap_position_in_sequence, gap_ext_flag):
        return self.gap_weight * super()._gap(i, j, paths, sequence_length, gap_position_in_sequence, gap_ext_flag)


def get_c_alpha_coords_matched_residues(alignment, pp1, pp2):
    """ in case of a protein complex, run this with each of the selected polypeptide pairs. And then run the RMSD minimization with all the
    point pairs combined
    """

    P = []
    Q = []

    i = 0
    j = 0

    for a, b in zip(*alignment.alignment_tuple):
        if '-' not in (a, b):
            P.append(pp1[i]["CA"].coord)
            Q.append(pp2[j]["CA"].coord)

        if a != '-':
            i += 1

        if b != '-':
            j += 1

    return P, Q


def center_and_kabsch(P, Q):
    """ returns centroid of P, Q and a rotational matrix (in subsequent analyses, you may want to subtract centroidP/Q from your point sets)
    that minimizes RMSD between corresponding points in P and Q, *when applied to P* """

    P = np.array(P)
    Q = np.array(Q)

    centroidP = rmsd.centroid(P)
    centroidQ = rmsd.centroid(Q)

    P = P - centroidP
    Q = Q - centroidQ

    return centroidP, centroidQ, rmsd.kabsch(P, Q)


class SuperpositionCAlphas:
    def __init__(self, c_alpha_points_1, c_alpha_points_2):
        self.c_alpha_points_1 = c_alpha_points_1
        self.c_alpha_points_2 = c_alpha_points_2

        assert len(c_alpha_points_1) == len(c_alpha_points_2)

    def calculate_rmsd(self):
        return rmsd.rmsd(self.c_alpha_points_1, self.c_alpha_points_2)

    @property
    def matched_residues_count(self):
        return len(self.c_alpha_points_1)


def chain_to_one_pp(chain):
    ppb = PPBuilder()

    polypeptides = ppb.build_peptides(chain)

    if len(polypeptides) != 1:
        print('warning ', len(polypeptides), ' polypeptides from one chain, extending first pp')

        for pp in polypeptides[1:]:
            polypeptides[0].extend(pp)

    return polypeptides[0]


def superimpose_chain(chain1, pp1, pp2, alignment):
    """ Superimposes chain1 on pp2

    chain1 and pp1 have to be made of residues that are the same objects in memory for the last line to work
    """

    P, Q = get_c_alpha_coords_matched_residues(alignment, pp1, pp2)

    assert len(P) == len(Q) > 0

    # run Kabsch algorithm (minimizing RMSD=root-mean-square deviation of atom positions)
    centroidP, centroidQ, rotation_matrix = center_and_kabsch(P, Q)

    for atom in chain1.get_atoms():
        atom.coord = np.dot(rotation_matrix, atom.coord - centroidP) + centroidQ
        # side effect to pp1, coordinates will change there too, as the object has identical Atom objects

    # run `get_c_alpha_coords_matched_residues` again to get coords of *kabsch-transformed* c_alphas
    # chain1 coords have changed. Therefore, also residues in pp1 have the coords changed (share same objects)
    # maybe better to change this design
    return SuperpositionCAlphas(*get_c_alpha_coords_matched_residues(alignment, pp1, pp2))


def save_structure(structure, file_name):
    io = PDBIO()
    io.set_structure(structure)
    io.save(file_name)


def get_rmsd_sequence_of_superposition(superposition, min_c_alphas=1):
    """ returns best RMSD with n c_alpha pairs, where n in [min_c_alphas, all superposition c_alphas] """
    V = superposition.c_alpha_points_1
    W = superposition.c_alpha_points_2

    point_pairs = zip(V, W)

    pair_sds = []  # square deviations of atom pairs

    for v,w in point_pairs:
        pair_sds.append(sum((v[i] - w[i]) ** 2.0 for i in range(3)))  # assuming 3 dimensions

    pair_sds.sort()

    pair_sds = iter(pair_sds)

    sd = sum(pair_sd for pair_sd in itertools.islice(pair_sds, min_c_alphas - 1))

    for N, pair_sd in enumerate(pair_sds, start=min_c_alphas):
        sd += pair_sd

        yield np.sqrt(sd / N)


def single_chain_structure(chain, name='superposition'):
    from Bio.PDB.Structure import Structure
    from Bio.PDB.Model import Model

    structure = Structure(name)
    model = Model(0)
    structure.add(model)

    model.add(chain)

    return structure


if __name__ == '__main__':
    import os

    from Bio.PDB import PDBParser, PPBuilder
    from Bio.PDB.PDBIO import PDBIO

    parser = PDBParser()
    chain1 = parser.get_structure('a2a', os.path.dirname(__file__) + os.path.sep + 'test_data/A2a_receptor.pdb')[0]['A']
    chain2 = parser.get_structure('3eml', os.path.dirname(__file__) + os.path.sep + 'test_data/3eml.pdb')[0]['A']

    superposition_rmsds = []

    def superposition_with_blosum(chain1, chain2, alignment):
        chain1 = chain1.copy()  # will change atom coordinates, so create a new object
        pp1, pp2 = map(chain_to_one_pp, (chain1, chain2))

        superposition_calphas = superimpose_chain(chain1, pp1, pp2, alignment)

        save_structure(single_chain_structure(chain1), 'a2a_blosum_alignment.pdb')
        print(superposition_calphas.calculate_rmsd(), superposition_calphas.matched_residues_count)

        superposition_rmsds.append(('blosum_al', list(get_rmsd_sequence_of_superposition(superposition_calphas))))

    pp1, pp2 = map(chain_to_one_pp, (chain1, chain2))
    als = Blosum62AlignmentComputer(pp1.get_sequence(), pp2.get_sequence()).compute()

    # there might be thousand of optimal alignments, take only one (the next ones could be very similar, so it doesn't make sense to
    # include just n first, not representative anyway)
    al = next(als)

    superposition_with_blosum(chain1, chain2, al)

    def superposition_with_blosum_and_structural_features(chain1, chain2):
        polypeptide1, polypeptide2 = map(chain_to_one_pp, (chain1, chain2))

        residues_min_neighbor_count = residues_get_min_neighbor_atoms_count(polypeptide1, 8.5)
        residues_min_neighbor_count.update(residues_get_min_neighbor_atoms_count(polypeptide2, 8.5))

        min_neighbor_res_count_metric = StructuralFeatureMetric(residues_min_neighbor_count)

        for structural_score_weight in (0.5, 1, 2, 3, 4, 10, 2000):
            print('processing ', structural_score_weight)

            als = StructuralFeaturesSequenceAlignmentComputer(polypeptide1.get_sequence(),
                                                              polypeptide2.get_sequence(),
                                                              polypeptide1,
                                                              polypeptide2,
                                                              [(min_neighbor_res_count_metric, structural_score_weight)],
                                                              1).compute()
            als = list(als)
            print('no of als: ', len(als))
            alignment = als[0]

            ch1 = chain1.copy()  # deep copying  chain1, because we will alter the coordinates to form the superposition
            pp1, pp2 = map(chain_to_one_pp, (ch1, chain2))

            superposition_calphas = superimpose_chain(ch1, pp1, pp2, alignment)

            save_structure(single_chain_structure(ch1), f'a2a_blosum_alignment_structural_features_{structural_score_weight}.pdb')
            print(superposition_calphas.calculate_rmsd(), superposition_calphas.matched_residues_count)

            superposition_rmsds.append((f'blosum_struct_al{structural_score_weight}',
                                        list(get_rmsd_sequence_of_superposition(superposition_calphas))))


    superposition_with_blosum_and_structural_features(chain1, chain2)

    import seaborn as sns
    from matplotlib import pyplot as plt

    x, y = np.transpose(list(itertools.chain(*(enumerate(values, start=1) for name, values in superposition_rmsds))))

    sns.lineplot(x=x, y=y, hue=[name for name, values in superposition_rmsds for _ in values])
    plt.show()
