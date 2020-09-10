import itertools

from Bio import pairwise2
from Bio.Data.SCOPData import protein_letters_3to1
from Bio.SubsMat import MatrixInfo
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

    # def __init__(self, match_fn, gap_fn):
    #     self.match_fn = match_fn
    #     self.gap_fn = gap_fn

    def match(self, i, j, paths):
        a = self.s1[i]
        b = self.s2[j]

        try:
            return MatrixInfo.blosum62[a, b]
        except KeyError:
            return MatrixInfo.blosum62[b, a]

    # parameter `gap_position_in_sequence` is before which letter the gap is introduced (zero-based). E.g. 0 means leading gap (before 0th
    # letter)
    def _gap(self, i, j, paths, sequence_length, gap_position_in_sequence, gap_ext_flag):
        # assuming gap extension penalties <= gap open penalties ( >= for `scores`)

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
    def __init__(self, aa_feature_dict):
        sorted_features = sorted(aa_feature_dict.items(), key=lambda aa_and_value: aa_and_value[1])

        scale_max_distance_to_one = len(sorted_features) - 1

        self.distance_index = {}
        for i, (aa, value) in enumerate(sorted_features):
            self.distance_index[aa] = i / scale_max_distance_to_one

    def __call__(self, aa1, aa2):
        return abs(self.distance_index[aa1] - self.distance_index[aa2])


# similar structural neighborhood - better match (higher match score - additive?, worse gap score - multiplicative (1-2)x)
# asi bych se musel podívat na rozdělení structural neighborhood
# na gap se vyprdnu, tam stejne nevim, co bych hodnotil. Jakože jsou jiný ty insertlý aa, než ty krajní na gapu? Asi ne...
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


def test():
    min_neighbor_res_count_metrics = StructuralFeatureMetric(residues_get_min_neighbor_residues_count(model.aa_residues, 8.5))

    computer = StructuralFeaturesSequenceAlignmentComputer()



# todo taking into account also structural features of the aligned residues. These can be, for example, number of residues in a sphere
#  centered in given residuum. Sequence information can be used as well.
# takže zkusim samotny + kombinaci
# rozmyslet, jak to napisu a jaky vsechny featury pouziju
# jak otestuji kvalitu? Pomocí RMSD

# todo test data


# in case of a protein complex, run this with each of the selected polypeptide pairs. And then run the RMSD minimization with all the
# point pairs combined
def get_c_alpha_coords_matched_residues(alignment, pp1, pp2):
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


# returns centroid of P, Q and a rotational matrix (in subsequent analyses, you may want to subtract centroidP/Q from your point sets)
# that minimizes RMSD between corresponding points in P and Q, when applied to P
def center_and_kabsch(P, Q):
    P = np.array(P)
    Q = np.array(Q)

    centroidP = rmsd.centroid(P)
    centroidQ = rmsd.centroid(Q)

    P = P - centroidP
    Q = Q - centroidQ

    return centroidP, centroidQ, rmsd.kabsch(P, Q)



class Superposition:
    # add also the alignment?
    def __init__(self, superposition_structure, c_alpha_points_1, c_alpha_points_2):
        self.structure = superposition_structure
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


def superimpose_chains(chain1, chain2, pp1, pp2, alignment):
    """ chains and their corresponding polypeptides (polypeptide must be a subset of chain's residues, of identical objects in memory) """
    P, Q = get_c_alpha_coords_matched_residues(alignment, pp1, pp2)

    assert len(P) == len(Q) > 0

    # run Kabsch algorithm (minimizing RMSD=root-mean-square deviation)
    centroidP, centroidQ, rotation_matrix = center_and_kabsch(P, Q)

    for atom in chain1.get_atoms():
        atom.coord = np.dot(rotation_matrix, atom.coord - centroidP) + centroidQ

    chain1.id = 'U'
    chain2.id = 'V'

    from ..pdb.pdb import Structure, Model

    structure = Structure('superposition')
    model = Model(0)
    structure.add(model)

    model.add(chain1)
    model.add(chain2)

    # run `get_c_alpha_coords_matched_residues` again to get coords of *kabsch-transformed* c_alphas
    return Superposition(structure, *get_c_alpha_coords_matched_residues(alignment, pp1, pp2))


def save_structure(structure, file_name):
    io = PDBIO()
    io.set_structure(structure)
    io.save(file_name)


def get_superposition_rmsd_sequence(superposition, min_c_alphas=1):
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


if __name__ == '__main__':
    import os

    from Bio.PDB import PDBParser, PPBuilder
    from Bio.PDB.PDBIO import PDBIO

    parser = PDBParser()
    chain1 = parser.get_structure('a2a', os.path.dirname(__file__) + os.path.sep + 'test_data/A2a_receptor.pdb')[0]['A']
    chain2 = parser.get_structure('3eml', os.path.dirname(__file__) + os.path.sep + 'test_data/3eml.pdb')[0]['A']

    superposition_rmsds = []

    def superposition_with_blosum(chain1, chain2, alignment):
        pp1, pp2 = map(chain_to_one_pp, (chain1, chain2))

        superposition = superimpose_chains(chain1, chain2, pp1, pp2, alignment)

        save_structure(superposition.structure, 'blosum_alignment.pdb')
        print(superposition.calculate_rmsd(), superposition.matched_residues_count)

        superposition_rmsds.append(('blosum_al', list(get_superposition_rmsd_sequence(superposition))))

    pp1, pp2 = map(chain_to_one_pp, (chain1, chain2))
    als = Blosum62AlignmentComputer(pp1.get_sequence(), pp2.get_sequence()).compute()

    # there might be thousand of optimal alignments, take only one (the next ones could be very similar, so it doesn't make sense to
    # include just n first, not representative anyway)
    al = next(als)

    superposition_with_blosum(chain1.copy(), chain2.copy(), al)

    def superposition_with_blosum_and_structural_features(chain1, chain2):
        polypeptide1, polypeptide2 = map(chain_to_one_pp, (chain1, chain2))

        residues_min_neighbor_count = residues_get_min_neighbor_atoms_count(polypeptide1, 8.5)
        residues_min_neighbor_count.update(residues_get_min_neighbor_atoms_count(polypeptide2, 8.5))

        min_neighbor_res_count_metrics = StructuralFeatureMetric(residues_min_neighbor_count)

        for structural_score_weight in (0.5, 1, 2, 3, 4, 10, 2000):
            print('processing ', structural_score_weight)

            als = StructuralFeaturesSequenceAlignmentComputer(polypeptide1.get_sequence(),
                                                              polypeptide2.get_sequence(),
                                                              polypeptide1,
                                                              polypeptide2,
                                                              [(min_neighbor_res_count_metrics, structural_score_weight)],
                                                              1).compute()
            als = list(als)
            print('no of als: ', len(als))
            alignment = als[0]

            ch1, ch2 = chain1.copy(), chain2.copy()  # deep copying the chains, because we will alter the coordinates to form the
            # superposition
            pp1, pp2 = map(chain_to_one_pp, (ch1, ch2))  # chain and polypeptide have to be formed with the same Residue objects (in
            # memory) for superimpose_chains to work

            superposition = superimpose_chains(ch1, ch2, pp1, pp2, alignment)

            save_structure(superposition.structure, f'blosum_alignment_structural_features_{structural_score_weight}.pdb')
            print(superposition.calculate_rmsd(), superposition.matched_residues_count)

            superposition_rmsds.append((f'blosum_struct_al{structural_score_weight}',
                                        list(get_superposition_rmsd_sequence(superposition))))


    superposition_with_blosum_and_structural_features(chain1.copy(), chain2.copy())

    import seaborn as sns
    from matplotlib import pyplot as plt

    print(np.transpose(list(itertools.chain(*(enumerate(values, start=1) for name, values in superposition_rmsds)))))
    x, y = np.transpose(list(itertools.chain(*(enumerate(values, start=1) for name, values in superposition_rmsds))))

    sns.lineplot(x=x, y=y, hue=[name for name, values in superposition_rmsds for _ in values])
    plt.show()

# todo otestovat i s ruznymi gap penaltami?
# todo ukládat do separate filu pdb (1) je naming v pymolu jasny (2) kdyz budu superimposovat vic chainu, nebude v tom bordel,
#  bude to konzistentni




# structure scoring. Použít, co už mám v
# vyhodnocení? RMSD, ale pokazdy je jinej počet párů, tak jak? Normalizace na počet atomů celý struktury? Ne to nejde
# ale šlo by, že páry insercí by byly ke stejným krajním aa, mezi kterými jsou gapy?

# co mě zajímá: m = #matched residues, R=RMSD. pokud m1>~m2, pak je lepší 1. superposition <=> R1 < R27
# vycházím z toho, že alignment se nezmění, je to součást té superposition. Takže v situaci m1< m2 můžu buď přidávat zbytky k 1 (s
# nejmenším RMSD nebo od 2 zbytky s největším RMSD odebrat. Problém, které přidávat (ke krajním aa gapů, nebo i jinam - složitější).
# Takže spíš odebírat od 2 ty nejhorší (a třeba nastane R2>R1)? Co když ale odebráním ještě jednoho nejhoršího v 1 najednou zas bylo R1 <
# R2? Můžu udělat graf s postupným odebírání z obou superpozicí... A ten pak vyhodnotit nějak. Např. pokud bude stabilně nižší RMSD u
# ekvivalentních délek je to prostě lepší superposition
# dobře, tak to bude evaluace.

# jak udělat structural features sequence alignment?
# najít rozdělení random podobných aa featur (1d), vzdálenost - počet aa mezi nimi. To počtem normalizuju na jedničku a pak můžu použít
# aditivně/multiplikativně,  (Jakkoliv škálované) pro jakoukoliv feature. Někde jsou ale aa hojné a bude to dávat šum a horší výsledek?
# Šlo by to nějak jinak, pravděpodobnostně?
# pravděpodobnost, že by
