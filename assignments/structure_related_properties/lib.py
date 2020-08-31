#! /bin/env python3
import math
import operator
from collections import defaultdict
from functools import reduce
from itertools import chain, groupby

from Bio.Data.SCOPData import protein_letters_3to1
from Bio.PDB import PDBParser, NeighborSearch
from Bio.PDB.ResidueDepth import ResidueDepth, min_dist, get_surface

from assignments.pdb.pdb import Residue
from ..pdb import pdb

import numpy as np

# tdo move to ipynb
np.set_printoptions(linewidth=200)

# is_buried
# is_exposed
# is_core (= is_buried)
# could be defined, by the ratio of amino of the chain (or model) vs other residues (but the first test data has no hetatms)
# no – try it by the count, (and compare it with ratio residue/water where there is water in the pdb)

# 3) visualize residue/water ratio bar chart with residue count bar chart (e.g. inner capsule with water should give differences,
# but is that rare?)
# (2) compute every CA to every other CA and just by that, that's ok, because what, 1000 calphas is just a million computations..
# 1) compare speed to NeighborSearch.search_all(radius, level='R'), visualize on bar chart, then define cutoff values for buried/exposed.


def residues_get_neighboring_water_count(structure, residues, radius_a):
    water_or_aa_atoms = []
    neighboring_water_count = dict()

    for r in structure.get_residues():
        if r.is_water:
            water_or_aa_atoms.extend(r.get_atoms())

        if r in residues:
            water_or_aa_atoms.extend(r.get_atoms())
            neighboring_water_count[r] = 0

    ns = NeighborSearch(water_or_aa_atoms)
    close_residue_pairs = ns.search_all(radius_a, level='R')  # undirected edges, unique

    for r1, r2 in close_residue_pairs:
        # move water to r1
        if r2.is_water:
            temp = r1
            r1 = r2
            r2 = temp

        if not r1.is_water:
            continue

        # now r1 water and r2 something

        # make r2 only aa residue
        if not r2.is_aa:
            continue

        neighboring_water_count[r2] += 1

    return neighboring_water_count


# todo neighbor_count = {r: 0 for r in residues} ?!!
# this does correlate positively with water count, because in core, probably the calphas are more far apart, as the need to fit
# side-chains too in there! Solution -> Count atoms in the vicinity of residue, not residues
def residues_get_neighbor_count(structure, residues, radius_a):
    ns = NeighborSearch(list(structure.get_atoms()))

    close_residue_pairs = ns.search_all(radius_a, level='R')  # undirected edges, unique

    neighbor_count = defaultdict(int)

    for r1r2 in close_residue_pairs:
        for r in r1r2:
            # ignore hetero atoms and water ('W') (possible residues arg doesn't contain them)
            if r not in residues:
                continue

            neighbor_count[r] += 1

    return neighbor_count


# for every atom in a residue,
#   choose the atom that has minimum number of neighboring entities
#   return the value for that residue
def residues_get_min_neighbor_entities(residues, radius_a, neighbor_function):
    residues = sorted(residues)

    assert len(residues) > 0

    residue_atoms = []
    result = dict()

    for r in residues:
        residue_atoms.extend(r.get_atoms())
        result[r] = 0

    ns = NeighborSearch(residue_atoms)
    close_atom_pairs = ns.search_all(radius_a, level='A')

    # 0) remove intra-residue atom pairs

    close_atom_pairs = list(filter(lambda pair: pair[0].get_parent() != pair[1].get_parent(), close_atom_pairs))

    # 1) make close_atom_pairs symmetric
    for i in range(len(close_atom_pairs)):
        a1, a2 = close_atom_pairs[i]
        close_atom_pairs.append((a2, a1))

    # 2) sort by residue and atom
    close_atom_pairs.sort(key=lambda pair: (pair[0].get_parent(), pair[0]))

    # 3) group_by residue and group_by atom

    gb_residue = groupby(close_atom_pairs, lambda pair: pair[0].get_parent())

    for residue, atoms_close_to_residue_atoms in gb_residue:
        min_neighbors = math.inf

        for residue_atom, atoms_close_to_atom in groupby(atoms_close_to_residue_atoms, lambda pair: pair[0]):
            neighbor_count = sum(1 for _ in neighbor_function(atoms_close_to_atom))

            min_neighbors = min(min_neighbors, neighbor_count)

        # 4) save a minimum (for residue's atoms)
        result[residue] = min_neighbors  # never infinity

    return result


def residues_get_min_neighbor_atoms_count(residues, radius_a):
    return residues_get_min_neighbor_entities(residues, radius_a, lambda pairs: pairs)


def residues_get_min_neighbor_residues_count(residues, radius_a):
    return residues_get_min_neighbor_entities(residues, radius_a, lambda pairs: set(a2.get_parent() for a1, a2 in pairs))



# 1) make close_atom_pairs symmetric
# 2) group_by residue and group_by atom
#       remove same-residue pairs
# 3) count atom's subgroup (either count all, or unique residues.
# 4) save a possible minimum (for residue's atoms)
# 5) all other residues should get 0

# minimum of all residue's atom


# def residues_get_min_neighbor_atoms_count(residues, radius_a):
#     residues = sorted(residues)
#
#     assert len(residues) > 0
#
#     residue_atoms = []
#     result = dict()
#
#     for r in residues:
#         residue_atoms.extend(r.get_atoms())
#         result[r] = 0
#
#     ns = NeighborSearch(residue_atoms)
#
#     close_atom_pairs = ns.search_all(radius_a, level='A')
#
#     close_atom_pairs.sort(lambda pair: (pair[0].get_parent(), pair[0]))
#
#     i_residues = iter(residues)
#     cur_residue = 'placeholder'
#
#     last_atom = None
#     neighboring_atoms = 0
#     min_for_cur_residue = math.inf
#
#     for atom, neighbor_atom in close_atom_pairs:
#         if cur_residue != atom.get_parent():
#             result[cur_residue] = neighboring_atoms
#             cur_residue = next(i_residues)
#
#         while cur_residue != atom.get_parent():
#             result[cur_residue] = 0
#             cur_residue = next(i_residues)
#
#         if last_atom != atom:
#             if neighboring_atoms < min_for_cur_residue:
#                 min_for_cur_residue = neighboring_atoms
#             neighboring_atoms = 0
#
#         neighboring_atoms += 1
#
#     if neighboring_atoms < min_for_cur_residue:
#         min_for_cur_residue = neighboring_atoms
#
#     result[cur_residue] = min_for_cur_residue
#
#     for r in i_residues:
#         result[r] = 0
#
#     del result['placeholder']
#
#     return result


# for every residue
#   count unique number of atoms within a radius from the residue
#       (however ignoring atoms within the residue itself)
def residues_get_neighbor_atom_count(residues, radius_a):
    residue_atoms = []
    neighbor_atom_count = dict()
    already_counted_atoms = defaultdict(set)

    for r in residues:
        residue_atoms.extend(r.get_atoms())
        neighbor_atom_count[r] = 0

    ns = NeighborSearch(residue_atoms)
    close_atom_pairs = ns.search_all(radius_a, level='A')

    for a1, a2 in close_atom_pairs:
        r1 = a1.get_parent()
        r2 = a2.get_parent()

        if r1 == r2:
            continue

        for neigh_a, r in ((a2, r1), (a1, r2)):
            if neigh_a in already_counted_atoms[r]:
                continue

            neighbor_atom_count[r] += 1
            already_counted_atoms[r].add(neigh_a)

    return neighbor_atom_count


def histogram_counts(counts, title):
    plt.hist(counts, bins=range(max(counts) + 2))
    plt.title(title)
    plt.show()


# https://stackoverflow.com/a/46328797
# {k:v1} {k:v2} -> {k: (v1, v2)}
def zip_mappings(*mappings):
    keys_sets = map(set, mappings)

    first_keys = next(keys_sets)
    common_keys = first_keys.copy()
    all_keys = first_keys.copy()

    for keys_set in keys_sets:
        common_keys &= keys_set
        all_keys |= keys_set

    if common_keys != all_keys:
        raise AssertionError('keys do not correspond')

    for key in all_keys:
        yield key, tuple(map(operator.itemgetter(key), mappings))


def scatter(dict1, dict2, title=''):
    x, y = np.transpose([vals for k, vals in zip_mappings(dict1, dict2)])

    plt.scatter(x, y, s=4)
    plt.title(title)
    plt.show()




def is_suface(residues_min_dist_to_surface, residue):
    # hodnota nastavena podle toho, ze jsem se podival do pymolu, co tak zhruba je na povrchu a co ne
    return residues_min_dist_to_surface[residue] < 1.85

def is_buried(residues_min_dist_to_surface, residue):
    return not is_suface(residues_min_dist_to_surface, residue)


def get_ratio_surface_to_buried(model, is_surface_fn, is_buried_fn=None):
    if is_buried_fn is None:
        is_buried_fn = lambda aa: not is_surface_fn(aa)

    return (sum((1 for aa in model.aa_residues if is_surface_fn(aa)))
     / sum((1 for aa in model.aa_residues
            if is_buried_fn(aa)))
     )


# returns dict {aa_three_letter_code: counts}
def get_aa_counts(aas):
    aa_counts = defaultdict(int)

    for aa in aas:
        aa_counts[aa.resname] += 1

    return aa_counts


def count_polar_aas(aa_counts):
    return sum((count for aa_code, count in aa_counts.items() if hp_model_tab[protein_letters_3to1[aa_code]] == 'P'))


if __name__ == '__main__':
    import os

    # structure = PDBParser(QUIET=True).get_structure('1tup', os.path.dirname(__file__) + '/../pdb/test_data/1tup.pdb')

    # structure = PDBParser(QUIET=True).get_structure('A2a receptor', os.path.dirname(__file__) + '/test_data/A2a_receptor.pdb')
    structure = PDBParser(QUIET=True).get_structure('hemoglobin', os.path.dirname(__file__) + '/test_data/hemoglobin(1b0b).pdb')

    model = next(structure.get_models())


    # for i in range(2, 15):
    #     histogram_counts(residues_get_neighboring_water_count(structure, i).values())
    #
    # quit()


    import matplotlib
    import matplotlib.pyplot as plt

    residue_depth = ResidueDepth(model)

    avg_residue_depth = {r: r_d for r, (r_d, calpha_d) in residue_depth}
    calpha_residue_depth = {r: calpha_d for r, (r_d, calpha_d) in residue_depth}


    surface = get_surface(model)
    min_dist_to_surface = {r: min(min_dist(atom.get_coord(), surface) for atom in r) for r in model.aa_residues}

    plt.hist(min_dist_to_surface.values())
    plt.title('avg residue depth')
    plt.show()


    # ratio of buried/exposed residues
    is_surface_fn = lambda aa: is_suface(min_dist_to_surface, aa)

    ratio_surface_to_buried = get_ratio_surface_to_buried(model, is_surface_fn)
    print('surface/buried: ', ratio_surface_to_buried)

    # histogram of buried/exposed amino acid types
    aa_types_surface = get_aa_counts(filter(is_surface_fn, model.aa_residues))
    aa_types_buried = get_aa_counts(filter(lambda aa: not is_surface_fn(aa), model.aa_residues))

    import seaborn as sns

    sns.barplot(list(aa_types_surface.keys()), list(aa_types_surface.values())).set_title('surface')
    plt.show()

    sns.barplot(list(aa_types_buried.keys()), list(aa_types_buried.values())).set_title('buried')
    plt.show()

    # percentage of polar in buried/exposed
    from Bio.Alphabet.Reduced import hp_model_tab

    for aa_counts, polar_non_polar_str in ((aa_types_buried, 'buried'), (aa_types_surface, 'surface')):
        polar_count = count_polar_aas(aa_counts)

        print(f'Ratio of polar aas {polar_non_polar_str}: ', polar_count/sum((count for count in aa_counts.values())))

    # quit()
    # compare two structures
    def closure():
        a2a = PDBParser(QUIET=True).get_structure('A2a receptor', os.path.dirname(__file__) + '/test_data/A2a_receptor.pdb')
        hemo = PDBParser(QUIET=True).get_structure('hemoglobin', os.path.dirname(__file__) + '/test_data/hemoglobin(1b0b).pdb')

        results = defaultdict(dict)

        for struct in (a2a, hemo):
            model = struct[0]

            surface = get_surface(model)
            min_dist_to_surface = {r: min(min_dist(atom.get_coord(), surface) for atom in r) for r in model.aa_residues}
            is_surface_fn = lambda aa: is_suface(min_dist_to_surface, aa)

            results[struct]['ratio surface to buried'] = get_ratio_surface_to_buried(model, lambda aa: is_suface(min_dist_to_surface, aa))

            aa_types_surface = get_aa_counts(filter(is_surface_fn, model.aa_residues))
            aa_types_buried = get_aa_counts(filter(lambda aa: not is_surface_fn(aa), model.aa_residues))

            for aa_counts, polar_non_polar_str in ((aa_types_buried, 'buried'), (aa_types_surface, 'surface')):
                polar_count = count_polar_aas(aa_counts)

                results[struct][f'Ratio of polar aas {polar_non_polar_str}'] =  polar_count / sum((count for count in aa_counts.values()))

        print(results)

    closure()

    # hemoglobin I z mušle LUCINA PECTINATA je monomerní rozpustný protein. V organismu je rozpuštěn v polární tekutině (hemolymfa),
    # má tedy na povrchu v naprosté většíně polární aminokyseliny
    # A2a je receptor spřažený s G proteinem. Obsahuje 7 transmembránových alfa-helixů. Jeho povrch je tedy ve značné míře v
    # kontaktu s alifatickými řetězci fosfolipidů lipidcké dvojvrstvy. V těchto hydrofobních místech jsou na povrchu proteinu převážně
    # kompatibilní hydrofobní aminokyseliny. Celkově je na povrchu proteinu zhruba polovina hydrofobních a polovina polárních aminokyselin.

    # hemoglobin má nižší poměr povrchových aminokyselin, protože jde o globulární protein (kulovitý tvar, nízký poměr
    # povrch:objem). A2a má protáhlejší, sudovitý tvar (7 agregovaných TM helixů) a má navíc kratší a delší netransmembránový helix,
    # což opět zvyšuje poměr povrchových aminokyselin.




    # quit()

    # residues_min_neighbor_atoms_count = residues_get_min_neighbor_atoms_count(model.aa_residues, 9)

    for radius_a in np.arange(3,10.5, 0.5):
        residues_min_neighbor_atoms_count = residues_get_min_neighbor_atoms_count(model.aa_residues, radius_a)
        # scatter(min_dist_to_surface, residues_min_neighbor_atoms_count, f'min_neighbor_atoms {radius_a}/depth')

    for radius_a in np.arange(3,10.5, 0.5):
        residues_neighbor_atoms_count = residues_get_neighbor_atom_count(model.aa_residues, radius_a)
        # scatter(min_dist_to_surface, residues_neighbor_atoms_count, f'neighbor_atoms {radius_a}/depth')


    neighboring_water_counts = residues_get_neighboring_water_count(model, model.aa_residues, 3)
    histogram_counts(neighboring_water_counts.values(), 'neighboring water')

    print(max(neighboring_water_counts.items(), key=lambda r_and_count: r_and_count[1]))
    res, count = max(neighboring_water_counts.items(), key=lambda t: t[1])

    neighbor_counts = residues_get_neighbor_count(model, model.aa_residues, 5)
    histogram_counts(neighbor_counts.values(), 'aa neighbors')

    neighbor_atom_counts = residues_get_neighbor_atom_count(model.aa_residues, 4.5)
    histogram_counts(neighbor_atom_counts.values(), 'aa atom neighbors')

    ratio_water_aas = {res: wc/ac for res, (wc, ac) in zip_mappings(neighboring_water_counts, neighbor_counts)}
    ratio_water_aa_atoms = {res: wc/aac for res, (wc, aac) in zip_mappings(neighboring_water_counts, neighbor_atom_counts)}

    normalized_neighbor_atom_counts = {r: ac/r.count_atoms for r, ac in neighbor_atom_counts.items()}  # tohle je stejne naprd,
    # zalezi na polomeru ze kt. beru neighbor atoms...
    ratio_water_normalized_aa_atoms = {k: wc/aac for k, (wc, aac) in zip_mappings(neighboring_water_counts, normalized_neighbor_atom_counts)}



    residues_min_neighbor_residues_count = residues_get_min_neighbor_residues_count(model.aa_residues, 8.5) #7
    residues_min_neighbor_atoms_count = residues_get_min_neighbor_atoms_count(model.aa_residues, 8.5)

    vars = np.transpose([values for k, values in zip_mappings(
        # based on surface computing program
        avg_residue_depth,
        calpha_residue_depth,
        min_dist_to_surface,

        # my variables
        neighboring_water_counts,
        neighbor_counts,
        neighbor_atom_counts,
        ratio_water_aas,
        ratio_water_aa_atoms,
        normalized_neighbor_atom_counts,
        ratio_water_normalized_aa_atoms,

        residues_min_neighbor_residues_count,
        residues_min_neighbor_atoms_count
    )])

    print(np.corrcoef(vars))  # prezentovat, ze to celkem funguje (korelace 0.75, i kdyz ta neni az tak dobra metrika, mozna clusterovani
    # nejak? Mohl bych pouzit 2-means na tom scatteru, ale s cim to budu porovnavat a jak skorovat?)
    # mozna mi 2-means urci cutoff hranice automaticky (s tim, ze tam maji byt dva clustery) a pak spocitam, false positives,
    # false negatives

    plt.hist(avg_residue_depth.values())
    plt.title('avg residue depth')
    plt.show()

    scatter(avg_residue_depth, neighboring_water_counts, 'water_count/depth')
    scatter(avg_residue_depth, neighbor_counts, 'aa_neighbors/depth')
    scatter(avg_residue_depth, normalized_neighbor_atom_counts, 'normalized_atom_neighbors/depth')
    scatter(avg_residue_depth, neighbor_atom_counts, 'neighbor_atoms/depth')
    scatter(avg_residue_depth, residues_min_neighbor_residues_count, 'min_neighbor_res/depth')
    scatter(avg_residue_depth, residues_min_neighbor_atoms_count, 'min_neighbor_atoms/depth')




    # pro zkoumani divnych ak (vysoka hloubka podle programu surface, ale já ukazuju málo neighbors a obráceně),
    # pouzit seaborn a nejak interaktivni chart, jestli jde, abych se na ty ak rychle podival v pymolu?


    # outliers = [
    #     (k, (x, y)) for k, (x, y)
    #         in sorted(zip_mappings(avg_residue_depth, neighbor_counts),
    #                   key=lambda k_xy: k_xy[1][0])
    #
    #     if y <= (9/6 * x + 5)
    # ]
    #
    # print([(k.get_parent(), k) for k, (x, y) in outliers])
    # print('f')


# todo map atom to residue
# against neighbor_atom_counts:
#   Maybe surface aa are larger and have more atoms -> can result in more neighbor atom counts
#   -> better to use ratio_water_aa_atoms

#y = 9/6 x + 5

# try min from res. atoms atom count (and res. count)
# jeste se podivat na ty sousedni zbytky, kolik maji sousedu? a to taky vzit v potaz
# limit >=85, někdy >=84 (leu B137)

# plot is on surface vs  min atoms

# set limits for buried/ surface (let's say surface up to 4?
# ratio sur/bur
# hist aa composition sur/bur
# polar in core vs sur
# two test data struct sur/bur and portion of pol (also on sur/bur and as a whole) Interpret differences