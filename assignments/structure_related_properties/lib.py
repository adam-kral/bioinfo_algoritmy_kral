import itertools
import math
import operator
from collections import defaultdict
from itertools import groupby

import matplotlib.pyplot as plt
import numpy as np
from Bio.Alphabet.Reduced import hp_model_tab
from Bio.Data.SCOPData import protein_letters_3to1
from Bio.PDB import PDBParser, NeighborSearch
from Bio.PDB.ResidueDepth import ResidueDepth, min_dist, get_surface


def residues_get_neighboring_water_count(structure, residues, radius_a):
    water_or_aa_atoms = []
    neighboring_water_count = dict()

    for r in residues:
        water_or_aa_atoms.extend(r.get_atoms())
        neighboring_water_count[r] = 0

    for r in structure.get_residues():
        if r.is_water:
            water_or_aa_atoms.extend(r.get_atoms())

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


def residues_get_neighbor_count(residues, radius_a):
    neighbor_count = {}
    atom_list = []

    for r in residues:
        neighbor_count[r] = 0
        atom_list.extend(r.get_atoms())

    ns = NeighborSearch(atom_list)

    close_residue_pairs = ns.search_all(radius_a, level='R')  # undirected edges, unique

    for r1r2 in close_residue_pairs:
        for r in r1r2:
            # ignore hetero atoms and water ('W') (possible residues arg doesn't contain them)
            if r not in neighbor_count:
                continue

            neighbor_count[r] += 1

    return neighbor_count


#
def residues_get_min_neighbor_entities(residues, radius_a, neighbor_function):
    """ for every atom in a residue,
           choose the atom that has minimum number of neighboring entities
               return the count for that residue
    """
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


def residues_get_neighbor_atom_count(residues, radius_a):
    """ for every residue
            count unique number of atoms within a radius from the residue however ignoring atoms within the residue itself """

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


def zip_mappings(*mappings):
    """ {k:v1} {k:v2} -> {k: (v1, v2)}

    mappings have to have the same key sets (nothing more, or less)
    based on https://stackoverflow.com/a/46328797
    """
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

    plt.scatter(x, y, s=2)
    plt.title(title)
    plt.show()


def is_surface(residues_min_dist_to_surface, residue):
    # hodnota nastavena podle toho, ze jsem se podival do pymolu, co tak zhruba je na povrchu a co ne
    return residues_min_dist_to_surface[residue] < 2.3

def is_buried(residues_min_dist_to_surface, residue):
    return not is_surface(residues_min_dist_to_surface, residue)


def get_ratio_surface_to_buried(model, is_surface_fn, is_buried_fn=None):
    if is_buried_fn is None:
        is_buried_fn = lambda aa: not is_surface_fn(aa)

    return (sum((1 for aa in model.aa_residues if is_surface_fn(aa)))
     / sum((1 for aa in model.aa_residues
            if is_buried_fn(aa)))
     )


def get_aa_counts(aas):
    """ Returns dict {aa_three_letter_code: counts} """
    aa_counts = defaultdict(int)

    for aa in aas:
        aa_counts[aa.resname] += 1

    return aa_counts


def count_polar_aas(aa_counts):
    return sum((count for aa_code, count in aa_counts.items() if hp_model_tab[protein_letters_3to1[aa_code]] == 'P'))


def get_model_without_waters(model):
    """ model without hoh e.g. for surface analysis (not interested in surface of solvation (water) layer!, only protein) """
    model_without_waters = model.copy()

    for chain in model_without_waters:
        residue_ids_to_delete = []
        for residue in chain:
            if residue.is_water:
                residue_ids_to_delete.append(residue.id)

        for r_id in residue_ids_to_delete:
            del chain[r_id]

    return model_without_waters


if __name__ == '__main__':
    import os


    # structure = PDBParser(QUIET=True).get_structure('A2a receptor', os.path.dirname(__file__) + '/test_data/A2a_receptor.pdb')
    from assignments.pdb.pdb import attach_custom_classes_to_pdb_structure_builder

    attach_custom_classes_to_pdb_structure_builder()
    structure = PDBParser(QUIET=True).get_structure('1tup', os.path.dirname(__file__) + '/test_data/1tup.pdb')

    model = next(structure.get_models())

    # model without hoh for surface analysis (not interested in surface of solvation (water) layer!, only protein)
    model_without_waters = get_model_without_waters(model)
    cp_to_original_res = lambda r: model[r.get_parent().id][r.id]  # gets the corresponding residue in the original model

    residue_depth = ResidueDepth(model_without_waters)

    avg_residue_depth = {cp_to_original_res(r): r_d for r, (r_d, calpha_d) in residue_depth}
    calpha_residue_depth = {cp_to_original_res(r): calpha_d for r, (r_d, calpha_d) in residue_depth}

    surface = get_surface(model_without_waters)
    min_dist_to_surface = {r: min(min_dist(atom.get_coord(), surface) for atom in r) for r in model.aa_residues}

    plt.hist(min_dist_to_surface.values())
    plt.title('avg residue depth')
    plt.show()

    # ratio of buried/exposed residues
    is_surface_fn = lambda aa: is_surface(min_dist_to_surface, aa)

    ratio_surface_to_buried = get_ratio_surface_to_buried(model, is_surface_fn)
    print('surface/buried: ', ratio_surface_to_buried)

    # histogram of buried/exposed amino acid types
    aa_types_surface = get_aa_counts(filter(is_surface_fn, model.aa_residues))
    aa_types_buried = get_aa_counts(filter(lambda aa: not is_surface_fn(aa), model.aa_residues))

    import seaborn as sns

    aa_types = set(itertools.chain(aa_types_surface, aa_types_buried))
    palette = dict(zip(aa_types, sns.color_palette('husl', len(aa_types))))

    sns.barplot(list(aa_types_surface.keys()), list(aa_types_surface.values()), palette=palette).set_title('surface')
    plt.show()

    sns.barplot(list(aa_types_buried.keys()), list(aa_types_buried.values()), palette=palette).set_title('buried')
    plt.show()

    # percentage of polar in buried/exposed

    for aa_counts, location_str in ((aa_types_buried, 'buried'), (aa_types_surface, 'surface')):
        polar_count = count_polar_aas(aa_counts)

        print(f'Ratio of polar aas {location_str}: ', polar_count/sum((count for count in aa_counts.values())))

    # compare two structures
    def closure():
        a2a = PDBParser(QUIET=True).get_structure('A2a receptor', os.path.dirname(__file__) + '/test_data/A2a_receptor.pdb')
        hemo = PDBParser(QUIET=True).get_structure('hemoglobin', os.path.dirname(__file__) + '/test_data/hemoglobin(1b0b).pdb')

        results = defaultdict(dict)

        for struct in (a2a, hemo):
            model = struct[0]

            surface = get_surface(get_model_without_waters(model))
            min_dist_to_surface = {r: min(min_dist(atom.get_coord(), surface) for atom in r) for r in model.aa_residues}
            is_surface_fn = lambda aa: is_surface(min_dist_to_surface, aa)

            results[struct]['ratio surface to buried'] = get_ratio_surface_to_buried(model, lambda aa: is_surface(min_dist_to_surface, aa))

            aa_types_surface = get_aa_counts(filter(is_surface_fn, model.aa_residues))
            aa_types_buried = get_aa_counts(filter(lambda aa: not is_surface_fn(aa), model.aa_residues))

            for aa_counts, location_str in ((aa_types_buried, 'buried'), (aa_types_surface, 'surface')):
                polar_count = count_polar_aas(aa_counts)

                results[struct][f'Ratio of polar aas {location_str}'] =  polar_count / sum((count for count in aa_counts.values()))

        print(results)

    closure()

    residues_min_neighbor_atoms_count = residues_get_min_neighbor_atoms_count(model.aa_residues, 9)

    for radius_a in np.arange(3, 10.5, 0.5):
        residues_min_neighbor_atoms_count = residues_get_min_neighbor_atoms_count(model.aa_residues, radius_a)
        scatter(min_dist_to_surface, residues_min_neighbor_atoms_count, f'min_neighbor_atoms {radius_a}/depth')

    for radius_a in np.arange(3, 10.5, 0.5):
        residues_neighbor_atoms_count = residues_get_neighbor_atom_count(model.aa_residues, radius_a)
        scatter(min_dist_to_surface, residues_neighbor_atoms_count, f'neighbor_atoms {radius_a}/depth')

    neighboring_water_counts = residues_get_neighboring_water_count(model, model.aa_residues, 3.5)
    histogram_counts(neighboring_water_counts.values(), 'neighboring water')

    neighbor_counts = residues_get_neighbor_count(model.aa_residues, 5)
    histogram_counts(neighbor_counts.values(), 'aa neighbors')

    neighbor_atom_counts = residues_get_neighbor_atom_count(model.aa_residues, 4.5)
    histogram_counts(neighbor_atom_counts.values(), 'aa atom neighbors')

    ratio_water_aas = {res: wc/ac for res, (wc, ac) in zip_mappings(neighboring_water_counts, neighbor_counts)}
    ratio_water_aa_atoms = {res: wc/aac for res, (wc, aac) in zip_mappings(neighboring_water_counts, neighbor_atom_counts)}

    normalized_neighbor_atom_counts = {r: ac/r.count_atoms for r, ac in neighbor_atom_counts.items()}  # tohle je stejne naprd,
    # neni to linearni s poctem atomu (spis afinní funkce) a ještě záleží na tvaru ak
    ratio_water_normalized_aa_atoms = {k: wc/aac for k, (wc, aac) in zip_mappings(neighboring_water_counts, normalized_neighbor_atom_counts)}



    residues_min_neighbor_residues_count = residues_get_min_neighbor_residues_count(model.aa_residues, 8.5)
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

    np.set_printoptions(linewidth=200)  # for the large correlation matrix
    print(np.corrcoef(vars))

    plt.hist(avg_residue_depth.values())
    plt.title('avg residue depth')
    plt.show()

    scatter(min_dist_to_surface, neighboring_water_counts, 'water_count/depth')
    scatter(min_dist_to_surface, neighbor_counts, 'aa_neighbors/depth')
    scatter(min_dist_to_surface, normalized_neighbor_atom_counts, 'normalized_atom_neighbors/depth')
    scatter(min_dist_to_surface, neighbor_atom_counts, 'neighbor_atoms/depth')
    scatter(min_dist_to_surface, residues_min_neighbor_residues_count, 'min_neighbor_res/depth')
    scatter(min_dist_to_surface, residues_min_neighbor_atoms_count, 'min_neighbor_atoms/depth')
