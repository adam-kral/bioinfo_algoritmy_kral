#! /bin/env python3
from collections import defaultdict
from itertools import chain

from Bio.PDB import PDBParser, NeighborSearch

from ..pdb import pdb

# is_buried
# is_exposed
# is_core (= is_buried)
# could be defined, by the ratio of amino of the chain (or model) vs other residues (but the first test data has no hetatms)
# no – try it by the count, (and compare it with ratio residue/water where there is water in the pdb)

# 3) visualize residue/water ratio bar chart with residue count bar chart (e.g. inner capsule with water should give differences,
# but is that rare?)
# (2) compute every CA to every other CA and just by that, that's ok, because what, 1000 calphas is just a million computations..
# 1) compare speed to NeighborSearch.search_all(radius, level='R'), visualize on bar chart, then define cutoff values for buried/exposed.



def residues_get_neighboring_water_count(structure, radius_a):
    water_or_aa_atoms = []
    neighboring_water_count = dict()

    for r in structure.get_residues():
        if r.is_water:
            water_or_aa_atoms.extend(r.get_atoms())

        if r.is_aa_residue:
            water_or_aa_atoms.extend(r.get_atoms())
            neighboring_water_count[r] = 0


    ns = NeighborSearch(water_or_aa_atoms)
    close_residue_pairs = ns.search_all(radius_a, level='R')  # undirected edges, unique

    for r1, r2 in close_residue_pairs:
        # water ('W') and aa residue
        if r2.is_water:
            temp = r1
            r1 = r2
            r2 = temp

        if not r1.is_water:
            continue

        # now r1 water and r2 something

        # make r2 only aa residue
        if not r2.is_aa_residue:
            continue

        neighboring_water_count[r2] += 1

    return neighboring_water_count




def residues_get_neighbor_count(structure, radius_a):
    ns = NeighborSearch(list(structure.get_atoms()))
    close_residue_pairs = ns.search_all(radius_a, level='R')  # undirected edges, unique

    neighbor_count = defaultdict(int)

    for r1r2 in close_residue_pairs:
        for r in r1r2:
            # ignore hetero atoms and water ('W')
            if r.id[0].strip():  # is set hetero-flag
                continue

            neighbor_count[r] += 1

    return neighbor_count


def histogram_counts(counts, title):
    plt.hist(counts, bins=range(max(counts) + 2))
    plt.title(title)
    plt.show()


if __name__ == '__main__':
    import os
    structure = PDBParser(QUIET=True).get_structure('1tup', os.path.dirname(__file__) + '/../pdb/test_data/1tup.pdb')

    # for i in range(2, 15):
    #     histogram_counts(residues_get_neighboring_water_count(structure, i).values())
    #
    # quit()



    import matplotlib
    import matplotlib.pyplot as plt
    import numpy as np


    neighboring_water_counts = residues_get_neighboring_water_count(structure, 3.5)
    histogram_counts(neighboring_water_counts.values(), 'neighboring water')



    neighbor_counts = residues_get_neighbor_count(structure, 3.5)
    histogram_counts(neighbor_counts.values(), 'aa neighbors')

    xy = [[],[]]
    for aa in neighboring_water_counts.keys():
        xy[0].append(neighboring_water_counts[aa])
        xy[1].append(neighbor_counts[aa])

    print(np.corrcoef(xy))

