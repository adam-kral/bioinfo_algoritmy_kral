#! /bin/env python3

from Bio import SeqIO


def parse_fasta(handle_or_filename):
    return SeqIO.parse(handle_or_filename, 'fasta')


if __name__ == '__main__':
    import os

    # Read in a FASTA files with an arbitrary number of molecules.
    for seq_record in parse_fasta(os.path.dirname(__file__) + os.path.sep + 'test_data/seq1.fasta'):
        # Obtain description/sequence of any molecules.
        print(seq_record.description)
        print(seq_record.seq)

        # Return sequence length for given sequence.
        print(len(seq_record.seq))

        # Return subsequence of given sequence.
        print(seq_record.seq[2:10])
