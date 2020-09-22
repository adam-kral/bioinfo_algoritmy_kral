class IncompatibleLengthException(Exception):
    pass


def hamming_distance(seq1, seq2):
    """ :param seq1/2, indexable container, its elements comparable with != """
    if len(seq1) != len(seq2):
        raise IncompatibleLengthException()

    d = 0

    for i in range(len(seq1)):
        if seq1[i] != seq2[i]:
            d += 1

    return d


if __name__ == '__main__':
    print(hamming_distance('abcd', 'abcd'))
    print(hamming_distance('abcd', 'bbcw'))
    print(hamming_distance('clock', 'lacks'))

    import os
    from ..fasta.fasta import parse_fasta

    seq_record = next(parse_fasta(os.path.dirname(__file__) + os.path.sep + '../fasta/test_data/seq1.fasta'))
    print(hamming_distance(seq_record.seq, seq_record.seq))

