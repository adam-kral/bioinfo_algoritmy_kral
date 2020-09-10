#! /bin/env python3

from io import StringIO

from Bio.Align import MultipleSeqAlignment
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

from Bio.SubsMat import MatrixInfo

# parser based on http://meme-suite.org/doc/clustalw-format.html


# block = block of text where a line contains a chunk of the aligned sequence
class MSAFromBlocksBuilder:
    def __init__(self, first_block):
        self.seqs = []

        for seq_name, seq_symbols in first_block:
            self.seqs.append((seq_name, StringIO(seq_symbols)))

    def add_block(self, block):
        for i, (seq_name, seq_symbols) in enumerate(block):
            if self.seqs[i][0] != seq_name:
                raise InvalidFileFormatException('unexpected sequence ID for a line')

            self.seqs[i][1].write(seq_symbols)

    def buildMSA(self):
        return MSA(records=(SeqRecord(Seq(seq_io.getvalue()), name, name, name) for name, seq_io in self.seqs))


class InvalidFileFormatException(Exception):
    pass


def parse_block_line(line):
    parts = line.split()

    if len(parts) < 2:
        raise InvalidFileFormatException(f'invalid sequence line, name and symbols must be present: `{line}` ')

    seq_name = parts[0]
    seq_symbols = parts[1]
    # ignoring potential cumulative count of symbols (optional parts[2])
    # it's probably mainly for the user if he wishes to orient himself in the text file
    # BioPython::MultipleSeqAlignment checks that all seqs have the same length, so an accidental user induced file change is covered

    return seq_name, seq_symbols


def read_line_skip_empty_lines(file):
    line = file.readline()
    while line == '\n':
        line = file.readline()

    return line


def parse_block(file):
    block_lines = []

    line = read_line_skip_empty_lines(file)

    # there is no block
    if line == '':
        return

    while line.strip('*:. ') != '\n':   # stripping potential degree-of-conservation line
        block_line_tuple = parse_block_line(line)
        block_lines.append(block_line_tuple)

        line = file.readline()

    # todo resolve,or is resolved with class PECF?
    # # line after degree-of-conservation line must be empty
    # if line != '\n' and file.readline() != '\n':
    #     raise InvalidFileFormatException('line dividing blocks must be empty')
    # ^^ apparently, not at the end of the file :(

    return block_lines


class ProperlyEndedClustalFile:
    """ Wrapper ensuring, that even after the last seq block there is an empty line (thus adding 2 lines, sometimes a file does not even
    end with a newline!) """

    def __init__(self, file):
        self.file = file
        self.added_newlines = 0

    def readline(self):
        line = self.file.readline()

        if not line.endswith('\n') and self.added_newlines < 2:
            self.added_newlines += 1
            return line + '\n'

        return line


def parse_clustal(file_handle):
    file = ProperlyEndedClustalFile(file_handle)

    if not file.readline().startswith('CLUSTAL'):
        raise InvalidFileFormatException('first line should start with `CLUSTAL`')

    if file.readline() != '\n':
        raise InvalidFileFormatException('second line must be empty')

    block = parse_block(file)

    if not block:
        raise InvalidFileFormatException('no sequence blocks')

    builder = MSAFromBlocksBuilder(block)

    while block:
        builder.add_block(block)
        block = parse_block(file)

    return builder.buildMSA()


class MSA(MultipleSeqAlignment):
    def sum_of_pairs_column(self, column_index, scoring_matrix, gap_score=0):
        """
        :param dict scoring_matrix: key = pair of symbols, value = score for that pair
            example: {('A', 'G'): 5, ...}
            can be only triangular
        """
        sum = 0

        for i in range(len(self)):
            for j in range(i+1, len(self)):
                pair = self[i, column_index], self[j, column_index]

                if '-' in pair:
                    sum += gap_score
                else:
                    if pair not in scoring_matrix:
                        pair = pair[1], pair[0]  # switch pair elements, in case matrix is just triangular

                    sum += scoring_matrix[pair]

        return sum

    def sum_of_pairs(self, scoring_matrix, gap_score=0, start_pos=0, end_pos=None):
        if end_pos is None:
            end_pos = len(self[0])

        return sum((self.sum_of_pairs_column(i, scoring_matrix, gap_score) for i in range(start_pos, end_pos)))

    def top_n_scoring_positions(self, n, scoring_matrix, gap_score=0):
        """ returns list of tuples (position, score) """
        column_scores = []

        for i in range(self.get_alignment_length()):
            column_scores.append((i, self.sum_of_pairs_column(i, scoring_matrix, gap_score)))

        column_scores.sort(key=lambda t: t[1], reverse=True)

        return column_scores[:n]


if __name__ == '__main__':
    import os

    # Read and parse MSA.
    with open(os.path.dirname(__file__) + os.path.sep + 'test_data/clustal_assigned') as f:
        msa = parse_clustal(f)

    print(msa)

    # Retrieve sequence by its position (but not by ID which is not unique in the assigned input data)
    print(msa[0])

    # Retrieve given column from the MSA
    print(msa[:,60])

    # Retrieve sum of pairs score of a column and whole MSA with respect to given scoring matrix.
    scoring_matrix = MatrixInfo.blosum62

    print(msa.sum_of_pairs_column(60, scoring_matrix))
    print(msa.sum_of_pairs(scoring_matrix))

    # You should be able to specify a sequence in the MSA and the method would return conservation scores for the positions of the MSA.
    print('Sum of pairs:', msa.sum_of_pairs(scoring_matrix, start_pos=50, end_pos=100))

    # Moreover, one should be able to identify top N scoring positions in the MSA.
    print(msa.top_n_scoring_positions(10, scoring_matrix))

    for seq_record in msa:
        print(seq_record.id)

    # Read and parse MSA.
    with open(os.path.dirname(__file__) + os.path.sep + 'test_data/meme-suit_clustalw-format') as f:
        msa = parse_clustal(f)

    print(msa)