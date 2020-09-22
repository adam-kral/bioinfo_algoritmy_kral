from collections import namedtuple
from io import StringIO
import numpy as np

# this only file is not an assignment by itself, it just defines an abstract class AlignmentComputer, used by more assignments

class Alignment:
    """ Represents one result of alignment computation, there can be more results -- optimal alignments """
    def __init__(self, score, alignment_tuple):
        self.score = score
        self.alignment_tuple = alignment_tuple

    def __str__(self):
        return '\n'.join(self.alignment_tuple)


class AlignmentComputer:
    ALIGNMENT_CLS = Alignment

    def match(self, i, j, matrix):
        raise NotImplementedError

    def gap_s1(self, i, j, matrix):
        raise NotImplementedError

    def gap_s2(self, i, j, matrix):
        raise NotImplementedError

    # flags for alignment recovery (paths matrix)
    LEFT_FLAG = 1
    CORNER_FLAG = 1 << 1
    TOP_FLAG = 1 << 2

    def __init__(self, s1, s2):
        self.s1 = s1
        self.s2 = s2

    # I save the path (steps) when computing
    #  (there's no other option, if the score is so much variable, it is not possible to reconstruct the alignment
    #  in O(m+n) time if I did not save the paths -- example -> knowing if previous step was gap extension/open (different score values)
    #  requires knowing if the step preceding to the previous one was gap extension/open, and so forth till the beginning)

    def _compute_matrix(self):
        s1 = self.s1
        s2 = self.s2

        matrix = np.full((len(s1) + 1, len(s2) + 1), -np.inf)

        paths = np.full((len(s1) + 1, len(s2) + 1), 0, dtype=np.uint8)  # the first row and column is not theoretically needed (is trivial)
        # but having them leads to simpler code and less code. paths[0,0] has no flags (== 0) and I use it in my code

        matrix[0, 0] = 0

        for i in range(1, len(s1) + 1):
            matrix[i, 0] = matrix[i - 1, 0] + self.gap_s2(i - 1, 0, paths)
            paths[i, 0] = self.TOP_FLAG

        for j in range(1, len(s2) + 1):
            matrix[0, j] = matrix[0, j - 1] + self.gap_s1(0, j - 1, paths)
            paths[0, j] = self.LEFT_FLAG

        for i in range(1, len(s1) + 1):
            for j in range(1, len(s2) + 1):
                scores = (
                    matrix[i - 1, j - 1] + self.match(i - 1, j - 1, paths),
                    matrix[i - 1, j] + self.gap_s2(i - 1, j, paths),
                    matrix[i, j - 1] + self.gap_s1(i, j - 1, paths),
                )

                max_score = max(scores)

                for score, from_direction_flag in zip(scores, (self.CORNER_FLAG, self.TOP_FLAG, self.LEFT_FLAG)):
                    if score == max_score:
                        paths[i, j] |= from_direction_flag  # note that more flags can be present

                matrix[i, j] = max_score

        return matrix, paths

    def _compute_alignments_from_paths(self, score, paths):
        s1, s2 = self.s1, self.s2

        # alignment extension functions
        INSERT = lambda i, j: ('-', s2[j])  # insertion in the second sequence
        DELETE = lambda i, j: (s1[i], '-')  # deletion in the second sequence
        MATCH = lambda i, j: (s1[i], s2[j])

        alignment = (StringIO(), StringIO())

        # stack replacing recursion
        stack = []
        StepTo = namedtuple('StepTo', ['i', 'j', 'alignment_extend_operation', 'alignment_length'])

        # depth first search branching
        def dfs_get_branches_from(i, j, alignment_len):
            for predecessor_ij, direction_flag, predecessor_alignment_extend_operation in (
                    ((i, j - 1), self.LEFT_FLAG, INSERT),
                    ((i - 1, j - 1), self.CORNER_FLAG, MATCH),
                    ((i - 1, j), self.TOP_FLAG, DELETE)):

                if paths[i, j] & direction_flag:
                    yield StepTo(*predecessor_ij, predecessor_alignment_extend_operation, alignment_len)

        stack.extend(dfs_get_branches_from(len(s1), len(s2), 0))

        while stack:
            i, j, alignment_extend_operation, alignment_length = stack.pop()

            # extend alignment
            for str_io, next_letter in zip(alignment, alignment_extend_operation(i,j)):
                str_io.write(next_letter)

            if (i, j) == (0, 0):
                # yield current alignment
                al1_al2 = tuple(map(lambda str_io: ''.join(reversed(str_io.getvalue())), alignment))
                yield self.ALIGNMENT_CLS(score, al1_al2)

                # seek alignment StringIOs according to last branching of optimal alignments, if any
                if stack:
                    last_branch = stack[-1]

                    for str_io in alignment:
                        str_io.seek(last_branch.alignment_length)
                        str_io.truncate()

                continue

            # add all optimal branches from this state (alignment suffix) to the stack
            stack.extend(dfs_get_branches_from(i, j, alignment_length + 1))

    def compute(self):
        matrix, paths = self._compute_matrix()

        score = matrix[len(self.s1), len(self.s2)]

        yield from self._compute_alignments_from_paths(score, paths)
