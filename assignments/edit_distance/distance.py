#! /bin/env python3

import math
from operator import add
from io import StringIO


# represents the result of Edit Distance calculation
# is created from the computed matrix
# for now, returns only one optimal alignment, even if there are many
class EditDistance:
    def __init__(self, distance, alignment):
        self.distance = distance
        self.alignment = alignment

    # wtf this too long, I should've built the alignment gradually in EditDistanceComputer.compute... if I knew it before
    @classmethod
    def from_matrix(cls, matrix, s1, s2):
        alignment = (StringIO(), StringIO())

        INSERT = lambda a, b: ('-', b)
        DELETE = lambda a, b: (a, '-')
        MATCH = lambda a, b: (a, b)

        i, j  = len(s1), len(s2)

        while (i,j) != (0,j) and (i,j) != (i,0):
            minimum = math.inf

            for (candidate_i, candidate_j), op in (((i-1,j-1), MATCH), ((i-1,j), DELETE), ((i,j-1), INSERT)):
                candidate_minimum = matrix[candidate_i][candidate_j]

                if candidate_minimum < minimum:
                    minimum = candidate_minimum
                    minimum_ij = (candidate_i, candidate_j)
                    minimum_op = op

            for s, a in zip(alignment, minimum_op(s1[i-1], s2[j-1])):
                s.write(a)

            (i, j) = minimum_ij

        # edge cases (edge of the matrix)

        while i > 0:
            for s, a in zip(alignment, DELETE(s1[i - 1], s2[j - 1])):
                s.write(a)
            i -= 1

        while j > 0:
            for s, a in zip(alignment, INSERT(s1[i - 1], s2[j - 1])):
                s.write(a)
            j -= 1

        # reverse buffer, unfortunately there is no easy way to read IO backwards, so I'm reversing the whole string
        alignment = tuple(map(lambda str_io: ''.join(reversed(str_io.getvalue())), alignment))

        return cls(matrix[len(s1)][len(s2)], alignment)


class EditDistanceComputer:
    W_INSERT = 1
    W_DELETE = 1
    W_UPDATE = 1

    def compute(self, s1, s2):
        matrix = [[None for i in range(len(s2)+1)] for i in range(len(s1)+1)]

        matrix[0][0] = 0

        for i in range(1, len(s1) + 1):
            matrix[i][0] = i * self.W_DELETE

        for j in range(1, len(s2) + 1):
            matrix[0][j] = j * self.W_INSERT

        for i in range(1, len(s1) + 1):
            for j in range(1, len(s2) + 1):
                matrix[i][j] = min(
                    matrix[i-1][j-1] + (0 if s1[i-1] == s2[j-1] else self.W_UPDATE),
                    matrix[i-1][j] + self.W_DELETE,
                    matrix[i][j-1] + self.W_INSERT,
                )

        return EditDistance.from_matrix(matrix, s1, s2)


def pprint(alignment):
    print(alignment[0])
    print(alignment[1])


if __name__ == '__main__':
    d = EditDistanceComputer().compute("clock", "overclocking")
    print('distance:', d.distance)
    pprint(d.alignment)
    print()

    d = EditDistanceComputer().compute("clock", "lacks")
    print('distance:', d.distance)
    pprint(d.alignment)
    print()

    d = EditDistanceComputer().compute("pitbull", "mr.worldwide")
    print('distance:', d.distance)
    pprint(d.alignment)
    print()
