#! /bin/env python3


from assignments.alignment import AlignmentComputer, Alignment


# represents the result of Edit Distance calculation
# is created from the computed matrix
class EditDistanceAlignment(Alignment):
    # Alignment class works with scores, not distances. The higher score, the better match (the lower distance)
    @property
    def distance(self):
        return -self.score


class EditDistanceComputer(AlignmentComputer):
    ALIGNMENT_CLS = EditDistanceAlignment

    # AlignmentComputer works with scores, not distances. The higher score, the better match (the lower distance)
    SCORE_INSERT = -1
    SCORE_DELETE = -1
    SCORE_UPDATE = -1

    def match(self, i, j, paths):
        return 0 if self.s1[i] == self.s2[j] else self.SCORE_UPDATE

    def gap_s1(self, i, j, paths):
        return self.SCORE_DELETE

    def gap_s2(self, i, j, paths):
        return self.SCORE_INSERT


if __name__ == '__main__':
    alignments = EditDistanceComputer("clock", "overclocking").compute()
    for al in alignments:
        print('distance:', al.distance)
        print(al)
        print()

    alignments = EditDistanceComputer("clock", "lacks").compute()
    for al in alignments:
        print('distance:', al.distance)
        print(al)
        print()

    alignments = EditDistanceComputer("pitbull", "mr.worldwide").compute()
    for al in alignments:
        print('distance:', al.distance)
        print(al)
        print()

