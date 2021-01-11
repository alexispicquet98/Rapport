from Matrix import *
from abc import ABC, abstractmethod

class Alignment(ABC):
    def __init__(self, seq1, seq2, I, E, submat, k):
        self.seq1 = seq1
        self.seq2 = seq2
        self.I = I
        self.E = E
        self.submat = submat
        self.k = k
        self.solution = []

        self.S = Matrix(len(seq1)+1, len(seq2)+1)
        self.V = Matrix(len(seq1)+1, len(seq2)+1)
        self.W = Matrix(len(seq1)+1, len(seq2)+1)

        self.S.set_value(0, 0, 0)
        for i in range(1, self.S.get_num_cols()):
            self.S.set_value(0, i, - I - ((i-1)*E))
        for j in range(1, self.S.get_num_rows()):
            self.S.set_value(j, 0, - I - ((j-1)*E))

        for i in range(1, self.V.get_num_rows()):
            self.V.set_value(i, 0, 0)
        for j in range(self.V.get_num_cols()):
            self.V.set_value(0, j, float("-inf"))

        for i in range(1, self.W.get_num_cols()):
            self.W.set_value(0, i, 0)
        for j in range(self.W.get_num_rows()):
            self.W.set_value(j, 0, float("-inf"))


    @abstractmethod
    def calculate_score(self, i, j):
        """
        Calculate the score of the cell in position i,j in the Score Matrix S, according to the values in S, V and W
        :param i: int, row index
        :param j: int, col index
        :return: float, the score of the cell
        """
        self.V.set_value(i, j, max((self.S[i-1, j]-self.I, self.V[i-1, j]-self.E)))
        self.W.set_value(i, j, max((self.S[i, j-1]-self.I, self.W[i, j-1]-self.E)))

        return max(self.S[i-1, j-1]+self.submat[self.seq1[i-1], self.seq2[j-1]], self.W[i, j], self.V[i, j])

    @abstractmethod
    def backtrack(self, i, j):
        """
        Use the Score Matrices to find the path taken and form the two strings of the alignment
        :param i: int, row index of current position
        :param j: int, col index of current position
        :return: list of tuples of size 5, (str, str, float, (int,int), (int,int)) for seq1 with gaps, seq2 with the
         inserted gaps and score of the alignment and the positions of beginning and end of the alignements sequences
        """
        pass

    def get_solution(self):
        """
        :return: list of tuples of size 5, (str, str, float, (int,int), (int,int)) for seq1 with gaps, seq2 with the
         inserted gaps and score of the alignment and the positions of beginning and end of the alignements sequences
        """
        return self.solution

    def compute_scores(self, alignment):
        """
        Compute the identity, similarity and number of gaps of an alignment
        :param alignment: list of tuples of size 5, (str, str, float, (int,int), (int,int)) for seq1 with gaps, seq2
        with the inserted gaps and score of the alignment and the positions of beginning and end of the alignments sequences
        :return: list of tuples of 3 floats (rounded two decimal places) respectively identity, similarity and gaps rates (in %)
        """
        solution = []
        for i in range(len(alignment)):
            identity = 0
            for j in range(len(alignment[i][0])):
                if (alignment[i][0][j] == alignment[i][1][j]):
                    identity += 1
            identity = (identity/len(alignment[i][0]))*100

            similarity = 0
            for j in range(len(alignment[i][0])):
                if (alignment[i][0][j] != "-" and alignment[i][1][j] != "-"):
                    if (self.submat[alignment[i][0][j], alignment[i][1][j]] > 0):
                        similarity += 1
            similarity = (similarity/len(alignment[i][0]))*100

            gaps_rates = 0
            for j in range(len(alignment[i][0])):
                if (alignment[i][0][j] == "-"):
                    gaps_rates += 1
                elif (alignment[i][1][j] == "-"):
                    gaps_rates += 1
            gaps_rates = (gaps_rates/(2*len(alignment[0][0])))*100

            solution.append((float("{:.2f}".format(identity)), float("{:.2f}".format(similarity)), float("{:.2f}".format(gaps_rates))))

        return solution

    @abstractmethod
    def run(self):
        """
        Run the alignment algorithm according to the parameters
        :return:
        """
        pass


class SmithWaterman(Alignment):
    def __init__(self, seq1, seq2, I, E, submat, k):
        """
        Local alignment algorithm
        :param seq1: str, first sequence of amino acids to align
        :param seq2: str, second sequence of amino acids to align
        :param I: float, initial gap value
        :param E: float, extension gap value
        :param submat: SubstitutionMatrix object
        :param k: int, maximum number of solutions to find
        """
        self.seq1 = seq1
        self.seq2 = seq2
        self.I = I
        self.E = E
        self.submat = submat
        self.k = k
        self.solution = []
        self.path = []

        self.S = Matrix(len(seq1)+1, self.submat.get_num_cols()+1)
        self.V = Matrix(len(seq1)+1, self.submat.get_num_cols()+1)
        self.W = Matrix(len(seq1)+1, self.submat.get_num_cols()+1)

        self.S.set_value(0, 0, 0)
        for i in range(1, self.S.get_num_cols()):
            self.S.set_value(0, i, 0)
        for j in range(1, self.S.get_num_rows()):
            self.S.set_value(j, 0, 0)

        for i in range(1, self.V.get_num_rows()):
            self.V.set_value(i, 0, 0)
        for j in range(self.V.get_num_cols()):
            self.V.set_value(0, j, float("-inf"))

        for i in range(1, self.W.get_num_cols()):
            self.W.set_value(0, i, 0)
        for j in range(self.W.get_num_rows()):
            self.W.set_value(j, 0, float("-inf"))

    def calculate_score(self, i, j):
        scoreV = max((self.S[i-1, j]-self.I, self.V[i-1, j]-self.E))
        if (scoreV < 0):
            self.V.set_value(i, j, 0)
        else:
            self.V.set_value(i, j, scoreV)

        scoreW = max((self.S[i, j-1]-self.I, self.W[i, j-1]-self.E))
        if (scoreW < 0):
            self.W.set_value(i, j, 0)
        else:
            self.W.set_value(i, j, scoreW)

        return max(self.S[i-1, j-1]+self.submat[self.seq1[i-1], j-1], self.W[i, j], self.V[i, j])

    def set_score(self):
        for i in range(1, self.S.get_num_rows()):
            for j in range(1, self.S.get_num_cols()):
                score = self.calculate_score(i,j)
                if (score < 0):
                    self.S.set_value(i, j, 0)
                else:
                    self.S.set_value(i, j, score)

    def recursive(self, i, j, seq1, solution):
        if (self.S[i, j] == 0):
            self.solution.append((seq1, round(self.S.get_max()[2], 2), (i+1, self.S.get_max()[0])))
            return
        self.path.append([i, j])
        if (i == 0 and j == 0):
            self.solution.append((self.seq1[i-1]+seq1, round(self.S.get_max()[2], 2), (i+1, self.S.get_max()[0])))
            return
        elif (i == 0):
            self.solution.append((self.seq1[i-1]+seq1, round(self.S.get_max()[2], 2), (i+1, self.S.get_max()[0])))
            return
        elif (j == 0):
            self.solution.append((seq1, round(self.S.get_max()[2], 2), (i+1, self.S.get_max()[0])))
            return
        elif (self.S[i, j] == self.S[i-1, j-1]+self.submat[self.seq1[i-1], j-1]):
            self.recursive(i-1, j-1, self.seq1[i-1]+seq1, solution)
        elif (self.S[i, j] == self.V[i, j]):
            self.recursive(i-1, j, self.seq1[i-1]+seq1, solution)
        elif (self.S[i, j] == self.W[i, j]):
            self.recursive(i, j-1, "-"+seq1, solution)

    def backtrack(self, i, j):
        solution = []
        self.recursive(i, j, "", solution)
        return solution

    def recalculate(self):
        """
        The path taken must be erased (values put to 0 and all the values in the matrix below the last cell of the path
        must be recomputed. The values at 0 must stay at 0.
        :return: None, but the ScoreMatrix S has been modified accordingly
        """
        for i in range(len(self.path)):
            self.S.set_value(self.path[i][0], self.path[i][1], 0)
            self.V.set_value(self.path[i][0], self.path[i][1], 0)
            self.W.set_value(self.path[i][0], self.path[i][1], 0)

        for i in range(self.path[-1][0], self.S.get_num_rows()):
            for j in range(self.path[-1][1], self.S.get_num_cols()):
                if ([i, j] in self.path):
                    continue
                score = self.calculate_score(i,j)
                if (score < 0):
                    self.S.set_value(i, j, 0)
                else:
                    self.S.set_value(i, j, score)

    def run(self):
        self.set_score()
        for i in range(self.k):
            self.backtrack(self.S.get_max()[0], self.S.get_max()[1])
            self.recalculate()


class SmithWatermanProfile(SmithWaterman):
    def __init__(self, seq, I, E, PSSM, k):
        """
        Local alignment algorithm
        :param seq: str, sequence of amino acids to find matches for the profile matrix
        :param I: float, initial gap value
        :param E: float, extension gap value
        :param PSSM: PSSM object
        :param k: int, maximum number of solutions to find
        """
        super().__init__(seq, '', I, E, PSSM, k) # to change as you prefer

    def get_solution(self):
        """
        :return: list of tuples of size 3, (str, float, (int,int)) for hit of seq1 with gaps, score of the
        alignment and a tuple for the positions of beginning and end of the alignments sequences
        """
        return self.solution
