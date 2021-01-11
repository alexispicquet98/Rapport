from math import log2, sqrt
import numpy as np

class Matrix:
    dict = {'A':0, 'R': 1, 'N':2, 'D':3, 'C':4, 'Q':5, 'E':6, 'G':7, 'H': 8, 'I':9,
            'L':10, 'K':11, 'M':12, 'F':13, 'P':14, 'S':15, 'T':16, 'W':17, 'Y':18, 'V':19}

    def __init__(self, nrows = 0, ncols = 0, value = None):
        self.rows = nrows
        self.cols = ncols
        self.mat = np.full((self.rows, self.cols), value, float)

    def __getitem__(self, item):
        #tranform the amino acid into index
        index = list(item)
        if index[0] in Matrix.dict:
            index[0] = Matrix.dict[index[0]]
        if index[1] in Matrix.dict:
            index[1] = Matrix.dict[index[1]]
        item = tuple(index)

        return self.mat[item]

    def get_num_cols(self):
        return self.cols

    def get_num_rows(self):
        return self.rows

    def get_max(self):
        """
        :return: tuple of size 3 of the row index, col index and value of the cell with the maximum value in the matrix
        If several occurences of the maximum value are found, the lowest indices are picked.
        """
        index = np.unravel_index(np.argmax(self.mat), self.mat.shape)
        return (index[0], index[1], np.amax(self.mat))

    def set_value(self, i, j, value):
        self.mat[i, j] = value

    def __str__(self):
        print(self.mat)

class SubstitutionMatrix(Matrix):
    """
    Format : free
    """
    def __init__(self, file):
        super().__init__(20, 20)
        self.parse_file(file)

    def parse_file(self, file):
        with open(file) as f:
            #skip comments
            for line in f:
                if line.startswith('#'):
                    continue
                else:
                    break
            #read lines but skip column 1
            for i in range(20):
                content = next(f).strip().split()[1:]
                for j in range(20):
                    self.mat[i, j] = int(content[j])


class ScoreMatrix(Matrix):
    """
    Format : list of list of float values
    """
    def __init__(self, seq1, seq2, type = "S"):
        super().__init__()


class PSSM(Matrix):
    PROB_AA = {'A': 0.0828, 'Q': 0.0394, 'L': 0.0967, 'S': 0.065,
               'R': 0.0553, 'E': 0.0676, 'K': 0.0585, 'T': 0.0532, 'N': 0.0405,
               'G': 0.0709, 'M': 0.0243, 'W': 0.0107, 'D': 0.0545, 'H': 0.0227,
               'F': 0.0386, 'Y': 0.0291, 'C': 0.0136, 'I': 0.0599, 'P': 0.0468,
               'V': 0.0687}

    def __init__(self, multiple_alignment):
        with open(multiple_alignment, 'r') as file:
            data = file.read().split('\n')
        self.N = len(data)/2 #number of sequences
        super().__init__(20, len(data[1]), 0)
        self.count_occurences(data)
        self.PWM()
        self.calculate_probabilities()

    def count_occurences(self, data):
        """
        Count the occurences of amino acids at each position
        """
        for i in range(1, len(data), 2):
            for j in range(len(data[i])):
                if (data[i][j] != "-" and data[i][j] != " "):
                    self.mat[self.dict[data[i][j]], j] += 1

    def PWM(self):
        """
        Position-weight matrix is a matrix using the occurences of amino acids at specific position and
        normalize the counts according to the number of sequences to obtain the frequency of the amino acid at
        a position. We use here also pseudocount (beta) to avoid a null frequency at a position. The gaps are ignored.
        Corresponds to the q(u,a) formula
        """
        alpha = self.N-1
        beta = sqrt(self.N)
        for i in range(self.get_num_rows()):
            a = [key for key, value in self.dict.items() if value == i][0] #find the right key
            for j in range(self.get_num_cols()):
                frequency = self.mat[i, j]/self.N
                self.mat[i,j] = ((alpha*frequency)+(beta*self.PROB_AA[a]))/(alpha+beta)

    def calculate_probabilities(self):
        """
        To do the transition from PWM to PSSM, you have to transform the normalized frequency to a probability, with
        the help of the log-odds. This is the logarithm between the probability of an event and
        of this event occurring randomly. The random event here is the probability of finding the amino acid anywhere,
        the probability of the event is our observation of having this amino acid at this specific position.
        Corresponds to m(u,a) formula
        """
        for i in range(self.get_num_rows()):
            a = [key for key, value in self.dict.items() if value == i][0] #find the right key
            for j in range(self.get_num_cols()):
                self.mat[i,j] = round(log2(self.mat[i, j]/self.PROB_AA[a]), 2)
