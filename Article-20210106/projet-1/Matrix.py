import numpy as np

class Matrix:
    def __init__(self, nrows = 0, ncols = 0, value = None):
        self.rows = nrows
        self.cols = ncols
        self.mat = np.full((self.rows, self.cols), value)

    def __getitem__(self, item):
        dict = {'A':0, 'R': 1, 'N':2, 'D':3, 'C':4, 'Q':5, 'E':6, 'G':7, 'H': 8, 'I':9,
                'L':10, 'K':11, 'M':12, 'F':13, 'P':14, 'S':15, 'T':16, 'W':17, 'Y':18, 'V':19}
        #tranform the amino acid into index
        index = list(item)
        if index[0] in dict:
            index[0] = dict[index[0]]
        if index[1] in dict:
            index[1] = dict[index[1]]
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


## Proposition of implementation, but free for you to use it or not, not included in tests
class ScoreMatrix(Matrix):
    """
    Format : list of list of float values
    """
    def __init__(self, seq1, seq2, I, E, local=False, type = "S"):
        super().__init__()
