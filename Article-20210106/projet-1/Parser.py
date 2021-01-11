class DSSPParser:
    file = "/home/antho/Desktop/Bio info/projet 3/project-3-secondary-structure-prediction-anthzhou/ressources/dataset/CATH_info.txt"
    def __init__(self, directory, trainset_size = 3000):
        """
        Read the DSSP files and retrieve the amino acid sequences of the chains and their structures
        Create a list of the sequences and structures for the training set and a list for the testing set
        :param directory: str, destination folder where we can find the dssp files
        :param trainset_size: number of proteins in the training set vs the proteins in the testing set
        """
        self.directory = directory
        self.trainset_size = trainset_size
        self.trainset = []
        self.testset = []
        self.parse()

    def parse(self):
        """
        Extract the information for each identifiers in the CATH info file. Each identifier contains the id of the
        protein, as well as the chain of the desired sequence. The amino acid X, Z, B will be ignored, and all
        lower caps amino acids are converted to Cystein (C). The structures are converted into C (C, T, S and " "),
        H (H, G, I) or E (E and B).
        :return: None, but create a list of tuples of size 2 (str, str), respectively the amino acid sequence and the
        sequence of the structure
        """
        count = [] #count for trainset_size
        with open(self.file) as f:
            for line in f:
                data = line.split(" ")[0]
                filename = data[:-1]
                id = data[-1:]
                if (filename not in count):
                    count.append(filename)

                acid = ""
                structure = ""
                with open(self.directory+"/"+filename+".dssp") as dssp:
                    for i in range(28): #skip lines we don't need
                        next(dssp)
                    for line in dssp:
                        if (line[9] != " " and line[10] == " " and line[11] == id and line[13] not in ("*","!","B","Z","X")):
                            #amino acid sequence
                            if (line[13].islower()):
                                acid += "C"
                            else:
                                acid += line[13]

                            #sequence of the structure
                            if (line[16] in ("H","G","I")):
                                structure += "H"
                            elif (line[16] in ("E","B")):
                                structure += "E"
                            else:
                                structure += "C"

                if (len(count) > self.trainset_size):
                    self.testset.append((acid,structure))
                else:
                    self.trainset.append((acid,structure))


    def get_trainset(self):
        #print(len(self.trainset))
        #print(len(self.testset))
        return self.trainset

    def get_testset(self):
        return self.testset
