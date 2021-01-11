from math import log,sqrt,ceil
import numpy as np

class GOR:
    def __init__(self, Parser):
        """
        GOR III implementation, divided in a training part and a prediction. The training creates the
        counters/frequency used afterwards to predict which structure an amino acid is part of according to the
        propensity of this amino acid R to be part of a structure S (I(S,R)) and its neighbourhood if it is part
        of this structure (I(S,R,Rj) with -8<=j<=8).
        :param Parser: DSSPParser object
        """
        self.trainset = Parser.get_trainset()
        self.fs = {}            #{structure : number of appearance}
        self.fsr = {}           #{structure : {acid : number of appearance}}
        self.fsrneighbor = {}   #{structure : {acid : {acid_j+m : {position : number of appearance}}}}

        structure = ["H","E","C"] #list of all structure
        acid = ["A","R","N","D","C","E","Q","G","H","I","L","K","M","F","P","S","T","W","Y","V"] #list of all amino acid
        position = [-8,-7,-6,-5,-4,-3,-2,-1,1,2,3,4,5,6,7,8] #list of all neighbours position

        #initialise all counters to 0
        for i in structure:
            self.fs[i] = 0
            self.fsr[i] = {}
            self.fsrneighbor[i] = {}
            for j in acid:
                self.fsr[i][j] = 0
                self.fsrneighbor[i][j] = {}
                for k in acid:
                    self.fsrneighbor[i][j][k] = {}
                    for l in position:
                        self.fsrneighbor[i][j][k][l] = 0

        for i in range(len(self.trainset)):
            self.update_counters(self.trainset[i][0], self.trainset[i][1])


    def update_counters(self, sequence, structure):
        """
        Update the counters as asked : structure in which is the central amino acid, as well as its neighbourhood.
        :param sequence: str, amino acid sequence
        :param structure: str, sequence of the related structure
        :return: None, update the counters
        """
        if (len(sequence) != 0): #avoid empty string
            for index in range(len(sequence)):
                self.fs[structure[index]] += 1

                self.fsr[structure[index]][sequence[index]] += 1

                for m in [-8,-7,-6,-5,-4,-3,-2,-1,1,2,3,4,5,6,7,8]:
                    if (index+m < 0 or index+m >= len(sequence)):
                        continue
                    else:
                        self.fsrneighbor[structure[index]][sequence[index]][sequence[index+m]][m] += 1

    def predict(self, sequence):
        """
        Compute the probability of the three structures for each amino acid according to the counters and
        select the highest as the prediction (I(S,R) + somme(I(S,R,Rj) with -8<=j<=8 and j=!0))
        :param sequence: str, sequence to predict
        :return: str, sequence of the predicted structure
        """
        predicted_structure = ""
        probability = {} #probability of becoming the corresponding structure
        for i in range(len(sequence)):
            for j in self.fs:
                structure = j

                fn = sum(self.fs.values())
                fn_s = fn - self.fs[j]

                fnr = sum(self.fsr[k][sequence[i]] for k in self.fsr)
                fn_sr = fnr - self.fsr[j][sequence[i]]

                Isr = log(self.fsr[j][sequence[i]]/fn_sr) + log(fn_s/self.fs[j])

                Isrr_j = 0
                for m in [-8,-7,-6,-5,-4,-3,-2,-1,1,2,3,4,5,6,7,8]:
                    if (i+m < 0 or i+m >= len(sequence)):
                        continue
                    else:
                        fnrr_j = 0
                        for k in self.fsrneighbor:
                            fnrr_j += self.fsrneighbor[k][sequence[i]][sequence[i+m]][m]

                        fn_srr_j = fnrr_j - self.fsrneighbor[j][sequence[i]][sequence[i+m]][m]

                        Isrr_j += (log(self.fsrneighbor[j][sequence[i]][sequence[i+m]][m]/fn_srr_j) + log(fn_sr/self.fsr[j][sequence[i]]))

                I = Isr + Isrr_j
                probability[structure] = I
            predicted_structure += max(probability, key=probability.get)
        return predicted_structure

    def validate(self, test_set):
        """
        Specific function to test with a test set of variable length
        :param test_set: list of tuples of size 2 (str, str), respectively amino acid sequence and structure sequence
        :return: tuple of size 4 each of size 2((float, float), (float, float), (float, float), (float, float)),
        respectively mean Q3, MCC-H, MMC-E and MCC-C and their standard deviation, all rounded at 2 decimal places
        """
        Q3 = []
        H = []
        E = []
        C = []

        for i in range(len(test_set)):
            predicted_structure = self.predict(test_set[i][0])
            Q3.append(self.Q3(test_set[i][1], predicted_structure))
            H.append(self.MCC(test_set[i][1], predicted_structure, "H"))
            E.append(self.MCC(test_set[i][1], predicted_structure, "E"))
            C.append(self.MCC(test_set[i][1], predicted_structure, "C"))

        q3 = (round(np.mean(Q3),2), round(np.std(Q3),2))
        mcc_H = (round(np.mean(list(filter(None, H))),2), round(np.std(list(filter(None, H))),2))
        mcc_E = (round(np.mean(list(filter(None, E))),2), round(np.std(list(filter(None, E))),2))
        mcc_C = (round(np.mean(list(filter(None, C))),2), round(np.std(list(filter(None, C))),2))

        return (q3, mcc_H, mcc_E, mcc_C)

    def Q3(self, real, predicted):
        """
        Q3 computation
        :param real: str, real structure sequence
        :param predicted: str, predicted structure sequence
        :return: float rounded at 2 decimal, percentage of good predictions overall
        """
        residue_predicted = 0
        for i in range(len(real)):
            if (real[i] == predicted[i]):
                residue_predicted += 1

        return round((residue_predicted/len(real))*100,2)

    def MCC(self, real, predicted, structure):
        """
        Matthew coefficient correlation. Evaluates the performance of a classifier (here, we use the binary variation)
        Needs to have the True Positives, True Negatives, False Positives and False Negatives computed according to
        the structure given in parameter. Returns None in case of the denominator is non valid (division by 0)
        :param real: str, real structure sequence
        :param predicted: str, predicted structure sequence
        :param structure: str, structure you want to evaluate the MCC of
        :return: float rounded at 2 decimal, MCC value, or None if division error
        """
        TP = 0
        FP = 0
        FN = 0
        TN = 0
        for i in range(len(real)):
            if (real[i] == structure and predicted[i] == structure):
                TP += 1
            elif (real[i] != structure and predicted[i] == structure):
                FP += 1
            elif (real[i] == structure and predicted[i] != structure):
                FN += 1
            else:
                TN += 1
        denominator = sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN))
        if (denominator == 0):
            return None
        else:
            return round((((TP*TN)-(FP*FN))/denominator)*100,2)
