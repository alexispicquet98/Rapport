import unittest
from parameterized import parameterized

from Alignment import *

# match, gap dans sequence horizontale, gap dans sequence verticale

DIR = "ressources" # must be changed when testing locally (annoying relative path in Github workflow)

class TestNW(unittest.TestCase):
    def setUp(self):
        self.submat = SubstitutionMatrix(f"{DIR}/blosum62.txt")
        self.seq1short = "GGVTTF"
        self.seq2short = "MGGETFA"
        self.seq1long = "VPLPAGWEMAKTSSGQRYFLNHIDQTTTWQDPRK"
        self.seq2long = "GPLPPGWEERTHTDGRIFYINHNIKRTQWEDPRL"

    def test_linear_short(self):
        align = NeedlemanWunsch(self.seq1short, self.seq2short, 4, 4, self.submat, 1)
        align.run()
        self.assertEqual(align.get_solution(), [('-GGVTTF-', 'MGG-ETFA', 10, (0, 0), (6, 7))])
        self.assertEqual(align.compute_scores(align.get_solution()), [(50.0, 50.0, 18.75)])


    def test_affine_short(self):
        align = NeedlemanWunsch(self.seq1short, self.seq2short, 6, 1, self.submat, 1)
        align.run()
        self.assertEqual(align.get_solution(), [('-GGVTTF', 'MGGETFA', 5, (0, 0), (6, 7))])
        self.assertEqual(align.compute_scores(align.get_solution()), [(42.86, 42.86, 7.14)])


    def test_affine_multiple(self):
        align = NeedlemanWunsch(self.seq1long, self.seq2long, 4, 1, self.submat, 5)
        align.run()
        self.assertEqual(set(align.get_solution()), {('VPLPAGWEMAKT-SSGQRYFLNH-IDQTTTWQDPRK',
                                                      'GPLPPGWE-ERTHTDGRIFYINHNI-KRTQWEDPRL', 90.0, (0, 0), (34, 34)),
                                                     ('VPLPAGWEMAKT-SSGQRYFLNH-IDQTTTWQDPRK',
                                                      'GPLPPGWE-ERTHTDGRIFYINHNIKR-TQWEDPRL', 90.0, (0, 0), (34, 34)),
                                                     ('VPLPAGWEMAKT-SSGQRYFLNH-IDQTTTWQDPRK',
                                                      'GPLPPGWE-ERTHTDGRIFYINHNIKRT-QWEDPRL', 90.0, (0, 0), (34, 34)),
                                                     ('VPLPAGWEMAKT-SSGQRYFLNH-IDQTTTWQDPRK',
                                                      'GPLPPGWE-ERTHTDGRIFYINHNIKRTQ-WEDPRL', 90.0, (0, 0), (34, 34))})


class TestSW(unittest.TestCase):
    def setUp(self):
        self.submat = SubstitutionMatrix(f"{DIR}/blosum62.txt")
        self.seq1short = "THISLINE"
        self.seq2short = "ISALIGNED"
        self.seq1long = "MADEEKLPPGWEKRMSRSSGRVYYFNHITNASQWERPSGNSSSGGKNGQGEPARVRCSHLLVKHSQSRRPSSWRQEKITRTKEEALELINGYIQKIKSGEEDFESLASQFSDCSSAKARGDLGAFSRGQMQKPFEDASFALRTGEMSGPVFTDSGIHIILRTE"
        self.seq2long = "MDPGQQPPPQPAPQGQGQPPSQPPQGQGPPSGPGQPAPAATQAAPQAPPAGHQIVHVRGDSETDLEALFNAVMNPKTANVPQTVPMRLRKLPDSFFKPPEPKSHSRQASTDAGTAGALTPQHVRAHSSPASLQLGAVSPGTLTPTGVVSGPAATPTAQHLRQSSFEIPDDVPLPAGWEMAKTSSGQRYFLNHIDQTTTWQDPRKAMLSQMNVTAPTSPPVQQNMMNSASGPLPDGWEQAMTQDGEIYYINHKNKTTSWLDPRLDPRFAMNQRISQSAPVKQPPPLAPQSPQGGVMGGSNSNQQQQMRLQQLQMEKERLRLKQQELLRQAMRNINPSTANSPKCQELALRSQLPTLEQDGGTQNPVSSPGMSQELRTMTTNSSDPFLNSGTYHSRDESTDSGLSMSSYSVPRTPDDFLNSVDEMDTGDTINQSTLPSQQNRFPDYLEAIPGTNVDLGTLEGDGMNIEGEELMPSLQEALSSDILNDMESVLAATKLDKESFLTWL"

    def test_linear_short(self):
        align = SmithWaterman(self.seq1short, self.seq2short, 4, 4, self.submat, 2)
        align.run()
        self.assertEqual(set(align.get_solution()), {('IS-LI-NE', 'ISALIGNE', 19.0, (3, 1), (8, 8)),
                                                     ('IN', 'IS', 5.0, (6, 1), (7, 2))}, msg = "Dont forget to add the"
                                                                                               "start and stop positions")

    def test_affine_short(self):
        align = SmithWaterman(self.seq1short, self.seq2short, 10, 1, self.submat, 4)
        align.run()
        self.assertEqual(set(align.get_solution()), {('NE', 'NE', 11.0, (7, 7), (8, 8)),
                                                     ('ISLI', 'ISAL', 9.0, (3, 1), (6, 4)),
                                                     ('SLI', 'ALI', 9.0, (4, 3), (6, 5)),
                                                     ('IN', 'IS', 5.0, (6, 1), (7, 2))})

    def test_affine_long(self):
        align = SmithWaterman(self.seq1long, self.seq2long, 10, 1, self.submat, 3)
        align.run()
        self.assertEqual(set(align.get_solution()), {('DEEKLPPGWEKRMSRSSGRVYYFNHITNASQWERP',
                                                      'DDVPLPAGWEMAKT-SSGQRYFLNHIDQTTTWQDP', 82.0, (3, 169), (37, 202)),
                                                     ('LPPGWEKRMSRSSGRVYYFNHITNASQWERP',
                                                      'LPDGWEQAMTQ-DGEIYYINHKNKTTSWLDP', 77.0, (7, 232), (37, 261)),
                                                     ('SLASQFSDCSSAKARGDLGAFSRGQMQKPFEDASFALRTGEMSGPVFTDSGIHI',
                                                      'ALTPQHVRAHSSPASLQLGAVSPGTL-------T-P--TGVVSGPAATPTAQHL',
                                                      44.0, (105, 117), (158, 160))})

if __name__ == '__main__':
    unittest.main()