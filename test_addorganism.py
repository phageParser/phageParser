from django.test import TestCase
from restapi.models import (
    Organism,
    Spacer,
    Repeat,
    LocusSpacerRepeat,
    AntiCRISPR,
    Locus
)
from add_organism import add_sequence_to_db

"""
Test Case for Adding Organism. Adds the test accession to the database
and compares the resulting database with the information from CrisprCASFinder
output (self.crisprs field). 

Run with:

$ python manage.py test

IMPORTANT: If getting the following error

>> Got an error creating the test database: permission denied to create database <<

Fix by:

1. Entering postgres shell:

$ psql phagedb

2. Authorize phageparser user to create databases:

==> ALTER USER phageparser CREATEDB;

"""

class AddOrganismTestCase(TestCase):

    test_accession = 'CP009273'
    entrez_email = 'wenxinchen11@hotmail.com'
    crisprs = [
    {
    "Name": "CP009273_1",
    "Start": 406532,
    "End": 406676,
    "DR_Consensus": "ATGCCTGATGCGACGCTTGCCGCGTCTTATCAGGCCTACAAAA",
    "Repeat_ID": "Unknown",
    "DR_Length": 43,
    "Spacers": 1,
    "Potential_Orientation": "ND",
    "CRISPRDirection": "ND",
    "Evidence_Level": 1,
    "Conservation_DRs": 97.6744186046512,
    "Conservation_Spacers": 0,
    "Regions": [
    {
    "Type": "LeftFLANK",
    "Start": 406432,
    "End": 406531,
    "Sequence": "AAACGCCGGTGCGTAAGGCGAAGCACGGTGATTCCAGCGGCGTACGCGGCGCTGCGTGGTTATGGCCACAAGAGTAAAAAACGTAGGCAATTGGCGCATC",
    "Leader": 0
    },
    {
    "Type": "DR",
    "Start": 406532,
    "End": 406574,
    "Sequence": "ATGCCTGATGCGACGCTTGCCGCGTCTTATCAGGCCTACAAAA"
    },
    {
    "Type": "Spacer",
    "Start": 406575,
    "End": 406633,
    "Sequence": "GGTGCCAGAACCGTAGGCCGGATAAGGCGTTCACGCCGCATCCGGCAATAAGTGCTCCG"
    },
    {
    "Type": "DR",
    "Start": 406634,
    "End": 406676,
    "Sequence": "ATGCCTGATGCGACGCTTGCCGCGTCTTATCAGGCCTGCAAAA"
    },
    {
    "Type": "RightFLANK",
    "Start": 406677,
    "End": 406776,
    "Sequence": "TGTGCCAGAACCGCGTAGGGCGGATAAGGCGTTCACGCCGCATCCGGCAATAAGTAATGAGCACCGAGACTATAACCTACCCCAGTGGTTTCGCCAGCAC",
    "Leader": 0
    }
    ]
    },
    {
    "Name": "CP009273_2",
    "Start": 546822,
    "End": 546918,
    "DR_Consensus": "GCTTGACGCGTCTTATCAGGCCTACAA",
    "Repeat_ID": "Unknown",
    "DR_Length": 27,
    "Spacers": 1,
    "Potential_Orientation": "ND",
    "CRISPRDirection": "ND",
    "Evidence_Level": 1,
    "Conservation_DRs": 96.2962962962963,
    "Conservation_Spacers": 0,
    "Regions": [
    {
    "Type": "LeftFLANK",
    "Start": 546722,
    "End": 546821,
    "Sequence": "GATTGGGGCGTTATCGCGAATTGAAGAGACGCTGGCGGGCGAAGCGGGGACCTGTATTTCGCTGTAGTCGTAGGCATTAGACATTTGTGCCTGATGCGAC",
    "Leader": 0
    },
    {
    "Type": "DR",
    "Start": 546822,
    "End": 546848,
    "Sequence": "GCTTGACGCGTCTTATCAGGCCTACAA"
    },
    {
    "Type": "Spacer",
    "Start": 546849,
    "End": 546891,
    "Sequence": "CCGGTGCCGCATCCGGCAATTGGTGCACAATGCCTGATGCGAT"
    },
    {
    "Type": "DR",
    "Start": 546892,
    "End": 546918,
    "Sequence": "GCTTGACGCATCTTATCAGGCCTACAA"
    },
    {
    "Type": "RightFLANK",
    "Start": 546919,
    "End": 547018,
    "Sequence": "TGGGTACCGGATCGGTAGGCCGGATAAGGCGTTTACGCCGCATCCGGCAAGAATAGAGCACCAGTTAACCGAACTTACTCTGCGCCCAAATCACGCCGCT",
    "Leader": 0
    }
    ]
    },
    {
    "Name": "CP009273_3",
    "Start": 2340680,
    "End": 2340806,
    "DR_Consensus": "TTTGTAGGCCTGATAAGACGCGCCAGCGTCGCATCAGGC",
    "Repeat_ID": "Unknown",
    "DR_Length": 39,
    "Spacers": 1,
    "Potential_Orientation": "-",
    "CRISPRDirection": "ND",
    "Evidence_Level": 1,
    "Conservation_DRs": 100,
    "Conservation_Spacers": 0,
    "Regions": [
    {
    "Type": "LeftFLANK",
    "Start": 2340580,
    "End": 2340679,
    "Sequence": "TGCCGTCAATCCAGGACGATGGCTGCGAAAGCGGCGCATGTAAGATCTGATATTGAGATGCCGGATGCGGCGTAAACGCCTTATCCGGCCTACGGCTCGG",
    "Leader": 0
    },
    {
    "Type": "DR",
    "Start": 2340680,
    "End": 2340718,
    "Sequence": "TTTGTAGGCCTGATAAGACGCGCCAGCGTCGCATCAGGC"
    },
    {
    "Type": "Spacer",
    "Start": 2340719,
    "End": 2340767,
    "Sequence": "TCCGGGTGCCGGATGCAGCGTGAACGCCTTATCCGGCCTACGGCTCGGA"
    },
    {
    "Type": "DR",
    "Start": 2340768,
    "End": 2340806,
    "Sequence": "TTTGTAGGCCTGATAAGACGCGCCAGCGTCGCATCAGGC"
    },
    {
    "Type": "RightFLANK",
    "Start": 2340807,
    "End": 2340906,
    "Sequence": "ACAGGATGCGGCGTAAAATGCCTTATCCGGCATTAAACTCCCAACAGGACACACTCATGGCATATACCACCTTTTCACAGACGAAAAATGATCAGCTCAA",
    "Leader": 1
    }
    ]
    },
    {
    "Name": "CP009273_4",
    "Start": 2871060,
    "End": 2871822,
    "DR_Consensus": "CGGTTTATCCCCGCTGGCGCGGGGAACTC",
    "Repeat_ID": "Unknown",
    "DR_Length": 29,
    "Spacers": 12,
    "Potential_Orientation": "-",
    "CRISPRDirection": "ND",
    "Evidence_Level": 4,
    "Conservation_DRs": 92.4677363820835,
    "Conservation_Spacers": 0,
    "Regions": [
    {
    "Type": "LeftFLANK",
    "Start": 2870960,
    "End": 2871059,
    "Sequence": "TTGGCGAAGGCGTCTTGATGGGTTTGAAAATGGGAGCTGGGAGTTCTACCGCAGAGGCGGGGGAACTCCAAGTGATATCCATCATCGCATCCAGTGCGCC",
    "Leader": 0
    },
    {
    "Type": "DR",
    "Start": 2871060,
    "End": 2871088,
    "Sequence": "CGGTTTATCCCCGCTGATGCGGGGAACAC"
    },
    {
    "Type": "Spacer",
    "Start": 2871089,
    "End": 2871120,
    "Sequence": "CAGCGTCAGGCGTGAAATCTCACCGTCGTTGC"
    },
    {
    "Type": "DR",
    "Start": 2871121,
    "End": 2871149,
    "Sequence": "CGGTTTATCCCTGCTGGCGCGGGGAACTC"
    },
    {
    "Type": "Spacer",
    "Start": 2871150,
    "End": 2871181,
    "Sequence": "TCGGTTCAGGCGTTGCAAACCTGGCTACCGGG"
    },
    {
    "Type": "DR",
    "Start": 2871182,
    "End": 2871210,
    "Sequence": "CGGTTTATCCCCGCTAACGCGGGGAACTC"
    },
    {
    "Type": "Spacer",
    "Start": 2871211,
    "End": 2871242,
    "Sequence": "GTAGTCCATCATTCCACCTATGTCTGAACTCC"
    },
    {
    "Type": "DR",
    "Start": 2871243,
    "End": 2871271,
    "Sequence": "CGGTTTATCCCCGCTGGCGCGGGGAACTC"
    },
    {
    "Type": "Spacer",
    "Start": 2871272,
    "End": 2871303,
    "Sequence": "CCGGGGGATAATGTTTACGGTCATGCGCCCCC"
    },
    {
    "Type": "DR",
    "Start": 2871304,
    "End": 2871332,
    "Sequence": "CGGTTTATCCCCGCTGGCGCGGGGAACTC"
    },
    {
    "Type": "Spacer",
    "Start": 2871333,
    "End": 2871364,
    "Sequence": "TGGGCGGCTTGCCTTGCAGCCAGCTCCAGCAG"
    },
    {
    "Type": "DR",
    "Start": 2871365,
    "End": 2871393,
    "Sequence": "CGGTTTATCCCCGCTGGCGCGGGGAACTC"
    },
    {
    "Type": "Spacer",
    "Start": 2871394,
    "End": 2871425,
    "Sequence": "AAGCTGGCTGGCAATCTCTTTCGGGGTGAGTC"
    },
    {
    "Type": "DR",
    "Start": 2871426,
    "End": 2871454,
    "Sequence": "CGGTTTATCCCCGCTGGCGCGGGGAACTC"
    },
    {
    "Type": "Spacer",
    "Start": 2871455,
    "End": 2871486,
    "Sequence": "TAGTTTCCGTATCTCCGGATTTATAAAGCTGA"
    },
    {
    "Type": "DR",
    "Start": 2871487,
    "End": 2871515,
    "Sequence": "CGGTTTATCCCCGCTGGCGCGGGGAACTC"
    },
    {
    "Type": "Spacer",
    "Start": 2871516,
    "End": 2871548,
    "Sequence": "GCAGGCGGCGACGCGCAGGGTATGCGCGATTCG"
    },
    {
    "Type": "DR",
    "Start": 2871549,
    "End": 2871577,
    "Sequence": "CGGTTTATCCCCGCTGGCGCGGGGAACTC"
    },
    {
    "Type": "Spacer",
    "Start": 2871578,
    "End": 2871610,
    "Sequence": "GCGACCGCTCAGAAATTCCAGACCCGATCCAAA"
    },
    {
    "Type": "DR",
    "Start": 2871611,
    "End": 2871639,
    "Sequence": "CGGTTTATCCCCGCTGGCGCGGGGAACTC"
    },
    {
    "Type": "Spacer",
    "Start": 2871640,
    "End": 2871671,
    "Sequence": "TCAACATTATCAATTACAACCGACAGGGAGCC"
    },
    {
    "Type": "DR",
    "Start": 2871672,
    "End": 2871700,
    "Sequence": "CGGTTTATCCCCGCTGGCGCGGGGAACTC"
    },
    {
    "Type": "Spacer",
    "Start": 2871701,
    "End": 2871732,
    "Sequence": "AGCGTGTTCGGCATCACCTTTGGCTTCGGCTG"
    },
    {
    "Type": "DR",
    "Start": 2871733,
    "End": 2871761,
    "Sequence": "CGGTTTATCCCCGCTGGCGCGGGGAACTC"
    },
    {
    "Type": "Spacer",
    "Start": 2871762,
    "End": 2871793,
    "Sequence": "TGCGTGAGCGTATCGCCGCGCGTCTGCGAAAG"
    },
    {
    "Type": "DR",
    "Start": 2871794,
    "End": 2871822,
    "Sequence": "CGGTTTATCCCCGCTGGCGCGGGGAACTC"
    },
    {
    "Type": "RightFLANK",
    "Start": 2871823,
    "End": 2871922,
    "Sequence": "TCTAAAAGTATACATTTGTTCTTAAAGCATTTTTTCCCATAAAAACAACCCACCAACCTTAATGTAACATTTCCTTATTATTAAAGATCAGCTAATTCTT",
    "Leader": 1
    }
    ]
    },
    {
    "Name": "CP009273_5",
    "Start": 2897373,
    "End": 2897766,
    "DR_Consensus": "GGTTTATCCCCGCTGGCGCGGGGAACAC",
    "Repeat_ID": "Unknown",
    "DR_Length": 28,
    "Spacers": 6,
    "Potential_Orientation": "-",
    "CRISPRDirection": "ND",
    "Evidence_Level": 4,
    "Conservation_DRs": 97.8868829336345,
    "Conservation_Spacers": 0,
    "Regions": [
    {
    "Type": "LeftFLANK",
    "Start": 2897273,
    "End": 2897372,
    "Sequence": "ATGGAGAACTATTTTGAACATGAGGTGTTACGTGGATATGTTGCTTATTACAAGTACTGCTAATATAAAAACTTGAGAAAGAGATAACGGGTTATATGGT",
    "Leader": 0
    },
    {
    "Type": "DR",
    "Start": 2897373,
    "End": 2897400,
    "Sequence": "GGTTTATCCCCGCTGGCGCGGGGAACTC"
    },
    {
    "Type": "Spacer",
    "Start": 2897401,
    "End": 2897433,
    "Sequence": "GACAGAACGGCCTCAGTAGTCTCGTCAGGCTCC"
    },
    {
    "Type": "DR",
    "Start": 2897434,
    "End": 2897461,
    "Sequence": "GGTTTATCCCCGCTGGCGCGGGGAACAC"
    },
    {
    "Type": "Spacer",
    "Start": 2897462,
    "End": 2897494,
    "Sequence": "CTGTTTTCGCAAATCTATGGACTATTGCTATTC"
    },
    {
    "Type": "DR",
    "Start": 2897495,
    "End": 2897522,
    "Sequence": "GGTTTATCCCCGCTGGCGCGGGGAACAC"
    },
    {
    "Type": "Spacer",
    "Start": 2897523,
    "End": 2897555,
    "Sequence": "GGGCGCACGGAATACAAAGCCGTGTATCTGCTC"
    },
    {
    "Type": "DR",
    "Start": 2897556,
    "End": 2897583,
    "Sequence": "GGTTTATCCCCGCTGGCGCGGGGAACAC"
    },
    {
    "Type": "Spacer",
    "Start": 2897584,
    "End": 2897616,
    "Sequence": "TGGCTCTGCAACAGCAGCACCCATGACCACGTC"
    },
    {
    "Type": "DR",
    "Start": 2897617,
    "End": 2897644,
    "Sequence": "GGTTTATCCCCGCTGGCGCGGGGAACAC"
    },
    {
    "Type": "Spacer",
    "Start": 2897645,
    "End": 2897677,
    "Sequence": "GAAATGCTGGTGAGCGTTAATGCCGCAAACACA"
    },
    {
    "Type": "DR",
    "Start": 2897678,
    "End": 2897705,
    "Sequence": "GGTTTATCCCCGCTGGCGCGGGGAACAC"
    },
    {
    "Type": "Spacer",
    "Start": 2897706,
    "End": 2897738,
    "Sequence": "ATTACGCCTTTTTGCGATTGCCCGGTTTTTGCC"
    },
    {
    "Type": "DR",
    "Start": 2897739,
    "End": 2897766,
    "Sequence": "GGTTTATCCCCGCTGGCGCGGGGAACAC"
    },
    {
    "Type": "RightFLANK",
    "Start": 2897767,
    "End": 2897866,
    "Sequence": "TCTAAACATAACCTATTATTAATTAATGATTTTTTAAGCCAGTCACAATCTACCAACTTTATAGTATCACACAAACAACACATCCATTATGTTAAAGAGC",
    "Leader": 1
    }
    ]
    }]

    # def setUp(self):
    #     test_accession = 'CP009273'
    #     entrez_email = 'wenxinchen11@hotmail.com'

    def test_add_organism(self):
        org = Organism.objects.filter(name=self.test_accession)
        self.assertEqual(org.count(), 0) #make sure organism is not in test database

        # add sequence to db
        add_sequence_to_db(self.test_accession, self.entrez_email)

        org = Organism.objects.filter(accession=self.test_accession)
        self.assertEqual(org.count(), 1) # make sure organism is added

        org = org.first()
        org_id = org.id
        loc_id = None

        for locus in self.crisprs:
            orientation = locus['Potential_Orientation'] == '+'
            if locus['Potential_Orientation'] == 'ND':
                orientation = None
            loci = Locus.objects.filter(organism_id=org_id, genomic_start=locus['Start'], genomic_end=locus['End'], consensus=locus['DR_Consensus'], orientation=orientation)
            self.assertEqual(loci .count(), 1)
            loc_id = loci.first().id
        
        seq_list = []
        seq = LocusSpacerRepeat.objects.filter(locus_id=loc_id).order_by('order')
        
        for s in seq.values():
            rep = Repeat.objects.filter(id=s['repeat_id'])
            self.assertEqual(rep.count(), 1)
            rep = rep.first()
            seq_list.append(rep.sequence)

            if s['spacer_id']:
                spac = Spacer.objects.filter(id=s['spacer_id'])
                self.assertEqual(spac.count(), 1)
                spac = spac.first()
                seq_list.append(spac.sequence)

        for i, region in enumerate(self.crisprs[-1]["Regions"][1:-1]):
            self.assertEqual(region['Sequence'], seq_list[i])

        




        