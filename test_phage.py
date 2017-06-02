import unittest
from phage import Phage
from parsers.find_accession import PhageFinder

class TestPhage(unittest.TestCase):

  phage_finder = PhageFinder('data/PhagesDB_Data.txt')

  def test_ena(self):
    ena1 = "ENA|AP013478|AP013478.1 Uncultured Mediterranean phage uvMED DNA, complete"
    ena1_phage = Phage(ena1, TestPhage.phage_finder)
    self.assertEqual("AP013478", ena1_phage.refseq)
    self.assertEqual("Uncultured Mediterranean phage uvMED DNA", ena1_phage.name)


  def test_ncbi(self):
    ncbi1 = "gi|526118413|ref|NC_021856.1| Bacillus phage phiNIT1 DNA, complete genome"
    ncbi1_phage = Phage(ncbi1, TestPhage.phage_finder)
    self.assertEqual("NC_021856.1", ncbi1_phage.refseq)
    self.assertEqual("Bacillus phage phiNIT1 DNA", ncbi1_phage.name)

  def test_ad(self):
    ad1 = "JoeDirt Final Sequence, 74914 bp including 10 bp 3' overhang (TCGAACGGCC), Cluster L1"
    ad1_phage = Phage(ad1, TestPhage.phage_finder)
    self.assertEqual("JF704108", ad1_phage.refseq)
    self.assertEqual("Mycobacteriophage JoeDirt", ad1_phage.name)


if __name__ == '__main__':
  unittest.main()
