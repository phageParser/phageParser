import pytest

from phage import Phage
from parsers.find_accession import PhageFinder


@pytest.fixture
def phage_finder():
    return PhageFinder('data/PhagesDB_Data.txt')


@pytest.mark.parametrize(
    'input, refseq, name',
    [
        ("ENA|AP013478|AP013478.1 Uncultured Mediterranean phage uvMED "
         "DNA, complete",
         "AP013478",
         "Uncultured Mediterranean phage uvMED DNA"),
        ("gi|526118413|ref|NC_021856.1| Bacillus phage phiNIT1 DNA, "
         "complete genome",
         "NC_021856.1",
         "Bacillus phage phiNIT1 DNA"),
        ("JoeDirt Final Sequence, 74914 bp including 10 bp 3' overhang "
         "(TCGAACGGCC), Cluster L1",
         "JF704108",
         "Mycobacteriophage JoeDirt"),
    ], ids=[
        'ena',
        'ncbi',
        'ad',
    ]
)
def test_phage_finder(phage_finder, phage_input, refseq, name):
    result = Phage(phage_input, phage_finder)
    assert refseq == result.refseq
    assert name == result.name
