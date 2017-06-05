import pytest

# Create your tests here.
from restapi.models import Organism, Repeat


@pytest.fixture
def organism():
    organism_name = "Organism Name"
    accession_number = "NC_000000"
    return Organism(name=organism_name, accession=accession_number)


@pytest.mark.django_db
def test_model_can_create_an_organism(organism):
    """Test the organism model can create an Organism."""
    old_count = Organism.objects.count()
    organism.save()
    new_count = Organism.objects.count()
    assert old_count != new_count


@pytest.fixture
def repeat():
    sequence = "AAGCT"
    return Repeat(sequence=sequence)


@pytest.mark.django_db
def test_model_can_create_a_repeat(repeat):
    """Test the repeat model can create an Repeat."""
    old_count = Repeat.objects.count()
    repeat.save()
    new_count = Repeat.objects.count()
    assert old_count != new_count
