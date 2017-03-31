import os

def populate_organism():
    def add_organism(name, accession):
        o, created= Organism.objects.get_or_create(name=name, accession=accession) #gets the object, this also checks for duplicates
        return o
    add_organism(name='BLA2', accession='BLA2')
if __name__ == '__main__':
    print "Starting orgnanism population script"
    os.environ.setdefault('DJANGO_SETTINGS_MODULE', 'phageAPI.settings')
    import django
    django.setup()
    from restapi.models import Organism
    populate_organism()
