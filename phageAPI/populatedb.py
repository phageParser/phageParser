from restapi.models import Repeat, Organism
from restapi.serializers import RepeatSerializer, OrganismSerializer
from rest_framework.renderers import JSONRenderer
from rest_framework.parsers import JSONParser

a1 = Repeat(sequence='A1')
a1.save()
a2 = Repeat(sequence='A2')
a2.save()
serializer = RepeatSerializer(a1)
serializer.data
