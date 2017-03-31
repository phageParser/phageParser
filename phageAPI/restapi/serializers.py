from restapi.models import Repeat, Organism
from rest_framework import serializers

class RepeatSerializer(serializers.ModelSerializer):
    class Meta:
        model = Repeat
        fields = ('id', 'sequence')
class OrganismSerializer(serializers.ModelSerializer):
    class Meta:
        model = Organism
        fields = ('id', 'name', 'accession', 'repeats')

