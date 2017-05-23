from restapi.models import Spacer, Repeat, Organism, OrganismSpacerRepeatPair
from rest_framework import serializers


class SpacerSerializer(serializers.ModelSerializer):
    class Meta:
        model = Spacer
        fields = ('id', 'sequence')


class RepeatSerializer(serializers.ModelSerializer):
    class Meta:
        model = Repeat
        fields = ('id', 'sequence')


class OrganismSerializer(serializers.ModelSerializer):
    class Meta:
        model = Organism
        fields = ('id', 'name', 'accession')


class OSRPairSerializer(serializers.ModelSerializer):
    class Meta:
        model = OrganismSpacerRepeatPair
        fields = ('id', 'organism', 'spacer', 'repeat',
                  'order', 'genomic_start', 'genomic_end')
