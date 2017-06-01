from restapi.models import Spacer, Repeat, Organism, OrganismSpacerRepeatPair, OrganismCas
from rest_framework import serializers
from dynamic_rest.serializers import DynamicModelSerializer, DynamicRelationField
from dynamic_rest.fields import DynamicComputedField


class GetSequenceLength(DynamicComputedField):
    def __init__(self, **kwargs):
        kwargs['field_type'] = int
        super(GetSequenceLength, self).__init__(**kwargs)

    def get_attribute(self, instance):
        return len(instance.sequence)

    def to_representation(self, value):
        return int(value)


class SequenceSerializer(DynamicModelSerializer):
    length = GetSequenceLength()


class SpacerSerializer(SequenceSerializer):
    class Meta:
        model = Spacer
        fields = ('id', 'sequence', 'length', 'sequence')


class RepeatSerializer(SequenceSerializer):
    class Meta:
        model = Repeat
        fields = ('id', 'length', 'sequence')


class OrganismSerializer(DynamicModelSerializer):
    class Meta:
        model = Organism
        fields = ('id', 'name', 'accession')


class OrganismCasSerializer(DynamicModelSerializer):
    class Meta:
        model = OrganismCas
        fields = ('organism_id', 'accession')


class OSRPairSerializer(DynamicModelSerializer):
    spacer = DynamicRelationField('SpacerSerializer')
    repeat = DynamicRelationField('RepeatSerializer')
    organism = DynamicRelationField('OrganismSerializer')

    class Meta:
        model = OrganismSpacerRepeatPair
        fields = ('id', 'organism', 'spacer', 'repeat',
                  'order', 'genomic_start', 'genomic_end')
