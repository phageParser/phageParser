from restapi.models import Spacer, Repeat, Organism, OrganismSpacerRepeatPair
from rest_framework import serializers


class SpacerSerializer(serializers.HyperlinkedModelSerializer):
    length = serializers.SerializerMethodField()

    def get_length(self, obj):
        return len(obj.sequence)

    class Meta:
        model = Spacer
        fields = ('url', 'sequence', 'length')


class RepeatSerializer(serializers.HyperlinkedModelSerializer):
    length = serializers.SerializerMethodField()

    def get_length(self, obj):
        return len(obj.sequence)

    class Meta:
        model = Repeat
        fields = ('url', 'sequence', 'length')


class OrganismSerializer(serializers.HyperlinkedModelSerializer):
    class Meta:
        model = Organism
        fields = ('url', 'name', 'accession')


class OSRPairSerializer(serializers.HyperlinkedModelSerializer):
    spacer = serializers.HyperlinkedRelatedField(
        many=False, read_only=True, view_name='spacer-detail')
    repeat = serializers.HyperlinkedRelatedField(
        many=False, read_only=True, view_name='repeat-detail')
    organism = serializers.HyperlinkedRelatedField(
        many=False, read_only=True, view_name='organism-detail')

    class Meta:
        model = OrganismSpacerRepeatPair
        fields = ('id', 'organism', 'spacer', 'repeat',
                  'order', 'genomic_start', 'genomic_end')
