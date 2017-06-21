from restapi.models import Spacer, Repeat, Organism, LocusSpacerRepeat, CasProtein, OrganismCasProtein, Locus
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
        fields = ('id', 'sequence', 'length', 'sequence', 'loci', 'repeats')
        deferred_fields = ('loci', 'repeats')
    loci = DynamicRelationField('LocusSerializer', many=True, embed=True)
    repeats = DynamicRelationField('RepeatSerializer', many=True)


class RepeatSerializer(SequenceSerializer):
    class Meta:
        model = Repeat
        fields = ('id', 'length', 'sequence', 'loci', 'spacers')
        deferred_fields = ('loci', 'spacers')
    loci = DynamicRelationField('LocusSerializer', many=True)
    spacers = DynamicRelationField('SpacerSerializer', many=True)


class OrganismSerializer(DynamicModelSerializer):
    class Meta:
        model = Organism
        fields = ('id', 'name', 'accession', 'cas_proteins',
                  'loci')
        deferred_fields = ('cas_proteins',
                           'loci')
    loci = DynamicRelationField(
        'LocusSerializer', source='locus_set', embed=True, many=True)
    cas_proteins = DynamicRelationField('CasProteinSerializer', many=True)


class LSRSerializer(DynamicModelSerializer):
    class Meta:
        model = LocusSpacerRepeat
        fields = ('id', 'locus', 'spacer', 'repeat',
                  'order')
    spacer = DynamicRelationField('SpacerSerializer')
    repeat = DynamicRelationField('RepeatSerializer')
    locus = DynamicRelationField('LocusSerializer')


class LocusSerializer(DynamicModelSerializer):

    class Meta:
        model = Locus
        fields = ('id', 'organism', 'genomic_start', 'genomic_end',
                  'spacerrepeats', 'spacers', 'repeats')
        deferred_fields = ('spacerrepeats', 'spacers', 'repeats')
    organism = DynamicRelationField('OrganismSerializer')
    spacers = DynamicRelationField('SpacerSerializer', embed=True, many=True)
    repeats = DynamicRelationField('RepeatSerializer', embed=True, many=True)


class CasProteinSerializer(DynamicModelSerializer):
    class Meta:
        model = CasProtein
        fields = ('id', 'profileID', 'function',
                  'gene', 'group', 'type_specificity')


class OCSerializer(DynamicModelSerializer):
    class Meta:
        model = OrganismCasProtein
        fields = ('id', 'organism', 'casprotein',
                  'genomic_start', 'genomic_end', 'evalue')

    casprotein = DynamicRelationField('CasProteinSerializer')
    organism = DynamicRelationField('OrganismSerializer')
