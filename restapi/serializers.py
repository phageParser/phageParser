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


class GetCRISPRType(DynamicComputedField):
    def __init__(self, **kwargs):
        kwargs['field_type'] = str
        super(GetCRISPRType, self).__init__(**kwargs)

    def get_attribute(self, instance):
        proteinlist = instance.cas_proteins.values('type_specificity')
        return list(proteinlist)

    def to_representation(self, value):
        return set([b for a in value for b in a['type_specificity'].split(',')])


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
        fields = ('id', 'name', 'accession', 'cas_proteins', 'crisprtypes')
    crisprtypes = GetCRISPRType()
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
        fields = ('organism', 'genomic_start', 'genomic_end', 'spacerrepeats')
        deferred_fields = ('spacerrepeats',)
    organism = DynamicRelationField('OrganismSerializer')


class CasProteinSerializer(DynamicModelSerializer):
    class Meta:
        model = CasProtein
        fields = ('profileID', 'function', 'gene', 'group', 'type_specificity')


class OCSerializer(DynamicModelSerializer):
    class Meta:
        model = OrganismCasProtein
        fields = ('organism', 'casprotein',
                  'genomic_start', 'genomic_end', 'evalue')

    casprotein = DynamicRelationField('CasProteinSerializer')
    organism = DynamicRelationField('OrganismSerializer')
