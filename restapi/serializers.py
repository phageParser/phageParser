from restapi.models import Spacer, Repeat, Organism, OrganismSpacerRepeat, CasProtein, OrganismCasProtein, Locus
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
    crisprtypes = GetCRISPRType()

    class Meta:
        model = Organism
        fields = ('id', 'name', 'accession', 'cas_proteins', 'crisprtypes')


class LocusSerializer(DynamicModelSerializer):
    organism = DynamicRelationField('Organism')
    class Meta:
        model = Locus
        fields = ('organism', 'genomic_start', 'genomic_end')


class CasProteinSerializer(DynamicModelSerializer):
    class Meta:
        model = CasProtein
        fields = ('profileID', 'function', 'gene', 'group', 'type_specificity')


class OCSerializer(DynamicModelSerializer):
    casprotein = DynamicRelationField('CasProteinSerializer')
    organism = DynamicRelationField('Organism')

    class Meta:
        model = OrganismCasProtein
        fields = ('organism', 'casprotein',
                  'genomic_start', 'genomic_end', 'evalue')


class OSRSerializer(DynamicModelSerializer):
    spacer = DynamicRelationField('SpacerSerializer')
    repeat = DynamicRelationField('RepeatSerializer')
    locus = DynamicRelationField('LocusSerializer')

    class Meta:
        model = OrganismSpacerRepeat
        fields = ('id', 'locus', 'spacer', 'repeat',
                  'order')
