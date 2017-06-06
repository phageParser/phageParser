from restapi.models import Spacer, Repeat, Organism, OrganismSpacerRepeatPair, CasProtein, OrganismCasPair
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


class GetCasProteins(DynamicComputedField):
    def __init__(self, **kwargs):
        kwargs['field_type'] = int
        super(GetCasProteins, self).__init__(**kwargs)

    def get_attribute(self, instance):
        pairlist = OrganismCasPair.objects.filter(
            organism=instance).values('casprotein')
        return list(CasProtein.objects.filter(id__in=pairlist))

    def to_representation(self, value):
        return [a.id for a in value]


class GetCRISPRType(DynamicComputedField):
    def __init__(self, **kwargs):
        kwargs['field_type'] = str
        super(GetCRISPRType, self).__init__(**kwargs)

    def get_attribute(self, instance):
        pairlist = OrganismCasPair.objects.filter(
            organism=instance).values('casprotein')
        proteinlist = CasProtein.objects.filter(id__in=pairlist).values('type_specificity')
        return list(proteinlist)

    def to_representation(self, value):
        return set([b for a in value for b in a['type_specificity'].split(',') ])
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
    casproteins = GetCasProteins()
    crisprtypes = GetCRISPRType()

    class Meta:
        model = Organism
        fields = ('id', 'name', 'accession', 'casproteins', 'crisprtypes')


class OCPairSerializer(DynamicModelSerializer):
    class Meta:
        model = OrganismCasPair
        fields = ('organism', 'casprotein',
                  'genomic_start', 'genomic_end', 'evalue')


class CasProteinSerializer(DynamicModelSerializer):
    class Meta:
        model = CasProtein
        fields = ('profileID', 'function', 'gene', 'group', 'type_specificity')


class OSRPairSerializer(DynamicModelSerializer):
    spacer = DynamicRelationField('SpacerSerializer')
    repeat = DynamicRelationField('RepeatSerializer')
    organism = DynamicRelationField('OrganismSerializer')

    class Meta:
        model = OrganismSpacerRepeatPair
        fields = ('id', 'organism', 'spacer', 'repeat',
                  'order', 'genomic_start', 'genomic_end')
