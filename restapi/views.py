from restapi.serializers import SpacerSerializer, RepeatSerializer, OrganismSerializer, OSRPairSerializer, OCPairSerializer, CasProteinSerializer
from restapi.models import Spacer, Repeat, Organism, OrganismSpacerRepeatPair, OrganismCasPair, CasProtein
from dynamic_rest.viewsets import DynamicModelViewSet


class SpacerViewSet(DynamicModelViewSet):
    """
    API endpoint that allows spacers to be viewed or edited.
    """
    queryset = Spacer.objects.all()
    serializer_class = SpacerSerializer


class RepeatViewSet(DynamicModelViewSet):
    """
    API endpoint that allows repeats to be viewed or edited.
    """
    queryset = Repeat.objects.all()
    serializer_class = RepeatSerializer


class OrganismViewSet(DynamicModelViewSet):
    """
    API endpoint that allows organisms to be viewed or edited.
    """
    queryset = Organism.objects.all()
    serializer_class = OrganismSerializer


class CasProteinViewSet(DynamicModelViewSet):
    """
    API endpoint that allows cas proteins to be viewed or edited.
    """
    queryset = CasProtein.objects.all()
    serializer_class = CasProteinSerializer


class OCPairViewSet(DynamicModelViewSet):
    """
    API endpoint that allows organism cas pairs to be viewed or edited.
    """
    queryset = OrganismCasPair.objects.all()
    serializer_class = OCPairSerializer


class OSRPairViewSet(DynamicModelViewSet):
    """
    API endpoint that allows organism spacer repeat pairs to be viewed or edited.
    """
    queryset = OrganismSpacerRepeatPair.objects.all()
    serializer_class = OSRPairSerializer
