from restapi.serializers import SpacerSerializer, RepeatSerializer, OrganismSerializer, OSRSerializer, OCSerializer, CasProteinSerializer, LocusSerializer
from restapi.models import Spacer, Repeat, Organism, OrganismSpacerRepeat, OrganismCasProtein, CasProtein, Locus
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


class LocusViewSet(DynamicModelViewSet):
    """
    API endpoint that allows loci to be viewed or edited.
    """
    queryset = Locus.objects.all()
    serializer_class = LocusSerializer


class CasProteinViewSet(DynamicModelViewSet):
    """
    API endpoint that allows cas proteins to be viewed or edited.
    """
    queryset = CasProtein.objects.all()
    serializer_class = CasProteinSerializer


class OCViewSet(DynamicModelViewSet):
    """
    API endpoint that allows organism cas pairs to be viewed or edited.
    """
    queryset = OrganismCasProtein.objects.all()
    serializer_class = OCSerializer


class OSRViewSet(DynamicModelViewSet):
    """
    API endpoint that allows organism spacer repeat pairs to be viewed or edited.
    """
    queryset = OrganismSpacerRepeat.objects.all()
    serializer_class = OSRSerializer
