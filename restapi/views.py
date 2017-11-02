from dynamic_rest.viewsets import DynamicModelViewSet

from restapi.models import (
    CasProtein,
    Locus,
    LocusSpacerRepeat,
    Organism,
    OrganismCasProtein,
    OrganismSelfSpacer,
    Repeat,
    Spacer
)
from restapi.serializers import (
    CasProteinSerializer,
    LSRSerializer,
    LocusSerializer,
    OCSerializer,
    OSSSerializer,
    OrganismSerializer,
    SpacerSerializer,
    RepeatSerializer
)


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


class OSSViewSet(DynamicModelViewSet):
    """
    API endpoint that allows organism self targeting spacer pairs to be
    viewed or edited.
    """
    queryset = OrganismSelfSpacer.objects.all()
    serializer_class = OSSSerializer


class LSRViewSet(DynamicModelViewSet):
    """
    API endpoint that allows locus spacer repeat pairs to be viewed or
    edited.
    """
    queryset = LocusSpacerRepeat.objects.all()
    serializer_class = LSRSerializer
