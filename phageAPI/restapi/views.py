from django.shortcuts import render
from rest_framework import viewsets
from restapi.serializers import SpacerSerializer, RepeatSerializer, OrganismSerializer, OSRPairSerializer, OrganismCasSerializer
from restapi.models import Spacer, Repeat, Organism, OrganismSpacerRepeatPair, OrganismCas
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


class OrganismCasViewSet(DynamicModelViewSet):
    """
    API endpoint that allows spacers to be viewed or edited.
    """
    queryset = OrganismCas.objects.all()
    serializer_class = OrganismCasSerializer


class OSRPairViewSet(DynamicModelViewSet):
    """
    API endpoint that allows repeats to be viewed or edited.
    """
    queryset = OrganismSpacerRepeatPair.objects.all()
    serializer_class = OSRPairSerializer
