from django.shortcuts import render
from rest_framework import viewsets
from restapi.serializers import SpacerSerializer, RepeatSerializer, OrganismSerializer, OSRPairSerializer
from restapi.models import Spacer, Repeat, Organism, OrganismSpacerRepeatPair


class SpacerViewSet(viewsets.ModelViewSet):
    """
    API endpoint that allows spacers to be viewed or edited.
    """
    queryset = Spacer.objects.all()
    serializer_class = SpacerSerializer


class RepeatViewSet(viewsets.ModelViewSet):
    """
    API endpoint that allows repeats to be viewed or edited.
    """
    queryset = Repeat.objects.all()
    serializer_class = RepeatSerializer


class OrganismViewSet(viewsets.ModelViewSet):
    """
    API endpoint that allows organisms to be viewed or edited.
    """
    queryset = Organism.objects.all()
    serializer_class = OrganismSerializer

    def get_queryset(self):
        """
        Optionally restricts the returned organisms to a given accession,
        by filtering against a `accession` query parameter in the URL.
        """
        queryset = Organism.objects.all()
        accession = self.request.query_params.get('accession', None)
        name = self.request.query_params.get('name', None)
        if accession is not None:
            queryset = queryset.filter(accession=accession)
        if name is not None:
            queryset = queryset.filter(name=name)
        return queryset


class OSRPairViewSet(viewsets.ModelViewSet):
    """
    API endpoint that allows repeats to be viewed or edited.
    """
    queryset = OrganismSpacerRepeatPair.objects.all()
    serializer_class = OSRPairSerializer
