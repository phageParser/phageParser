from django.shortcuts import render
from rest_framework import viewsets
from restapi.serializers import RepeatSerializer, OrganismSerializer
from restapi.models import  Repeat, Organism
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
        Optionally restricts the returned purchases to a given accession,
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
