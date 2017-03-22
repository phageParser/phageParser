from __future__ import unicode_literals

from django.db import models

import uuid
class Repeat(models.Model):
    id = models.UUIDField(primary_key=True, default=uuid.uuid4, editable=False)
    sequence = models.CharField(max_length=200, blank=True, default='')
    def __str__(self):
        return self.sequence
class Organism(models.Model):
    id = models.UUIDField(primary_key=True, default=uuid.uuid4, editable=False)
    name = models.CharField(max_length=50, blank=True, default='')
    accession = models.CharField(max_length=50, blank=True, default='')
    repeats = models.ManyToManyField(Repeat)
    def __str__(self):
        return self.name
class Spacer(models.Model):
    pass
class SpacerRepeatPair(models.Model):
    #Has 1 spacer and 1 repeat
    pass
class OrganismSpacerRepeatPair(models.Model):
    # Has 1 organism and a spacerrepeatpair list with start and end values for locus
    pass
