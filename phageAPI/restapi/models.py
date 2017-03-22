from __future__ import unicode_literals

from django.db import models

class Spacer(models.Model):
    sequence = models.CharField(max_length=2000, blank=True, default='')
class Repeat(models.Model):
    sequence = models.CharField(max_length=2000, blank=True, default='')
class AntiCRISPR(models.Model):
    sequence = models.CharField(max_length=2000, blank=True, default='')
class SpacerRepeatPair(models.Model):
    spacer = models.ForeignKey(Spacer, on_delete=models.CASCADE, null=True)
    repeat = models.ForeignKey(Repeat, on_delete=models.CASCADE, null=True)
class Organism(models.Model):
    name = models.CharField(max_length=50, blank=True, default='')
    accession = models.CharField(max_length=50, blank=True, default='')
    def __str__(self):
        return self.name
class OrganismSpacerRepeatPair(models.Model):
    # Has 1 organism and a spacerrepeatpair list with start and end values for locus
    organism = models.ForeignKey(Organism, on_delete=models.CASCADE, null=True)
    spacerrepeatpair = models.ForeignKey(SpacerRepeatPair, on_delete=models.CASCADE, null=True)
    order = models.PositiveIntegerField(default=0)
    genomic_start = models.PositiveIntegerField(default=0)
    genomic_end = models.PositiveIntegerField(default=0)
class OrganismAntiCRISPRPair(models.Model):
    organism = models.ForeignKey(Organism, on_delete=models.CASCADE, null=True)
    antiCRISPR = models.ForeignKey(AntiCRISPR, on_delete=models.CASCADE, null=True)
    genomic_start = models.PositiveIntegerField(default=0)
    genomic_end = models.PositiveIntegerField(default=0)
class Phage(models.Model):
    accession = models.CharField(max_length=50, blank=True, default='')
class PhageAntiCRISPRPair(models.Model):
    phage = models.ForeignKey(Phage, on_delete=models.CASCADE, null=True)
    antiCRISPR = models.ForeignKey(AntiCRISPR, on_delete=models.CASCADE, null=True)
    genomic_start = models.PositiveIntegerField(default=0)
    genomic_end = models.PositiveIntegerField(default=0)
class PhageSpacerPair(models.Model):
    phage = models.ForeignKey(Phage, on_delete=models.CASCADE, null=True)
    spacer = models.ForeignKey(Spacer, on_delete=models.CASCADE, null=True)
    genomic_start = models.PositiveIntegerField(default=0)
    genomic_end = models.PositiveIntegerField(default=0)
