from __future__ import unicode_literals

from django.db import models


class Spacer(models.Model):
    sequence = models.CharField(max_length=255, blank=True, default='')


class Repeat(models.Model):
    sequence = models.CharField(max_length=255, blank=True, default='')


class AntiCRISPR(models.Model):
    sequence = models.CharField(max_length=512, blank=True, default='')
    accession = models.CharField(max_length=32, blank=True, default='')


class Organism(models.Model):
    name = models.CharField(max_length=512, blank=True, default='')
    accession = models.CharField(max_length=32, blank=True, default='')

    def __str__(self):
        return self.name


class Phage(models.Model):
    accession = models.CharField(max_length=32, blank=True, default='')


class CasProtein(models.Model):
    profileID = models.CharField(max_length=32, blank=True, default='')
    function = models.CharField(max_length=32, blank=True, default='')
    gene = models.CharField(max_length=32, blank=True, default='')
    group = models.CharField(max_length=32, blank=True, default='')
    type_specificity = models.CharField(max_length=255, blank=True, default='')


class OrganismSpacerRepeatPair(models.Model):
    # Has 1 organism and a spacerrepeatpair list with start and end values in
    # locus
    organism = models.ForeignKey(Organism, on_delete=models.CASCADE, null=True)
    spacer = models.ForeignKey(Spacer, on_delete=models.CASCADE, null=True)
    repeat = models.ForeignKey(Repeat, on_delete=models.CASCADE, null=True)
    order = models.PositiveIntegerField(default=0)
    genomic_start = models.PositiveIntegerField(default=0)
    genomic_end = models.PositiveIntegerField(default=0)


class OrganismAntiCRISPRPair(models.Model):
    organism = models.ForeignKey(Organism, on_delete=models.CASCADE, null=True)
    antiCRISPR = models.ForeignKey(
        AntiCRISPR, on_delete=models.CASCADE, null=True)
    genomic_start = models.PositiveIntegerField(default=0)
    genomic_end = models.PositiveIntegerField(default=0)


class PhageAntiCRISPRPair(models.Model):
    phage = models.ForeignKey(Phage, on_delete=models.CASCADE, null=True)
    antiCRISPR = models.ForeignKey(
        AntiCRISPR, on_delete=models.CASCADE, null=True)
    genomic_start = models.PositiveIntegerField(default=0)
    genomic_end = models.PositiveIntegerField(default=0)


class PhageSpacerPair(models.Model):
    phage = models.ForeignKey(Phage, on_delete=models.CASCADE, null=True)
    spacer = models.ForeignKey(Spacer, on_delete=models.CASCADE, null=True)
    genomic_start = models.PositiveIntegerField(default=0)
    genomic_end = models.PositiveIntegerField(default=0)


class OrganismCasPair(models.Model):
    organism = models.ForeignKey(
        Organism, on_delete=models.CASCADE, null=True)
    casprotein = models.ForeignKey(
        CasProtein, on_delete=models.CASCADE, null=True)
    genomic_start = models.PositiveIntegerField(default=0)
    genomic_end = models.PositiveIntegerField(default=0)
    evalue = models.FloatField(default=0)
