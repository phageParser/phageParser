from __future__ import unicode_literals

from django.db import models


class Spacer(models.Model):
    sequence = models.CharField(max_length=255, blank=True, default='')
    loci = models.ManyToManyField('Locus', through='LocusSpacerRepeat')
    repeats = models.ManyToManyField('Repeat', through='LocusSpacerRepeat')


class Repeat(models.Model):
    sequence = models.CharField(max_length=255, blank=True, default='')
    loci = models.ManyToManyField('Locus', through='LocusSpacerRepeat')
    spacers = models.ManyToManyField('Spacer', through='LocusSpacerRepeat')


class AntiCRISPR(models.Model):
    sequence = models.CharField(max_length=512, blank=True, default='')
    accession = models.CharField(max_length=32, blank=True, default='')


class CasProtein(models.Model):
    profileID = models.CharField(max_length=32, blank=True, default='')
    function = models.CharField(max_length=32, blank=True, default='')
    gene = models.CharField(max_length=32, blank=True, default='')
    group = models.CharField(max_length=32, blank=True, default='')
    type_specificity = models.CharField(max_length=255, blank=True, default='')


class Organism(models.Model):
    name = models.CharField(max_length=512, blank=True, default='')
    accession = models.CharField(max_length=32, blank=True, default='')
    cas_proteins = models.ManyToManyField(
        CasProtein, through='OrganismCasProtein')
    self_spacers = models.ManyToManyField(Spacer, through='OrganismSelfSpacer')
    def __str__(self):
        return '{} with accession {}'.format(self.name, self.accession)


class Locus(models.Model):
    organism = models.ForeignKey(Organism, on_delete=models.CASCADE, null=True)
    genomic_start = models.PositiveIntegerField(default=0)
    genomic_end = models.PositiveIntegerField(default=0)
    spacers = models.ManyToManyField(Spacer, through='LocusSpacerRepeat')
    repeats = models.ManyToManyField(Repeat, through='LocusSpacerRepeat')


class Phage(models.Model):
    accession = models.CharField(max_length=32, blank=True, default='')


class LocusSpacerRepeat(models.Model):
    locus = models.ForeignKey(
        Locus, on_delete=models.CASCADE, null=True, related_name='spacerrepeats')
    spacer = models.ForeignKey('Spacer', on_delete=models.CASCADE, null=True)
    repeat = models.ForeignKey('Repeat', on_delete=models.CASCADE, null=True)
    order = models.PositiveIntegerField(default=0)


class OrganismAntiCRISPR(models.Model):
    organism = models.ForeignKey('Organism', on_delete=models.CASCADE, null=True)
    antiCRISPR = models.ForeignKey(
        'AntiCRISPR', on_delete=models.CASCADE, null=True)
    genomic_start = models.PositiveIntegerField(default=0)
    genomic_end = models.PositiveIntegerField(default=0)


class PhageAntiCRISPR(models.Model):
    phage = models.ForeignKey('Phage', on_delete=models.CASCADE, null=True)
    antiCRISPR = models.ForeignKey(
        'AntiCRISPR', on_delete=models.CASCADE, null=True)
    genomic_start = models.PositiveIntegerField(default=0)
    genomic_end = models.PositiveIntegerField(default=0)


class PhageSpacer(models.Model):
    phage = models.ForeignKey('Phage', on_delete=models.CASCADE, null=True)
    spacer = models.ForeignKey('Spacer', on_delete=models.CASCADE, null=True)
    genomic_start = models.PositiveIntegerField(default=0)
    genomic_end = models.PositiveIntegerField(default=0)

class OrganismSelfSpacer(models.Model):
    organism = models.ForeignKey(
        'Organism', on_delete=models.CASCADE, null=True)
    spacer = models.ForeignKey(
        'Spacer', on_delete=models.CASCADE, null=True)
    genomic_start = models.PositiveIntegerField(default=0)
    genomic_end = models.PositiveIntegerField(default=0)
    evalue = models.FloatField(default=0)
class OrganismCasProtein(models.Model):
    organism = models.ForeignKey(
        'Organism', on_delete=models.CASCADE, null=True)
    casprotein = models.ForeignKey(
        'CasProtein', on_delete=models.CASCADE, null=True)
    genomic_start = models.PositiveIntegerField(default=0)
    genomic_end = models.PositiveIntegerField(default=0)
    evalue = models.FloatField(default=0)
