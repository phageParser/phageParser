import glob
import json
import os


# !/usr/bin/env python3
# -*- coding: utf-8 -*-


def main(blastdir):
    print('Populating organism selftargeting spacer objects')
    for fn in glob.glob(blastdir + '/*.json'):
        # get loci intervals of organism
        accession = os.path.splitext(os.path.split(fn)[1])[0]
        q_org = Organism.objects.filter(accession=accession)
        if not q_org.exists():
            print('Organism with accession {} is not in db but blast '
                  'report exists'.format(accession))
            continue
        org = q_org[0]
        interval_loci = [
            (entry['genomic_start'], entry['genomic_end'])
            for entry in org.locus_set.all().values('genomic_start',
                                                    'genomic_end')
        ]
        with open(fn, 'r') as f:
            try:
                blastrec = json.loads(f.read())
            except Exception as e:
                print('Error on accession {}\n{}'.format(accession, e))
                continue
        for res in blastrec['BlastOutput2']:
            query = res['report']['results']['bl2seq'][0]
            spacerid = query['query_title']
            for hit in query['hits']:
                for hsps in hit['hsps']:
                    start_h, end_h = hsps['hit_from'], hsps['hit_to']
                    in_locus = any([start_h > start and end_h < end
                                    for start, end in interval_loci])
                    if in_locus:
                        continue
                    q_spacer = Spacer.objects.filter(id=int(spacerid))
                    if not q_spacer.exists():
                        print('Spacer with sequence {} for organism {} '
                              'not found in db'.format(hsps['qseq'],
                                                       org.accession))
                        continue
                    spacer = q_spacer[0]
                    evalue = float(hsps['evalue'])
                    oselftarget, _ = OrganismSelfSpacer.objects.get_or_create(
                        organism=org,
                        spacer=spacer,
                        evalue=evalue,
                        genomic_start=start_h,
                        genomic_end=end_h
                    )


if __name__ == '__main__':
    import django

    os.environ.setdefault("DJANGO_SETTINGS_MODULE", "phageAPI.settings")
    django.setup()
    from restapi.models import Organism, OrganismSelfSpacer, Spacer

    main('gbfiles/blastoutput')
