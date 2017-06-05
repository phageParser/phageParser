'''
pfam_db.py extracts and submits the translations from a GenBank input file (such as data/Genbank_example.txt) 
to PFAM, before writing results to a CSV file in the same directory.

Usage:
python pfam_db.py <infile> 

Ref - PR #73
'''

import re
import requests
from xml.etree import ElementTree
import time
import csv
import os
import sys

filePath = str(sys.argv[1])
filePathWithoutExt = os.path.splitext(filePath)[0]

def extractCDS(file):
    cds = []

    f = open(file, 'r')
    fileText = f.read()
    m = re.finditer('CDS\s{13}(?P<CDS>\S*)', fileText)

    for reg in m:
        element = reg.group('CDS')
        if element != "":
            cds.append(element)

    return cds

def extractTranslations(file):
    translations = []

    f = open(file, 'r')
    fileText = f.read()
    m = re.finditer('/translation="(?P<translation>[\w\n ]+)"', fileText, re.MULTILINE)

    for reg in m:
        translation = reg.group('translation')
        translation = translation.replace("\n", "")
        translation = translation.replace(" ", "")
        translations.append(translation)

    return translations

def getTranslationURL(translation):
    url = "http://pfam.xfam.org/search/sequence"
    data = {'seq': translation, 'output': 'xml'}

    response = requests.get(url, data = data)

    tree = ElementTree.fromstring(response.content)

    return tree[0][1].text

def getValues(url):
    sleep_time = 1
    tries = 60
    xml = ""
    result = []

    for ntry in range(tries):
        response = requests.get(url)

        if response.status_code != 202 and response.status_code != 200:
            print("Error: " + response.status_code)
            break

        if response.status_code == 200:
            xml = response.content
            if xml:
                break

        time.sleep(sleep_time)

    if ntry + 1 == tries:
        print("Reached max tries without response.")
        return


    try:
        tree = ElementTree.fromstring(xml)

        for i, match in enumerate(tree[0][0][0][0]):
            result.append(dict(match.attrib, **match[0].attrib))
    except:
        pass

    return result

def main():
    translations = extractTranslations(filePath)
    cds = extractCDS(filePath)

    for i, translation in enumerate(translations):
        print("Getting translation " + str(i + 1) + "/" + str(len(translations)) + "...")
        url = getTranslationURL(translation)
        result = getValues(url)

        for indexRow, row in enumerate(result):
            row["translation_number"] = i + 1
            row["cds"] = cds[i]
            if i == 0 and indexRow == 0:
                with open(filePathWithoutExt + '.csv', 'w') as f:
                    w = csv.DictWriter(f, row.keys())
                    w.writeheader()

            with open(filePathWithoutExt + '.csv', 'a') as f:
                w = csv.DictWriter(f, row.keys())
                w.writerow(row)

if __name__ == "__main__":
    main()
