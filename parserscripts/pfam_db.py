"""
pfam_db.py extracts and submits the translations from a GenBank input
file (such as data/Genbank_example.txt) to PFAM, before writing results
to a CSV file in the same directory.

Usage:
python pfam_db.py <infile>

Ref - PR #73
"""

import csv
import os
import re
import sys
import time
from xml.etree import ElementTree

import requests

file_path = str(sys.argv[1])
file_path_without_ext = os.path.splitext(file_path)[0]


def extract_cds(file):
    cds = []

    with open(file, 'r') as f:
        file_text = f.read()
    m = re.finditer('CDS\s{13}(?P<CDS>\S*)', file_text)

    for reg in m:
        element = reg.group('CDS')
        if element != "":
            cds.append(element)

    return cds


def extract_translations(file):
    translations = []

    with open(file, 'r') as f:
        file_text = f.read()
    m = re.finditer(
        '/translation="(?P<translation>[\w\n ]+)"',
        file_text,
        re.MULTILINE
    )

    for reg in m:
        translation = reg.group('translation')
        translation = translation.replace("\n", "")
        translation = translation.replace(" ", "")
        translations.append(translation)

    return translations


def get_translation_url(translation):
    url = "http://pfam.xfam.org/search/sequence"
    data = {'seq': translation, 'output': 'xml'}

    response = requests.get(url, data=data)

    tree = ElementTree.fromstring(response.content)

    return tree[0][1].text


def get_values(url):
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
    translations = extract_translations(file_path)
    cds = extract_cds(file_path)

    for i, translation in enumerate(translations):
        print("Getting translation " + str(i + 1) + "/"
              + str(len(translations)) + "...")
        url = get_translation_url(translation)
        result = get_values(url)

        for index_row, row in enumerate(result):
            row["translation_number"] = i + 1
            row["cds"] = cds[i]
            if i == 0 and index_row == 0:
                with open(file_path_without_ext + '.csv', 'w') as f:
                    w = csv.DictWriter(f, row.keys())
                    w.writeheader()

            with open(file_path_without_ext + '.csv', 'a') as f:
                w = csv.DictWriter(f, row.keys())
                w.writerow(row)


if __name__ == "__main__":
    main()
