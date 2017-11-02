#!/usr/bin/python
from __future__ import print_function

import sqlite3


def print_dict(d):
    for i in d:
        print(i, d[i])


def sql_add(elem, name, columns, db_name):
    conn = sqlite3.connect(db_name)
    c = conn.cursor()
    if isinstance(elem, list):
        # It is assumed that the length of list and of columns will
        # always match
        q = ("INSERT INTO " + str(name) + " (" + ",".join(columns)
             + ") VALUES ('" + "','".join(elem) + "')")
    else:
        q = "INSERT INTO " + str(name) + columns + " VALUES " + elem
    # print q
    try:
        c.execute(q)
        conn.commit()
    except sqlite3.IntegrityError:
        print("Already In Database", elem)
    c.close()
    conn.close()


def sql_search(elem, name, column, db_name):
    conn = sqlite3.connect(db_name)
    c = conn.cursor()
    if isinstance(elem, list):
        # It is assumed that the length of list and of columns will
        # always match
        q = "SELECT * FROM " + name + " " + " WHERE "
        for i in range(0, len(elem)):
            q += str(column[i]) + " = " + str(elem[i]) + " and "
        q = q[:-5]
    else:
        q = "SELECT * FROM " + name + " " + " WHERE " + column + " = " + elem
    # print q
    c.execute(q)
    r = c.fetchone()
    c.close()
    conn.close()
    return r


def read_repeat_file(filename):
    fin = open(filename, 'r')
    out_dict = {}
    flines = fin.readlines()
    it = iter(flines)
    tuple_list = zip(it, it)
    for entry in tuple_list:
        accessions = entry[0][1:].replace('\n', '').split('|')
        for acc in accessions:
            out_dict[acc] = entry[1].replace('\n', '')
    return out_dict


def get_id_and_add(elem, name, column, db_name):
    item_id = sql_search(elem, name, column, db_name)
    if item_id is None:
        sql_add(elem, name, column, db_name)
        item_id = sql_search(elem, name, column, db_name)
    return item_id[0]


def match_repeat_to_spacer(repeat_data, spacer_file_name, db_name):
    spacer_file = open(spacer_file_name, 'r').readlines()
    it = iter(spacer_file)
    tuple_list = zip(it, it)
    for entry in tuple_list:
        accessions = entry[0][1:].replace('\n', '').split('|')
        sequence = entry[1].replace('\n', '')
        # print entry[1]
        for acc in accessions:
            acc_match = "_".join(acc.split("_")[:-1])
            spacer_id = get_id_and_add(
                "('" + sequence + "')",
                'Spacer',
                '(SpacerSequence)',
                db_name
            )
            repeat_id = get_id_and_add(
                "('" + repeat_data[acc_match] + "')",
                'Repeat',
                '(RepeatSequence)',
                db_name
            )
            try:
                match = repeat_data[acc_match]
            except KeyError:
                print("Error: Wrong Accession code")
            get_id_and_add(
                [str(spacer_id), str(repeat_id)],
                'SpacerRepeatPair',
                ['SpacerID', 'RepeatID'],
                db_name
            )


spacer_file = 'data/spacerdatabase.txt'
repeat_file = 'data/DRdatabase.txt'
db_name = 'crispr.sqlite'
repeat_dict = read_repeat_file(repeat_file)
match_repeat_to_spacer(repeat_dict, spacer_file, db_name)
