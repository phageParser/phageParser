#!/usr/bin/python
import sqlite3


# THIS file will parse in two files named 'DRdatabase' and 'Spacerdatabase',
# and populate the 'crispr.sqlite' file tables Spacers and Repeats, as per 
# ISSUE #99 https://github.com/goyalsid/phageParser/issues/99

def read_input_file(filename):
    fin = open(filename, 'r')
    flines = fin.readlines()
    file_dict = dict()
    for f in flines:
        f = f.replace('\n', '')
        if f[0] == '>':
            accessions = f[1:].split('|')
        else:
            file_dict[f] = accessions
    return file_dict


def print_dict(d):
    for i in d:
        print(i, d[i])


def get_largest_id(c, name, id_name):
    # NOT USED
    # SELECT * FROM Spacer ORDER BY SpacerID DESC LIMIT 1
    q = "SELECT * FROM " + name + " ORDER BY " + id_name + " DESC LIMIT 1"
    c.execute(q)
    return c.fetchone()[0]


def sql_add(elem, name, db_name):
    conn = sqlite3.connect(db_name)
    c = conn.cursor()
    col_name = name + 'Sequence'
    q = ("INSERT INTO " + str(name) + " (" + str(col_name)
         + ") VALUES ('" + elem + "')")
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
    q = "SELECT * FROM " + name + "  WHERE " + column + " = '" + elem + "'"
    # print q
    c.execute(q)
    r = c.fetchone()
    c.close()
    conn.close()
    return r


def populate_db(db_name, table_name, file_name):
    input_dict = read_input_file(file_name)
    for d in input_dict:
        sql_add(d, table_name, db_name)


repeat_file = 'data/DRdatabase.txt'
spacer_file = 'data/spacerdatabase.txt'
db_name = 'crispr.sqlite'

populate_db(db_name, 'Repeat', repeat_file)
populate_db(db_name, 'Spacer', spacer_file)
