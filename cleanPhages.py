import time, os, operator, csv
from models import Phage
from parsers.find_accession import PhageFinder

def csv_to_list(csv_file, delimiter=","):
    with open (csv_file, 'r') as csv_con:
        reader = csv.reader(csv_con, delimiter=delimiter)
        return list(reader)
    
def convert_cells_to_floats(csv_cont):
    for row in range(len(csv_cont)):
        for cell in range(len(csv_cont[row])):
            try:
                csv_cont[row][cell] = float(csv_cont[row][cell])
            except ValueError:
                pass

def sort_by_column(csv_cont, col, reverse=False):
    header = csv_cont[0]
    body = csv_cont[1:]
    if isinstance(col,str):
        col_index = header.index(col)
    else:
        col_index = col
    body = sorted(body,
                  key=operator.itemgetter(col_index),
                  reverse=reverse)
    body.insert(0,header)
    return body

def get_phagename_and_refseq(row, phage_finder):
    phage = Phage(row, phage_finder)
    return phage.name, phage.refseq

def split_name(csv_list, phage_finder):
    csv_list[0].insert(1,"Accession")
    for row in csv_list[1:]:
        row.insert(1,'')
        name, refseq = get_phagename_and_refseq(row[1], phage_finder)
        if name is not None and refseq is not None:
            row[1], row[2] = name, ref_seq
    return csv_list

def compare_phages(csv_sorted):
    csv_new = []
    datas = [] 
    for i in range(len(csv_sorted)):
        if csv_sorted[i][3:] in datas:
            continue
        datas.append(csv_sorted[i][3:])
        csv_new.append(csv_sorted[i])
    return csv_new


def print_csv(csv_content):
    print (50*'-')
    for row in csv_content:
        row = [str(e) for e in row]
        print ('\t' .join(row))
    print (50*'_')

def write_csv(dest, csv_cont):
    with open(dest, 'w+') as out_file:
        writer = csv.writer(out_file, delimiter=',')
        for row in csv_cont:
            writer.writerow(row)

# -------------------------------------------------------------------           
directory = "output"
phage_finder = PhageFinder('data/PhagesDB_Data.txt')

sorted_dir = os.path.dirname("%s/" %directory) + "/sorted"
if not os.path.isdir(sorted_dir):
    os.makedirs(sorted_dir)

for fn in os.listdir("%s/" %directory):
    if fn == 'sorted':
        continue
    csv_cont = csv_to_list("%s/%s" %(directory,fn))
    convert_cells_to_floats(csv_cont)
    csv_sorted = sort_by_column(csv_cont, "Expect")
    csv_sorted = split_name(csv_sorted, phage_finder)
    csv_new = compare_phages(csv_sorted)
    write_csv("%s/sorted/sorted.%s" %(directory,fn), csv_new)
    #print_csv(csv_sorted)
    
