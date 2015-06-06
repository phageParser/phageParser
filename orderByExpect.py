import time, os, operator, csv


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

def print_csv(csv_content):
    print (50*'-')
    for row in csv_content:
        row=[str(e) for e in row]
        print ('\t' .join(row))
    print (50*'_')

def write_csv(dest, csv_cont):
    with open(dest, 'w') as out_file:
        writer= csv.writer(out_file, delimiter=',')
        for row in csv_cont:
            writer.writerow(row)
            

for fn in os.listdir("output/"):
    csv_cont = csv_to_list('output/'+fn)
    convert_cells_to_floats(csv_cont)
    csv_sorted = sort_by_column(csv_cont, "Expect")
    write_csv("output/sorted/"+"sorted."+fn, csv_sorted)
    #print_csv(csv_sorted)
    
