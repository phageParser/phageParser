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

def split_name(csv_list):
    csv_list[0].insert(1,"Accession") 
    for row in csv_list[1:]:
        phagename = row[1]
        row.insert(1,'')
        try:
            phagename.index('gi|')
            #print "gi Phage", phagename
            r = phagename.index('ref|')
            s = phagename[r+4:].index('|')
            refseq = phagename[r+4:r+4+s]
            row[1] = refseq
            n = phagename.index(' ')
            name = phagename[n+1:phagename.index(',')]
            row[2] = name
        except:
            pass
        try:
            phagename.index('ENA|')
            r = phagename[4:].index('|')
            n = phagename.index(' ')
            refseq = phagename[r+5:n]
            row[1] = refseq
            try:
                name = phagename[n+1:phagename.index(',')]
            except:
                name = phagename[n+1:]
            row[2] = name   
        except:
            pass
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
        row=[str(e) for e in row]
        print ('\t' .join(row))
    print (50*'_')

def write_csv(dest, csv_cont):
    with open(dest, 'w') as out_file:
        writer= csv.writer(out_file, delimiter=',')
        for row in csv_cont:
            writer.writerow(row)

# -------------------------------------------------------------------           
directory = "output"

for fn in os.listdir("%s/" %directory):
    if fn == 'sorted':
        continue
    csv_cont = csv_to_list("%s/%s" %(directory,fn))
    convert_cells_to_floats(csv_cont)
    csv_sorted = sort_by_column(csv_cont, "Expect")
    csv_sorted = split_name(csv_sorted)    
    csv_new = compare_phages(csv_sorted)
    write_csv("%s/sorted/sorted.%s" %(directory,fn), csv_new)
    #print_csv(csv_sorted)
    
