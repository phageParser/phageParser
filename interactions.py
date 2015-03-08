## Making interactions.txt file from sorted BLAST outputs

import json,os,csv

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
            


elements = []
phages = []
edges = []


for fn in os.listdir("output/sorted"):
    # made node for each bacteria and enter in elements
    bac_name = fn[7:16]
    obj = { 
        "group": 'nodes',
        "data": {"id":bac_name},
        "classes": 'bacteria'
    }
    elements.append(obj)
    csv_cont = csv_to_list("output/sorted/"+fn)[1:]
    for i in range(len(csv_cont)):
        #only make phage node if it doesn't already exists
        if (not(csv_cont[i][1] in phages)):
            obj = {
                "group": 'nodes',
                "data": { "id": csv_cont[i][1]},
                "classes": 'phage'
            }
            elements.append(obj)
            phages.append(csv_cont[i][1])
    
    for i in range(len(csv_cont)):
        if (not(bac_name+"_to_"+csv_cont[i][1] in edges)):
            obj ={
                "group": 'edges',
                "data": {"id": bac_name+"_to_"+csv_cont[i][1],
                         "source": bac_name,
                         "target": csv_cont[i][1]}
            }
            elements.append(obj)
            edges.append(bac_name+"_to_"+csv_cont[i][1])
            
            

    # open each sorted file in output/sorted and get the phages, if phage is not phages list, add and create node, create edge
    



"""
json_file = open("json.txt","w")

json_file.write("cytoscape({" +"\n")
json_file.write("container: document.getElementById('cy')," +"\n")

json_file.write("elements: [" + "\n")
"""
with open("json.txt","w") as outfile:
    json.dump(elements,outfile)
    
#for thing in elements:
#    json.write(thing)
    
