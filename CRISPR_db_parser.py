import linecache

# CRISPR db parser
# MB Mar 07 2015

filename = "spacerdatabase.txt" # File from CRISPRdb to sort 
spacer_db = open("%s" %(filename),"r") 

# places to dump accession numbers during execution
refseq_list = []
refseq_dict = {}

for num, line in enumerate(spacer_db,1):
    check = True # awkward while loop
    if line[0] == ">": # use the headers, indicated by >, to sort
        line = line[1:] # delete 1st character to make loop same each time around
        counter = 0
        while check == True: 
            counter +=1
            refseq = line[0:9] # this part of the header is the NCBI accession#
            if refseq not in refseq_list: # open new file if it's a new bacteria
                refseq_dict[refseq] = open("%s.fasta" %(refseq),"w")
                if "|" in line: # if more than one bacteria contain spacer
                    i = line.index("|")
                    writeline = line[10:i]  # include in header the locus identifier and spacer position identifier
                    writeline2 = writeline.replace('_','.')
                else: # if it's only one bacteria
                    writeline = line[10:] 
                    writeline2 = writeline.replace('_','.')
                # write header and spacer to file
                refseq_dict[refseq].write(">" + writeline2 + "\n")
                refseq_dict[refseq].write(linecache.getline("%s" %(filename), num+1))
                if counter == 1: # since the file is organized alphabetically by the first bacteria in the header, if we see a different first bacteria we can close the previous file to free up space. This might be buggy. 
                    try:
                        refseq_prev = linecache.getline("%s" %(filename), num-2)[1:10]
                        refseq_dict[refseq_prev].close() 
                    except:
                        pass # throws exception on the first time through, otherwise wouldn't
                refseq_list.append(refseq) 
            if refseq in refseq_list:
                if "|" in line:
                    i = line.index("|")
                    writeline = line[10:i]  # include in header the locus identifier and spacer position identifier
                    writeline2 = writeline.replace('_','.')
                else:
                    writeline = line[10:] 
                    writeline2 = writeline.replace('_','.')
                
                refseq_dict[refseq].write(">" + writeline2 + "\n")
                refseq_dict[refseq].write(linecache.getline("%s" %(filename), num+1))
            try:
                i = line.index("|")
                line = line[i+1:] # change the header so that the next bacteria is up for the loop
            except:
                check = False

for key in refseq_dict:
    if not refseq_dict[key].closed:
        refseq_dict[key].close()

spacer_db.close()
