## trimGenbankDNA - trims all the DNA sequence data from downloaded Genbank files to save space

import linecache

filename = "out3.gb"

# ----------------------------------------------------------------------------
i = filename.index('.')
nm = filename[0:i]

data = open(filename,"rb")
new_data = open("%s_cleaned.gb" %nm,"wb")

for num, line in enumerate(data,1):
    #print line

    if line[0:5] == 'LOCUS':
        new_data.write(line)
    else:
        continue
    
    check = True
    count = 1
    while check == True:
        genomeline = linecache.getline(filename,num+count)
        
        if genomeline[0:6] == 'ORIGIN':
            check = False
            break
        new_data.write(genomeline)
        count += 1
    
data.close()
new_data.close()

