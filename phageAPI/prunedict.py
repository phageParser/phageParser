import pickle
import collections
from IPython import embed
def findincompleterecords(gendict):
    incomplist = []
    for locid in gendict:
        locus = gendict[locid]
        if not set(locus.keys()) == {'Spacers', 'Start', 'RepeatSeq', 'Stop'}:
            incomplist.append(locid)
    return incomplist
def findsizeoffsets(gendict):
    offsetlist = []
    for locid in gendict:
        locus = gendict[locid]
        possize = int(locus['Stop']) - int(locus['Start'])
        spacersize = sum(len(s) for s in locus['Spacers'].values())
        repeatsize = len(locus['RepeatSeq']) * (len(locus['Spacers']) + 1)
        sizediff = possize - spacersize - repeatsize
        offsetlist.append((sizediff, locid))
    return offsetlist

def delkeys(dictionary, keylist):
    for key in keylist:
        if key in dictionary:
            del dictionary[key]
    return dictionary



def prunedict(gendict):
    #remove loci with incomplete fields
    incomplist = findincompleterecords(gendict)
    gendict = delkeys(gendict, incomplist)
    #remove loci with offsets greater than 0
    offsetlist = findsizeoffsets(gendict)
    locitodel = [x[1] for x in offsetlist if x[0] > 0]
    gendict = delkeys(gendict, locitodel)
    return gendict
if __name__ == '__main__':
    with open('gendict.pickle', 'rb') as f:
        gendict = pickle.load(f)
    # counter = collections.Counter(zip(*offsetlist)[0])
    gendict = prunedict(gendict)
    with open('gendictpruned.pickle', 'wb') as f:
        pickle.dump(gendict, f, protocol=pickle.HIGHEST_PROTOCOL)
