phageParser
===========

##Taming Pathogens

Pathogens have played a crucial role in human history. The most iconic is the  Black Death, which was (is thought to have been) caused by bacteria named Yersinia pestis, in the 14th century. A more recent example is AIDS caused by human immunodeficiency virus (HIV). Hence, understanding fundamentals of the host-pathogen interactions have been a central problem in epidemiology, which with the genomic age requires tools from quantitative fields to make sense of the large amount of sequencing data that is being generated.  

A virus-bacteria interaction provides the simplest host-pathogen pair and is crucial for functioning of human gut to marine ecosystems. Interestingly, amechanism for “adaptive” immune system in bacteria against the viruses infecting them (usually referred to as phages) was discovered just a few years ago. Curiously, much like the anti-virus software and intrusion detection systems that rely on detecting patters found in malicious code, bacteria keeps a dynamic library of small pieces of phage genomes (spacers) to detect and neutralize phage attacks.  

The basic problem of understanding how this immune system works is to understand the pattern of spacers on phage genomes: how many per phage genome, where on a phage genome, if the spacers containing regions are more or less dynamic compared to the rest of the phage genome etc. Since we have a large number of sequenced phages and a library of spacers from a variety of bacteria - ranging from deadly human pathogens such as tuberculosis to bacteria that live in our guts - we can attempt to aggregate this information to develop a more “complete” understanding of phage-bacteria interactions.

##Data Challenge

Happily, much of the existing data needed to understand bacteria / phage interaction has been released openly to the public and is available over the web; the current challenge is to help extract the relevant parts of that huge database, and automate the production of targeted datasets for these studies.  More details are in the issue tracker!

##Installation

This package depends on Biopython:

```
sudo pip install Biopython
```

Also, make yourself a directory `phageParser/output` - some data cleaning scripts will dump their results there.

##Usage

 - To get a phage dataset, take a fasta-formatted list of genes (example in `data/velvet-distinct-spacers.fasta`) and upload to http://phagesdb.org/blast/ - example result in `data/blast-phagesdb.txt`

 - To clean up the results returned from phagesdb.org, change the raw filename in `getAccession.py` from `data/PhagesDB_Data.txt` to whatever file contains the results from the BLAST search, then do

 `python filterByExpect.py`

 The result will be written to a file in `phageParser/output`, in a CSV formatted as
```
 Query, Name, Length, Score, Expect, QueryStart, QueryEnd, SubjectStart, SubjectEnd
```

 with one header row (see #1 for discussion and details)

 - To query NCBI for full genomes, do
```
 cat accessionNumber.txt | python acc2gb.py youremail@yourinstitution.org > NCBIresults.txt
```
 where `accessionNumber.txt` contains a list of accession numbers of interest; results will be dumped to `NCBIresults.txt`
