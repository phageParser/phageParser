phageParser
===========

##Taming Pathogens

Pathogens have played a crucial role in human history. The most iconic is the Black Death, which was (is thought to have been) caused by bacteria named Yersinia pestis in the 14th century. A more recent example is AIDS caused by human immunodeficiency virus (HIV). Hence, understanding fundamentals of the host-pathogen interactions has been a central problem in epidemiology, which with the genomic age requires tools from quantitative fields to make sense of the large amount of sequencing data that is being generated.  

A virus-bacteria interaction provides the simplest host-pathogen pair and is crucial for functioning of human gut to marine ecosystems. Interestingly, a mechanism for “adaptive” immune system in bacteria against the viruses infecting them (usually referred to as phages) was discovered just a few years ago, called CRISPR-Cas. Curiously, much like the anti-virus software and intrusion detection systems that rely on detecting patterns found in malicious code, bacteria keeps a dynamic library of small pieces of phage genomes (spacers) to detect and neutralize phage attacks.  

The basic problem of understanding how this immune system works is to understand the pattern of spacers on phage genomes: how many per phage genome, where on a phage genome, if the spacers containing regions are more or less dynamic compared to the rest of the phage genome, etc. Since we have a large number of sequenced phages and a library of spacers from a variety of bacteria - ranging from deadly human pathogens such as tuberculosis to bacteria that live in our guts - we can attempt to aggregate this information to develop a more “complete” understanding of phage-bacteria interactions.

##Data Challenge

Happily, much of the existing data needed to understand bacteria / phage interaction has been released openly to the public and is available over the web; the current challenge is to help extract the relevant parts of that huge database, and automate the production of targeted datasets for these studies.  More details are in the issue tracker!

##Relevant Literature

[CRISPR-Cas Systems: Prokaryotes Upgrade to Adaptive Immunity](http://www.cell.com/molecular-cell/abstract/S1097-2765%2814%2900216-0): a very good review paper on the CRISPR-cas system, the biological backdrop of this project.

##Installation

This package depends on Biopython:

```
sudo pip install Biopython
```

Also, make yourself a directory `phageParser/output` - some data cleaning scripts will dump their results there.

##Usage

 - To get a phage dataset, take a fasta-formatted list of genes (example in `data/velvet-distinct-spacers.fasta`) and upload to http://phagesdb.org/blast/ - example result in `data/blast-phagesdb.txt`

 - To clean up the results returned from phagesdb.org, change the raw filename in `filterByExpect.py` from `data/blast-phagesdb.txt` to whatever file contains the results from the BLAST search, then do

 `python filterByExpect.py`

 The result will be written to a file in `phageParser/output`, in a CSV formatted as
```
 Query, Name, Length, Score, Expect, QueryStart, QueryEnd, SubjectStart, SubjectEnd
```

 with one header row (see [#1](https://github.com/goyalsid/phageParser/issues/1) for discussion and details)

 - To query NCBI for full genomes, do
```
 cat accessionNumber.txt | python acc2gb.py youremail@yourinstitution.org > NCBIresults.txt
```
 where `accessionNumber.txt` contains a list of accession numbers of interest; results will be dumped to `NCBIresults.txt` - see [#2](https://github.com/goyalsid/phageParser/issues/2) for ongoing development here. 

##Alternative usage for code sprint materials

All of the following assumes you are using the reference CRISPR database set of spacers (file `spacerdatabase.txt`). Things *should* work with other spacer files; however there are several things hard-coded that might break. `filterByExpect.py` assumes the header line for each spacer is a number, for example, and `bac_name` is hardcoded in `interactions.py` as the 8th to 16th characters of the file name.

- To get individual spacer files for each bacteria species in the reference set, run `CRISPR_db_parser` on with the input file `spacerdatabase.txt` (downloaded from the [Utilities](http://crispr.u-psud.fr/crispr/BLAST/Spacer/Spacerdatabase) page of CRISPRdb). Put all the output files in a folder `/spacers` under `data`.

- Make folders `data/phages` and `/output`. The current files in `data/spacers` and `data/phages` are examples.

- Blasting of spacer-containing files against the phage database can be done locally (handy if you have many files to blast). Download a local version of blast (blast+) [here](http://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download) and find/follow instructions for your OS. (We used [these instructions](http://www.ncbi.nlm.nih.gov/books/NBK52637/) for Windows successfully.) Put the file `Mycobacteriophages-All.fasta` (in [data](https://github.com/goyalsid/phageParser/tree/master/data) folder) into the main blast+ directory and use the command `makeblastdb -in "Mycobacteriophages-All.fasta" -dbtype nucl -title PhageDatabase -out phagedb` to create a blast-ready database. It's possible to combine multiple databases into one blastable database by including more than one filename between the quotes in the `-in` command (i.e. the ENA phage database or NCBI virus database). Now you should be able to run the script `BLAST_loop.py`, but make sure directory names are correct - probably `BLAST_loop.py` will need to be run from inside wherever you installed blast+. 

- run `filterByExpectPhages.py`, which essentially runs filterByExpect.py on all files in the `/phages` folder. These will be saved to `/output`.

- make a directory called `sorted` under `output`. run `orderByExpect.py`, which rearranges the results of filterByExpectPhages in each file to be in order of lowest to highest expect value.

- run `interactions.py`, which makes a json file `json.txt` for visualization in cytoscape.js. 

Visualization

- paste the contents of `json.txt` into the `elements[]` field in the file `ui.js`. This creates the structure needed for cytoscape.js to plot stuff. Various style fields can be changed, see [cytoscape.js](http://js.cytoscape.org/) for documentation (or ask @MaxKFranz for help).

- paste the file `index.html` into a web browser. 
