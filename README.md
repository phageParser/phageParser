phageParser
===========

phageParser is a project to extract and organize CRISPR information from open genetic data.

## What is this tool?

Many bacterial and archaeal genomes have been sequenced, and [a large fraction of them have CRISPR systems](http://crispr.u-psud.fr/crispr/), ranging from deadly human pathogens to archaea living in the harshest environments on earth. Some CRISPR systems have been studied very well, and more is being discovered about CRISPR every day. phageParser is a tool to collect this growing pool of information and generate versatile and useful annotations. These are some of the annotations we include:
* Spacer matches to known phages and prophages
* Phage genome content near spacer matches
* Spacer self-matches to host genome
* *cas* gene content and inferred CRISPR type

These annotations will be collected in a database that can be queried through a GUI. Neither of these exist yet, and we are actively looking for contributors. 

This tool is currently in development, and it will always be possible to modify and enhance what is included as CRISPR research moves forward. We welcome suggestions for features or annotations you'd like to see! To suggest a feature, create an issue in our [issue tracker](https://github.com/goyalsid/phageParser/issues).

## Who is this for?

phageParser is for anyone interested in exploring what we know about CRISPR systems in nature. This includes researchers, educators, and the general public.

## Where can I get involved?

We need many different skills and areas of expertise to build this tool, and you can help! 
* [Good first bugs](https://github.com/goyalsid/phageParser/labels/good%20first%20bug) include documentation and coding tasks that are doable by a newcomer. Mentoring is available for these tasks. 
* Do you know about CRISPR biology? Issues labeled [science](https://github.com/goyalsid/phageParser/issues?q=is%3Aissue+is%3Aopen+label%3Ascience) are things we need people with science background to work on.
* Are you interested in contributing to project documentation? Any issues labeled [documentation](https://github.com/goyalsid/phageParser/issues?q=is%3Aissue+is%3Aopen+label%3Adocumentation) are ways to create or improve our docs. 
* Do you know about databases? We're just starting to think about how to structure our data - join the discussion in issue #64.
* Do you know about Python and/or developing code? Check out our [code](https://github.com/goyalsid/phageParser/labels/code)-specific issues.

## About the CRISPR system

**C**lustered **R**egularly **I**nter-spaced **S**hort **P**alindromic **R**epeats (CRISPR), and associated proteins (Cas) are part of the CRISPR-Cas system in bacteria. First observed in 1987 (Ishino et al., 1987), the CRISPR system acts like the immune system for bacteria.

In humans, when a virus enters the body, specialized immune cells are often quick to recognize the foreign pathogen and kill it. Bacteria do not have the benefit of a legion of cells to protect them against viruses, so they use the CRISPR-Cas mechanism that begins with the creation of *spacer* sequences from the invading virus' genetic material followed by the production of small interfering crRNAs. Finally, when the bacterium is invaded again, it uses the crRNAs to target the viral nucleic acid and prevent infection. 

The spacers acquired are always flanked by bacterial DNA forming a library of small pieces of phage genomes from which to form crRNAs to target future viral genetic material (see Figure below). 
![alt text](https://github.com/goyalsid/phageParser/blob/master/CRISPR_locus_diagram.png?raw=true "CRISPR Locus")

Amazingly, even though this process is not entirely common, when it does happen, the bacteria are essentially vaccinated against future attacks and can pass on their phage genome libraries to future generations.

However, more research is still needed to better understand the mechanisms that make this process possible. 

*Ishino, Y., Shinagawa, H., Makino, K., Amemura, M., and Nakata, A. (1987). Nucleotide sequence of the iap gene, responsible for alkaline phosphatase isozyme conversion in Escherichia coli, and identification of the gene product. J. Bacteriol. 169, 5429–5433.


##Relevant Literature

[CRISPR-Cas Systems: Prokaryotes Upgrade to Adaptive Immunity](http://www.cell.com/molecular-cell/abstract/S1097-2765%2814%2900216-0): a very good review paper on the CRISPR-cas system, the biological backdrop of this project.

##Installation

This package depends on Biopython, which can be installed as follows:

```
sudo make install-deps
```

##Usage

There are several usage options depending on what data outcome is desired.

###Usage - Blast individual bacteria files and get phage info from NCBI

 - To get a phage dataset, take a fasta-formatted list of genes (example in [`data/velvet-distinct-spacers.fasta`](data/velvet-distinct-spacers.fasta)) and upload to http://phagesdb.org/blast/ - example result in [`data/blast-phagesdb.txt`](data/blast-phagesdb.txt)

 - To clean up the results returned from phagesdb.org, you can call the Make target filter_by_expect, as in the example below.

 `make filter_by_expect infile=data/blast-phagesdb.txt output=output/ threshold=0.21`

 The result will be written to a file in `output/`, in a CSV formatted as
```
 Query, Name, Length, Score, Expect, QueryStart, QueryEnd, SubjectStart, SubjectEnd
```

 with one header row (see [#1](https://github.com/goyalsid/phageParser/issues/1) for discussion and details)

 - To query NCBI for full genomes, do
```
 cat accessionNumber.txt | python acc2gb.py youremail@yourinstitution.org > NCBIresults.txt
```
 where `accessionNumber.txt` contains a list of accession numbers of interest; results will be dumped to `NCBIresults.txt` - see [#2](https://github.com/goyalsid/phageParser/issues/2) for ongoing development here. 

###Usage - Blasting multiple bacteria files and visualizing interactions

All of the following assumes you are using the reference CRISPR database set of spacers (file `spacerdatabase.txt`). Things *should* work with other spacer files; however there are several things hard-coded that might break. [`filterByExpect.py`](filterByExpect.py) assumes the header line for each spacer is a number, for example, and `bac_name` is hardcoded in [`interactions.py`](interactions.py) as the 8th to 16th characters of the file name.

- To get individual spacer files for each bacteria species in the reference set, run `CRISPR_db_parser` on with the input file `spacerdatabase.txt` (downloaded from the [Utilities](http://crispr.u-psud.fr/crispr/BLAST/Spacer/Spacerdatabase) page of CRISPRdb). The output files will be saved in the folder `data/spacers`.

- Make folders `data/phages` and `/output`. The current files in `data/spacers` and `data/phages` are examples.

- Blasting of spacer-containing files against the phage database can be done locally (handy if you have many files to blast). Download a local version of blast (blast+) [here](http://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download) and find/follow instructions for your OS. (We used [these instructions](http://www.ncbi.nlm.nih.gov/books/NBK52637/) for Windows successfully.) Put the file `Mycobacteriophages-All.fasta` (in [data](https://github.com/goyalsid/phageParser/tree/master/data) folder) into the main blast+ directory and use the command `makeblastdb -in "Mycobacteriophages-All.fasta" -dbtype nucl -title PhageDatabase -out phagedb` to create a blast-ready database. It's possible to combine multiple databases into one blastable database by including more than one filename between the quotes in the `-in` command (i.e. the ENA phage database or NCBI virus database). Now you should be able to run the script `BLAST_loop.py`, but make sure directory names are correct - probably `BLAST_loop.py` will need to be run from inside wherever you installed blast+. 

- run [`filterByExpectPhages.py`](filterByExpectPhages.py), which essentially runs filterByExpect.py on all files in the `/phages` folder. These will be saved to `/output`.

- make a directory called `sorted` under `output`. run [`orderByExpect.py`](orderByExpect.py), which rearranges the results of filterByExpectPhages in each file to be in order of lowest to highest expect value.

- run [`interactions.py`](interactions.py), which makes a json file `json.txt` for visualization in cytoscape.js. 

Visualization

####Install Front-End Dependencies (to visualize in browser). 
- Install [node.js](https://nodejs.org/). Installing node will also install the node package manager (NPM).

- Install [bower](bower.io) 

- paste the contents of `json.txt` into the `elements[]` field in the file `ui.js`. This creates the structure needed for cytoscape.js to plot stuff. Various style fields can be changed, see [cytoscape.js](http://js.cytoscape.org/) for documentation (or ask @MaxKFranz for help).

- paste the file `index.html` into a web browser. 

###Usage - Detecting CRISPR type from bacterial genome metadata

- Start with a list of bacteria of interest - in this case, it's all the bacteria from CRISPRdb that had hits to a conglomerate of phage databases - `bac_accession_list.txt`. 

- Next is to fetch bacterial genome data from NCBI. Run the following:
```
 cat bac_accession_list.txt | python acc2gb.py youremail@yourinstitution.org > NCBIresults.txt
```
Be warned that this will take a long time (~1-2 hours) because the list is long. For testing, shorten the list to only a few accession numbers. 

- Run the script `trimGenbankDNA.py` to get rid of unnecessary data and make the file size more manageable. 

- Run [`cas_in_gb.pl`](cas_in_gb.pl) (it's in Perl) to detect which Cas genes are in each organism.
