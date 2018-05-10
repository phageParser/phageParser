phageParser
===========

phageParser is a project to extract and organize CRISPR information from open genetic data.

Come chat with us!    
[![Gitter chat](https://badges.gitter.im/gitterHQ/gitter.png)](https://gitter.im/phageParser/Lobby)

## What is this tool?

Many bacterial and archaeal genomes have been sequenced, and [a large fraction of them have CRISPR systems](http://crispr.i2bc.paris-saclay.fr/), ranging from deadly human pathogens to archaea living in the harshest environments on earth. Some CRISPR systems have been studied very well, and more is being discovered about CRISPR every day. phageParser is a tool to collect this growing pool of information and generate versatile and useful annotations. These are some of the annotations we include:
* Spacer matches to known phages and prophages
* Phage genome content near spacer matches
* Spacer self-matches to host genome
* *cas* gene content and inferred CRISPR type

We will collect these annotations in a database that can users can query through a GUI (graphical user interface). Neither of these exist yet, and we are looking for contributors!

This tool is currently in development, and it will always be possible to modify and enhance what is included as CRISPR research moves forward. We welcome suggestions for features or annotations you'd like to see! To suggest a feature, create an issue in our [issue tracker](https://github.com/goyalsid/phageParser/issues).

## Who is this for?

phageParser is for anyone interested in exploring what we know about CRISPR systems in nature. This includes researchers, educators, and the general public.

## Where can I get involved?

We need many different skills and areas of expertise to build this tool, and you can help!
* Check out the [open canvas](https://github.com/phageParser/phageParser/blob/master/documentation/Open_canvas.png) - this is a short outline of the project goals and plans.
* Check out the [Roadmap](https://github.com/goyalsid/phageParser/issues/112) for an overview of where we're going and when.
* [Good first bugs](https://github.com/goyalsid/phageParser/labels/good%20first%20issue) include documentation and coding tasks that are doable by a newcomer. Mentoring is available for these tasks.
* Do you know about CRISPR biology? Issues labeled [science](https://github.com/goyalsid/phageParser/issues?q=is%3Aissue+is%3Aopen+label%3Ascience) are things we need people with science background to work on.
* Are you interested in contributing to project documentation? Any issues labeled [documentation](https://github.com/goyalsid/phageParser/issues?q=is%3Aissue+is%3Aopen+label%3Adocumentation) are ways to create or improve our docs.
* Do you know about databases? We're just starting to think about how to structure our data - join the discussion in issue #64.
* Do you know about Python and/or developing code? Check out our [code](https://github.com/goyalsid/phageParser/labels/code)-specific issues.
* Check out our [code of conduct](https://github.com/goyalsid/phageParser/blob/master/CODE-OF-CONDUCT.md) which applies to all maintainers and contributors to this project.

## About the CRISPR system

**C**lustered **R**egularly **I**nterspaced **S**hort **P**alindromic **R**epeats (CRISPR), and associated proteins (Cas) are part of the CRISPR-Cas system in bacteria. First observed in 1987 (Ishino et al., 1987), the CRISPR system is an adaptive immune system for bacteria.

When a virus enters a human body, specialized immune cells are often quick to recognize the virus invader and kill it. Bacteria do not have the benefit of millions of immune cells to protect them against viruses, but they have something else: CRISPR-Cas. The CRISPR-Cas immune response begins with the creation of *spacer* sequences from the invading virus' DNA followed by the production of small interfering crRNAs. Finally, when the bacterium is invaded again, the crRNAs recognize and cut the viral DNA, preventing infection.

Bacteria store their acquired spacers in their own DNA. The spacers are flanked by short pieces of bacterial DNA called *repeats* (see figure below).

![CRISPR Locus Diagram](/documentation/CRISPR_locus_diagram.png?raw=true "CRISPR Locus")

Amazingly, CRISPR-Cas immunity is both *adaptive* and *hereditary*! After acquiring a spacer, bacteria are both protected against future virus attacks and they can pass on their spacer libraries to their descendants.

More research is needed to better understand how bacteria use their CRISPR systems in nature.

*Ishino, Y., Shinagawa, H., Makino, K., Amemura, M., and Nakata, A. (1987). Nucleotide sequence of the iap gene, responsible for alkaline phosphatase isozyme conversion in Escherichia coli, and identification of the gene product. J. Bacteriol. 169, 5429–5433.


## Relevant Literature

[CRISPR-Cas Systems: Prokaryotes Upgrade to Adaptive Immunity](http://www.cell.com/molecular-cell/abstract/S1097-2765%2814%2900216-0): a very good review paper on the CRISPR-cas system, the biological backdrop of this project.

## Installation

You can download the source code of the project by git:
`git clone https://github.com/phageParser/phageParser.git`

After getting the local copy of the project, it is generally a good idea
to create an isolated environment that belongs to the project and its specific
packages. For this, python has a tool called virtualenv that can help create a
python instance that has different packages than the system's version. To get started:

Make sure you have python3 in your system, if not,
you can download python3 via their [website](https://www.python.org/downloads/)

You can then install virtualenv package by pip

`pip install virtualenv`

For creating a virtualenv with a specific python version, you can supply the
path of the python binary as an argument. The virtual python instances are conventionally
kept in one place, usually in `~/.virtualenvs`. You can create the folder and make
an environment for phageParser as such:

`mkdir ~/.virtualenvs && cd "$_"`

`python3 -m venv ~/.virtualenvs/pparserdev`

You now have a separate environment which you can use to contribute
phageParser. Whenever you're developing for phageParser, use the following command
to activate the environment:

`source ~/.virtualenvs/pparserdev/bin/activate`

To install the required libraries for phageParser, after heading to the project folder
containing [`requirements.txt`](requirements.txt), activate the project environment and run the following command:

`pip install -r requirements.txt`

For viewing the database, we recommend the [Firefox SQLite Manager plugin](https://addons.mozilla.org/en-US/firefox/addon/sqlite-manager/). Once installed, launch it from the 'Tools' menu in Firefox. 
