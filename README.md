phageParser
===========

##Taming Pathogens

Pathogens have played a crucial role in human history. The most iconic is the  Black Death, which was (is thought to have been) caused by bacteria named Yersinia pestis, in the 14th century. A more recent example is AIDS caused by human immunodeficiency virus (HIV). Hence, understanding fundamentals of the host-pathogen interactions have been a central problem in epidemiology, which with the genomic age requires tools from quantitative fields to make sense of the large amount of sequencing data that is being generated.  

A virus-bacteria interaction provides the simplest host-pathogen pair and is crucial for functioning of human gut to marine ecosystems. Interestingly, amechanism for “adaptive” immune system in bacteria against the viruses infecting them (usually referred to as phages) was discovered just a few years ago. Curiously, much like the anti-virus software and intrusion detection systems that rely on detecting patters found in malicious code, bacteria keeps a dynamic library of small pieces of phage genomes (spacers) to detect and neutralize phage attacks.  

The basic problem of understanding how this immune system works is to understand the pattern of spacers on phage genomes: how many per phage genome, where on a phage genome, if the spacers containing regions are more or less dynamic compared to the rest of the phage genome etc. Since we have a large number of sequenced phages and a library of spacers from a variety of bacteria - ranging from deadly human pathogens such as tuberculosis to bacteria that live in our guts - we can attempt to aggregate this information to develop a more “complete” understanding of phage-bacteria interactions.

##Data Challenge

Happily, much of the existing data needed to understand bacteria / phage interaction has been released openly to the public and is available over the web; the current challenge is to help extract the relevant parts of that huge database, and automate the production of targeted datasets for these studies.  More details are in the issue tracker!
