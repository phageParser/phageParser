#!/usr/bin/perl

use:strict;
use:warnings;

#read in file
$file = $ARGV[0];

open(FILE, "<$file") or die "cannot open infile!\n";

#line by line, store data for each organism in hashes
my %organismhash; #storing info related to each unique organism
my $organism; #storing organism name
my %proteins; #hash of every unique protein product

while(my $line = <FILE>) {
    chomp $line;
    if ($line =~ m/ORGANISM/) { #if the line has the word "ORGANISM" in it
	my ($element1, $element2, @rest) = split("  ", $line);
	$organism=$element2;
	#print STDERR "$organism is the current organism\n";
    }
    elsif ($line =~ m/TAXONOMY/) { #if it's the taxonomy line, we want the first element (the kingdom)
	$line =~ m/TAXONOMY (.[A-Za-z]{1,});.*/;
	my $lineage = $1;
	$organismhash{$organism}{"kingdom"}=$lineage;
	#print STDERR "the organism: $organism belongs to the kingdom $lineage\n";
    }
    elsif($line =~ m/PRODUCT/i) { #if it's a product line, we just want the protein product name
	my @products = (); #array to store protein products
	my @productline = split('"', $line);
	#print STDERR "$productline[1]\n\n";
	$proteins{$productline[1]}=0;
	#push(@products, $productline[1]); #this line and the line below are the real 'proper' non hacky way to populate an array and store it in a hash, but i can't get the output in the foreach loop below to behave how i expect, so the line two down is another way that works
	#$organismhash{$organism}{"products"}=\@products;
        push(@{$organismhash{$organism}{"products"}}, $productline[1]);
	#print STDERR "the organism: $organism has a CAS protein product: $productline[1]\n";
    }
}


#now iterate over organisms again, printing (to STDOUT) info we collected for each
#print header for outfile
print STDOUT "species\tkingdom\t"; #print species and kingdom
#iterate over unique proteins in the protein hash, join by tab, print to the header
foreach my $proteinproduct (sort keys %proteins) {
    print STDOUT "$proteinproduct\t";
}
print STDOUT "\n";

foreach my $bacteria (sort keys %organismhash) {
    my $type=$organismhash{$bacteria}{"kingdom"};
    #print the species and the kingdom
    print STDOUT "$bacteria\t$type\t";
    foreach my $a (@{$organismhash{$bacteria}{"products"}}) { #iterate over the protein products for each organism to count them
	$proteins{$a}++;
	#print STDOUT "$a\tline58\n";
    }
    foreach my $proteinproducts (sort keys %proteins) { #now print out the numbers we've counted for each protein product
	print STDOUT "$proteins{$proteinproducts}\t";
    }
    print STDOUT "\n";
    foreach my $proteinproducts (keys %proteins) { #clear out the hash of any incrementation
	$proteins{$proteinproducts}=0;
    }
}
