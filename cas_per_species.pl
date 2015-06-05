#!/usr/bin/perl

use:strict;
use:warnings;

#read in file
$file = $ARGV[0];

open(FILE, "<$file") or die "cannot open infile!\n";

#line by line, store data for each organism in hashes
my %organismhash;
my $organism;
my %proteins; #hash of every unique protein product
my $kingdom=1; #placeholder array, used later in both loops

while(my $line = <FILE>) {
    chomp $line;
    if ($line =~ m/ORGANISM/) { #if the line has the word "ORGANISM" in it
	my ($element1, $element2, @rest) = split("  ", $line);
	$organism=$element2;
	print STDERR "$organism is the current organism\n";
    }
    elsif ($line =~ m/TAXONOMY/) { #if it's the taxonomy line, we want the first element (the kingdom)
	$line =~ m/TAXONOMY (.[A-Za-z]{1,});.*/;
	my $lineage = $1;
	$organismhash{$organism}{"kingdom"}=$lineage;
	print STDERR "the organism: $organism belongs to the kingdom $lineage\n";
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


#iterate over organisms, collecting information for each organism
#print header for outfile
print 
#print STDERR "species\tkingdom\t"; #print species and kingdom
#iterate over unique proteins in the protein hash, join by tab, print to the header

foreach my $proteinproducts (sort keys %proteins) {
    #print STDERR "$proteinproducts\t";
}
#print STDERR "\n";

foreach my $bacteria (sort keys %organismhash) {
    print "$bacteria\n";
    my $type=$organismhash{$bacteria}{$kingdom};
    foreach my $a (@{$organismhash{$bacteria}{"products"}}) {
	$proteins{$a}++;
	#print STDERR "$a\tline58\n";
    }
    foreach my $proteinproducts (sort keys %proteins) {
	#print STDERR "$proteins{$proteinproducts}\t";
    }
    #print STDERR "\n";
    foreach my $proteinproducts (keys %proteins) { #clear out the hash of any incrementation
	$proteins{$proteinproducts}=0;
    }
}

#print collected info to file
