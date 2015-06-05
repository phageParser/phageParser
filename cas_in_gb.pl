#!/usr/bin/env perl

use strict;
use warnings;

my $gbfile = shift || die "Usage: $0 <gb_file> \b";

open my $gb_fh, "<", $gbfile;

my $line = readline $gb_fh;

my $accession = "";
my $curr_taxonomy = "";
my $curr_location = "";

while (my $line = readline $gb_fh) {
	if ($line =~ m/^\s*(ORGANISM\s+.+)$/) {
		print "$1\n$accession\n";
		$curr_taxonomy = "TAXONOMY ";
	} elsif ($line =~ m/^\s*(ACCESSION\s+.+)$/) {
		$accession = "$1";
	} elsif (($line !~ m/^ {12}\S/) && ($curr_taxonomy ne "")) {
		print "$curr_taxonomy\n";
		$curr_taxonomy = "";
	} elsif ($curr_taxonomy ne "") {
		# we are still in the taxonomy line; concat.
		$line =~ s/\s//g;
		$line =~ s/\n//;
		$curr_taxonomy .= $line;
	} elsif ($line =~ m/^\s+\/(product=.*Cas.*)$/) {
		print "$1\t$curr_location\n";
		$curr_location = "";
	} elsif ($line =~ m/^\s+CDS\s+(.*)$/) {
		$curr_location = $1;
	}
}
