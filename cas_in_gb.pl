#!/usr/bin/env perl

my $gbfile = shift @ARGV;

open GB_FH, "<", $gbfile;

my $line = readline GB_FH;

my $in_organism = "";
my $curr_taxonomy = "";
my $curr_location = "";

while (my $line = readline GB_FH) {
	if ($line =~ m/^\s+(ORGANISM\s+.+)$/) {
		print "$1\n";
		$curr_taxonomy = "TAXONOMY ";
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
