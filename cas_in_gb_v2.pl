#!/usr/bin/env perl

use strict;
use warnings;

#Use with 
# ./cas_in_gb_v2.pl file_name.txt

my $gb_file = shift || "data/Genbank_example.txt";

open my $fh, "<", $gb_file;
my $file_contents = join( "", <$fh> );
close $fh;

my %cas = ( Cas3 => 'Type I',
	    Cas9 => 'Type II', 
	    Cas10 => 'Type III');

my $cas_str = join("|",keys %cas );
my $cas_regex = qr/($cas_str)/;

my @loci = split(/LOCUS\s+/s, $file_contents );

for my $l (@loci ) {
  next if !$l;
  my ($ac_code) = ($l =~ /^(\w+)/); 
  my ( @genes ) = split(/\n\s+gene\s{2,}/, $l );
  for my $g ( @genes ) {
    my ($seq ) = ($g =~ /(\d+\.\.\d+)/);
    next if $g !~ m{/product};
    my ($product) = ($g =~ m{/product="([^\"]+)"}g);
    if ( $product =~ $cas_regex ) {
      my $pattern = $1;
      print "->$ac_code: $seq ==> $cas{$pattern}\n";
    }
  }
}
