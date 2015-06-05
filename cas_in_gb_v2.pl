#!/usr/bin/env perl

use strict;
use warnings;

use v5.14;

use File::Slurp::Tiny 'read_file';

my $gb_file = shift || "data/Genbank_example.txt";

my $file_contents = read_file( $gb_file);


my @loci = split(/LOCUS\s+/s, $file_contents );

for my $l (@loci ) {
  next if ! ($l =~ m{/product});
  my ( @genes ) = split(/\n\s+gene\s{2,}/, $l );
  for my $g ( @genes ) {
    next if $g !~ m{/product};
    my ($seq ) = ($g =~ /(\d+\.\.\d+)/);
    my (@products) = ($l =~ m{/product="([^\"]+)"}g);
    say "-> $seq ==> ", join(" - ", @products );
  }
}
