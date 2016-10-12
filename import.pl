#!/usr/bin/env perl

use strict;
use warnings;

use Graph::Directed;

my $infile = $ARGV[0];

my $g = Graph::Directed->new;

open(FH, "<", $infile) || die "Unable to open file '$infile'\n";

while (<FH>)
{
    next unless (/^>/);
}

close(FH) || die "Unable to close file '$infile'\n";
