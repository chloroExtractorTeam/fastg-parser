#!/ussr/bin/env perl

use strict;
use warnings;

use Graph;

my $infile = $ARGV[0];

open(FH, "<", $infile) || die "Unable to open file '$infile'\n";

while (<FH>)
{
    next unless (/^>/);

    
}

close(FH) || die "Unable to close file '$infile'\n";
