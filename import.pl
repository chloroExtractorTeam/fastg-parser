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

    # extract Node information
    my @nodes = $_ =~ /(EDGE_[^:;,]+)/g;

    foreach my $node (@nodes)
    {

	my $nodename = $node;
	# estimate node basename

	my $nodebasename = $nodename;
	$nodename=~s/'$//; # delete a ' if necessary
	my $rev_nodename = $nodebasename."'";

	# check if the node is known
	unless ($g->has_vertex($nodebasename))
	{
	    $g->add_edge($nodebasename, $rev_nodename);
	    $g->add_edge($rev_nodename, $nodebasename);	    
	}
    }
}

print "$g";

close(FH) || die "Unable to close file '$infile'\n";
