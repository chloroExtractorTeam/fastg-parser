#!/usr/bin/env perl

use strict;
use warnings;

use Graph;
#use Graph::Directed;
#use Graph::Subgraph;

my $infile = $ARGV[0];

my $g = Graph->new(directed => 1);

open(FH, "<", $infile) || die "Unable to open file '$infile'\n";

while (<FH>)
{
    next unless (/^>/);

    my $line = $_;

    chomp($line);

    # delete a ; if necessary
    $line =~ s/;$//;
    # delete the leading >
    $line =~ s/^>//;

    # extract Node information
    my @nodes = $line =~ /((?:EDGE|NODE)_[^:;,]+)/g;

    foreach my $node (@nodes)
    {

	my $nodename = $node;
	# estimate node basename

	my $nodebasename = $nodename;
	$nodebasename =~ s/\'//g;
	my $rev_nodename = $nodebasename."'";

	# check if the node is known
	unless ($g->has_vertex($nodebasename))
	{
	    $g->add_edge($nodebasename, $rev_nodename);
	    $g->add_edge($rev_nodename, $nodebasename);
	    print STDERR "Created Node: $node ($nodebasename)\n";
	}
    }

    # do we have a graph?
    if ($line =~ /:/)
    {
	@nodes = split(/:/, $line);
	if (@nodes != 2)
	{
	    die "Something went wrong";
	}
	my ($from, $to) = @nodes;

	# check if $to contains more than one node (delimeter would be an ,)
	my @to_nodes = ($to);
	if ($to =~ /,/)
	{
	    @to_nodes = split(/,/, $to);
	}

	foreach my $to_node (@to_nodes)
	{
	    $g->add_edge($from, $to_node);
	}
    }
}

foreach my $v ($g->vertices)
{
    #print STDERR "Vertices: $v\n";
}

foreach my $wcc ($g->weakly_connected_components())
{
    #print STDERR Dumper($wcc); use Data::Dumper;

    next if (@{$wcc} == 2);

    my $c = subgraph($g, @{$wcc});

    next unless ($c->is_cyclic && $c->vertices > 2);

#    print Dumper($wcc); use Data::Dumper;

    print "$c", "\n";
}

#print $g;

close(FH) || die "Unable to close file '$infile'\n";

# from http://stackoverflow.com/questions/8729302/creating-a-subgraph-of-a-graph-induced-by-node-list

# subgraph ($graph, @node_list);
# return subgraph (with the same setup)
# induced by node list
sub subgraph {
    my $self = shift;
    my $new = $self->new;
    my @edges;
    foreach my $v(@_) {
        $self->has_vertex($v) or next;
        $new->add_vertex($v);
        foreach my $u(@_) {
            $self->has_edge($u, $v) and push @edges, $u, $v;
        };
    };
    $new->add_edges(@edges);
    return $new;
};
