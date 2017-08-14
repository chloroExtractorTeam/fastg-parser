#!/usr/bin/env perl

use strict;
use warnings;

use Graph;
use Term::ProgressBar;

my $infile = $ARGV[0];

my $g = Graph->new(directed => 1);

my %names = ();
my @seq = ();
my @digraphs = ();
my %seq2seqname = ();
my $lastname = "";

my $MINNODES = 3;
my $MAXNODES = 100;
my $MINSEQLEN = 25000;
my $MAXSEQLEN = 1000000;


open(FH, "<", $infile) || die "Unable to open file '$infile'\n";

while (<FH>)
{
    my $line = $_;

    next unless ($line =~ /^>/ || $lastname);

    chomp($line);

    if ($line =~ /^>/)
    {
	# delete a ; if necessary
	$line =~ s/;$//;
	# delete the leading >
	$line =~ s/^>//;

	# check if sequence contains a digraph
	my ($seqname, @digraph) = split(/:/, $line);

	if (@digraph)
	{
	    @digraph = split(/,/, join(",", @digraph));
	}

	# does the seqname contain a ' to indicate the reverse complement?
	my $reverse = 0; # assume non reverse sequence
	if ($seqname =~ /'/)
	{
	    $reverse = 1;
	    $seqname =~ s/'//g;
	}

	# if digraph information are given, store them
	push(@digraphs, { from => $seqname, from_rev => $reverse, to => \@digraph }) if (@digraph);

	# sequence already known?
	unless (exists $names{$seqname})
	{
	    $lastname = $seqname;

	    push(@seq, "");
	    $names{$seqname} = @seq-1;
	    $seq2seqname{@seq-1} = $seqname;
	} else {
	    # in case the sequence is already known
	    $lastname = undef;
	}

    } else {
	$seq[$names{$lastname}] .= $line;
    }
}

close(FH) || die "Unable to close file '$infile'\n";

# build graph information
foreach my $digraph (@digraphs)
{
    # create the starting node
    unless (exists $names{$digraph->{from}})
    {
	die "Unable to find the node '$digraph->{from}'\n";
    }
    my $from_node = $names{$digraph->{from}};
    $from_node .= "'" if ($digraph->{from_rev});

    # create the starting vertex
    unless ($g->has_vertex($from_node))
    {
	$g->add_vertex($from_node);
    }

    for (my $i = 0; $i < @{$digraph->{to}}; $i++)
    {
	my $connected_node = $digraph->{to}[$i];
	my $reverse = 0;

	if ($connected_node =~ /'/)
	{
	    $reverse = 1;
	    $connected_node =~ s/'//g;
	}

	# the node need to exist
	unless (exists $names{$connected_node})
	{
	    die "Unable to find the node '$connected_node'\n";
	}

	my $node = $names{$connected_node};
	$node .= "'" if ($reverse);

	# create the vertex
	unless ($g->has_vertex($node))
	{
	    $g->add_vertex($node);
	}

	# and the edge
	$g->add_edge($from_node, $node);
    }
}

    # 	# extract Node information
    # 	my @nodes = $line =~ /((?:EDGE|NODE)_[^:;,]+)/g;

    # 	foreach my $node (@nodes)
    # 	{

    # 	    my $nodename = $node;
    # 	    # estimate node basename

    # 	    my $nodebasename = $nodename;
    # 	    $nodebasename =~ s/\'//g;
    # 	    my $rev_nodename = $nodebasename."'";

    # 	    # check if the node is known
    # 	    unless ($g->has_vertex($nodebasename))
    # 	    {
    # 		$g->add_edge($nodebasename, $rev_nodename);
    # 		$g->add_edge($rev_nodename, $nodebasename);
    # 		print STDERR "Created Node: $node ($nodebasename)\n";
    # 	    }
    # 	}

    # # do we have a graph?
    # if ($line =~ /:/)
    # {
    # 	@nodes = split(/:/, $line);
    # 	if (@nodes != 2)
    # 	{
    # 	    die "Something went wrong";
    # 	}
    # 	my ($from, $to) = @nodes;

    # 	# check if $to contains more than one node (delimeter would be an ,)
    # 	my @to_nodes = ($to);
    # 	if ($to =~ /,/)
    # 	{
    # 	    @to_nodes = split(/,/, $to);
    # 	}

    # 	foreach my $to_node (@to_nodes)
    # 	{
    # 	    $g->add_edge($from, $to_node);
    # 	}
    # }


my $seqlen = 0;
foreach my $name (keys %names)
{
    $seqlen += length($seq[$names{$name}]);
}
printf STDERR "Found %d different nodes with %d bp total sequence length and %d digraph information packages\n", (keys %names)+0, $seqlen, @digraphs+0;

foreach my $v ($g->vertices)
{
    #print STDERR "Vertices: $v\n";
}

my @all_weakly_connected_components = $g->weakly_connected_components();

printf STDERR "Found %d weakly connected components\n", @all_weakly_connected_components+0;

my $max = @all_weakly_connected_components+0;
my $progress = Term::ProgressBar->new({name => 'WCC', count => $max, remove => 1, ETA   => 'linear'});
$progress->minor(0);
my $next_update = 0;

for(my $i = 0;  $i < @all_weakly_connected_components+0; $i++)
{
    $next_update = $progress->update($i) if $i >= $next_update;

    my $wcc = $all_weakly_connected_components[$i];

    next if (@{$wcc} < $MINNODES || @{$wcc} > $MAXNODES);

    my $c = subgraph($g, @{$wcc});

    my @V = $c->vertices();

    my $seqlen = 0;
    foreach my $pos (@V)
    {
	$pos =~ s/'//;
	$seqlen += length($seq[$pos]);
    }

    next unless ($c->is_cyclic && $seqlen >= $MINSEQLEN && $seqlen <= $MAXSEQLEN);

#    print Dumper($wcc); use Data::Dumper;

    print "$c", "\n";
}

$progress->update($max) if $max >= $next_update;

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
