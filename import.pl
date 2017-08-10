#!/usr/bin/env perl

use strict;
use warnings;

use Graph;

my $infile = $ARGV[0];

my $g = Graph->new(directed => 1);

my %names = ();
my @seq = ();
my $lastname = "";

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

	# sequence already known?
	unless (exists $names{$seqname})
	{
	    $lastname = $seqname;

	    push(@seq, "");
	    $names{$seqname} = @seq-1;
	} else {
	    # in case the sequence is already known
	    $lastname = undef;
	}

    } else {
	$seq[$names{$lastname}] .= $line;
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
