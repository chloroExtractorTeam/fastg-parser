#!/usr/bin/env perl

use strict;
use warnings;

use Graph;
use Term::ProgressBar;

use Getopt::Long qw(:config no_ignore_case);
use File::Temp;

my $infile = "";
my $outfile = "";
my $blastdbfile = '../chloroExtractor-github/data/cds.nr98.fa';

GetOptions(
    'i|infile=s' => \$infile,
    'o|outfile=s' => \$outfile,
    'b|blastdb=s' => \$blastdbfile
    );

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

my @contigs_with_blast_hits = ();

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

    # blast
    my %seen = ();
    foreach (@{$wcc})
    {
	$_ =~ s/'//g;

	$seen{$_}++;
    }

    my ($fh, $filename) = File::Temp::tempfile("tempXXXXX", SUFFIX => ".fa", UNLINK => 1);

    foreach my $seqindex (keys %seen)
    {
	print $fh ">", $seqindex, "\n";
	print $fh $seq[$seqindex], "\n";
    }
    close($fh) || die;

    my $output = qx(tblastx -db $blastdbfile -query $filename -evalue 1e-10 -outfmt 6 -num_alignments 1 -num_threads 4);

    if (length($output) > 0)
    {
	print STDERR "Found hits for cyclic graph: ", $c, "\n";;
	push(@contigs_with_blast_hits, $c);
    }
}

$progress->update($max) if $max >= $next_update;

# check if only a single thing is left
my $chloroplast_seq = "";
if (@contigs_with_blast_hits == 1)
{
    my $c = shift @contigs_with_blast_hits;

    # check the number of nodes
    my %nodes = ();

    # find the node with the highest connectivity. This should be the inverted repead
    my @degree = ();

    foreach my $v ($c->vertices)
    {
	push(@degree, { v => $v, val => $c->in_degree($v)+$c->out_degree($v) });

	my $node = "$v";
	$node =~ s/'//g;

	$nodes{$node}++;
    }

    printf STDERR "Found %d different nodes\n", (keys %nodes)+0;

    @degree = sort { $b->{val} <=> $a->{val} || $b->{v} cmp $a->{v} } @degree;

    my $inverted_repeat = "$degree[0]{v}";

    if ($inverted_repeat =~ /\'/)
    {
	$inverted_repeat =~ s/'//g;
    }

    printf STDERR "The inverted repeat is node number: %d\n", $inverted_repeat;
    delete $nodes{$inverted_repeat};

    # number of nodes only 2?
    if (keys %nodes == 2)
    {
	my ($lsc, $ssc) = (keys %nodes);

	# check if assumption is correct and nodes are assigned correct
	if (length($seq[$lsc]) < length($seq[$ssc]))
	{
	    ($lsc, $ssc) = ($ssc, $lsc);
	}

	printf STDERR "The LSC is node number: %d and the SSC is node number: %d\n", $lsc, $ssc;

	$chloroplast_seq = ">potential_chloroplast_sequence\n";

	# the order of lsc(0), inverted_repeat(1), ssc(2) is
	# 0-1,1-2,2-1' or 0-1',1'-2,2-1 but the orientation of the ssc
	# is guessed due to lack of long reads, but to identify 1-0 or
	# 1'-0 we need to find that edge
	if ($c->has_edge($inverted_repeat, $lsc) || $c->has_edge($lsc, $inverted_repeat."'") || $c->has_edge($lsc."'", $inverted_repeat) || $c->has_edge($lsc."'", $inverted_repeat."'") )
	{
	    $chloroplast_seq .= get_orig_sequence_by_number($lsc).get_orig_sequence_by_number($inverted_repeat."'").get_orig_sequence_by_number($ssc).get_orig_sequence_by_number($inverted_repeat)."\n";
	} else {
	    $chloroplast_seq .= get_orig_sequence_by_number($lsc).get_orig_sequence_by_number($inverted_repeat).get_orig_sequence_by_number($ssc).get_orig_sequence_by_number($inverted_repeat."'")"\n";
	}
    }
}

unless ($chloroplast_seq)
{
    foreach my $subgraphs_with_blast_hits (@contigs_with_blast_hits)
    {
	foreach my $v ($subgraphs_with_blast_hits->vertices)
	{
	    $chloroplast_seq .= ">potential_chloroplast_contig_".$v."\n".get_orig_sequence_by_number($v)."\n";
	}
    }
}

open(FH, ">", $outfile) || die "Unable to open '$outfile' for writing\n";
print FH $chloroplast_seq;
close(FH) || die "Unable to close '$outfile' after writing\n";

# return sequences
sub get_orig_sequence_by_number
{
    my ($number) = @_;

    my ($name, $sequence);

    my $reverse = 0;
    if ($number =~ /'/)
    {
	$number =~s/'//;
	$reverse = 1;
    }

    unless (exists $seq2seqname{$number})
    {
	die sprintf "Unable to find original name for sequence #%d\n", $number;
    }

    $name = $seq2seqname{$number};

    unless (defined $seq[$number])
    {
	die sprintf"Unable to find sequence information for #%d\n", $number;
    }

    $sequence = $seq[$number];
    if ($reverse)
    {
	$sequence = reverse $sequence;
	$sequence =~ tr/ACGTacgt/TGCAtgca/;
    }

    # return {name => $name, seq => $sequence};
    return $sequence;
}


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
