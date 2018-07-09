#!/usr/bin/env perl

use strict;
use warnings;

use LWP::Simple;
use Digest::MD5;

use Test::More;
use Test::Script::Run;

# prepare run command
my $input_file = "data/SRR1503755_assembly_graph.fastg";
my $output_file = "data/005.out";
my $expected_output_md5 = "dd3b6d3641d0df9b01d40a953489c1ff";

my @arg = ("-i", $input_file, "-o", $output_file, '-b', 'data/cds.nr98.fa');

my ($ret, $stdout, $stderr ) = run_script("import.pl", \@arg);
my $error_code = Test::Script::Run::last_script_exit_code();

diag($stdout);
diag($stderr);
diag("Error code: $error_code and return value $ret");

# calculate MD5SUM
my $md5 = Digest::MD5->new;
open(my $fh, "<", $output_file) || die "Unable to open file '$output_file' $!\n";
$md5->addfile($fh);
close($fh) || die "Unable to close file '$output_file': $!\n";

my $md5hex = $md5->hexdigest();

is($error_code, 0, 'Run of import returns 0 as error code');
ok(-e $output_file, 'Output file was created');
is($md5hex, $expected_output_md5, 'Output file contains expected result');
like($stderr, qr/\[fcg\] Single circular chloroplast seems to be found/, 'Expected output is part of the log');

done_testing;
