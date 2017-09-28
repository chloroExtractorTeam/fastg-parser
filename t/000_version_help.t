#!/usr/bin/env perl

use strict;
use warnings;

use LWP::Simple;
use Digest::MD5;

use Test::More;
use Test::Script::Run;

# prepare run command for version
my @arg = ("-v");

my ($ret, $stdout, $stderr ) = run_script("import.pl", \@arg);
my $error_code = Test::Script::Run::last_script_exit_code();

is($error_code, 0, 'Run of import returns 0 as error code');
like($stdout, qr/^v\d+\.\d+\.\d/, 'Version number is retured');

# prepare run command for help
@arg = ("-h");

($ret, $stdout, $stderr ) = run_script("import.pl", \@arg);
$error_code = Test::Script::Run::last_script_exit_code();

is($error_code, 1, 'Run of import.pl -h returns 0 as error code');

done_testing;
