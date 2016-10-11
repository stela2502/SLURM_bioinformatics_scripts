#! /usr/bin/perl
use strict;
use warnings;
use Test::More tests => 2;
BEGIN { use_ok 'stefans_libs::FastqFile' }
use stefans_libs::root;

use FindBin;
my $plugin_path = "$FindBin::Bin";

my ( $value, @values, $exp );
my $OBJ = stefans_libs::FastqFile -> new({'debug' => 1});
is_deeply ( ref($OBJ) , 'stefans_libs::FastqFile', 'simple test of function stefans_libs::FastqFile -> new() ');

#print "\$exp = ".root->print_perl_var_def($value ).";\n";
$value = $OBJ -> open_file ($plugin_path."/data/repeat.fastq.gz" );
is_deeply ( ref($value), 'GLOB', 'opened normal txt file' );
@values = <$value>;
close ( $value );
print "\$exp = ".root->print_perl_var_def([@values[0..3]]).";\n";
is_deeply ( [@values[0..3]], [], "the right sequnce has been read\n" );



