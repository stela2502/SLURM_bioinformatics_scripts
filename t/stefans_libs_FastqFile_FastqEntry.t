#! /usr/bin/perl
use strict;
use warnings;
use Test::More tests => 7;
use stefans_libs::root;
BEGIN { use_ok 'stefans_libs::FastqFile::FastqEntry' }

BEGIN { use_ok 'stefans_libs::FastqFile' }

use FindBin;
my $plugin_path = "$FindBin::Bin";

my ( $value, @values, $exp );
my $OBJ = stefans_libs::FastqFile::FastqEntry->new( { 'debug' => 1 } );
is_deeply(
	ref($OBJ),
	'stefans_libs::FastqFile::FastqEntry',
	'simple test of function stefans_libs::FastqFile::FastqEntry -> new() '
);

$OBJ->quality(
'!"#$%&\'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\]^_`abcdefghijklmnopqrstuvwxyz{|}~'
);
@values = $OBJ->quality( undef, 1 );

is_deeply( [ @values ], [ 0 .. 93 ],  "quality score from 0 to 93" );

$OBJ = stefans_libs::FastqFile->new( { 'debug' => 1 } );

is_deeply( ref($OBJ), 'stefans_libs::FastqFile',
	  'simple test of function stefans_libs::FastqFile -> new() ' );

@values = ();
my $get_objects = sub {
	  my ( $obj, $entry ) = @_;
	  push( @values, $entry->copy() );
	  $entry->clear();
  };

my $infile = "$plugin_path/data/Quality.fastq.gz";

$OBJ->filter_multiple_files( $get_objects, $infile );

ok( scalar(@values) == 10, "I got the 10 fastq entries" );

print join( "-", split( "", $values[0]->quality() ) ) . "\n";
print join( "-", $values[0]->quality( undef, 1 ) ) . "\n";

my @sequences = map{ substr($_->sequence,0,16) } @values;

#print join("\n",@sequences)."\n";
#$values[0]->{'debug'} = 1;
#print "\$exp = ".root->print_perl_var_def( {res => $values[0]->distance_to( $sequences[0], 0, 16 ), from => $sequences[0], to=>  $sequences[0] } ).";\n";



#print "\$exp = ".root->print_perl_var_def( [ $values[0]->distance_to( \@sequences, 0, 16 ) ] ) .";\n";
$exp = [ {
  'dist' => '0',
  'mm' => '0'
}, {
  'dist' => '441',
  'mm' => '13'
}, {
  'dist' => '443',
  'mm' => '13'
}, {
  'dist' => '333',
  'mm' => '10'
}, {
  'dist' => '356',
  'mm' => '11'
}, {
  'dist' => '392',
  'mm' => '12'
}, {
  'dist' => '396',
  'mm' => '12'
}, {
  'dist' => '405',
  'mm' => '12'
}, {
  'dist' => '396',
  'mm' => '12'
}, {
  'dist' => '468',
  'mm' => '14'
} ];
is_deeply(  [ $values[0]->distance_to( \@sequences, 0, 16 ) ], $exp, "right distance results" );


#print "\$exp = ".root->print_perl_var_def($value ).";\n";



