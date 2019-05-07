#! /usr/bin/perl
use strict;
use warnings;
use Test::More tests => 3;
BEGIN { use_ok 'stefans_libs::file_readers::transfaq::motive' }

use FindBin;
my $plugin_path = "$FindBin::Bin";


my ($test_object, $value, $exp, @values);

$test_object =  stefans_libs::file_readers::transfaq::motive-> new();
is_deeply ( ref($test_object) , 'stefans_libs::file_readers::transfaq::motive', 'simple test of function stefans_libs::file_readers::transfaq::motive -> new()' );

$value = $test_object->AddDataset ( {             'A' => 0,
            'C' => 1,
            'G' => 2,
            'T' => 3, } );
is_deeply( $value, 1, "we could add a sample dataset");
$test_object =  stefans_libs::file_readers::transfaq::motive-> new( );
$test_object->read_file ("$plugin_path/data/motive.xls" );
#print "\$exp = ".root->print_perl_var_def($test_object -> {'data'}).";\n";
$exp = [ 
[ '213', '66', '539', '223' ],
[ '84', '328', '541', '88' ], 
[ '204', '36', '9', '792' ], 
[ '108', '56', '32', '845' ], 
[ '195', '461', '25', '360' ], 
[ '286', '52', '283', '420' ] ];

is_deeply( $test_object -> {'data'}, $exp, "load from table file" );

ok ( $test_object->score_seq ('GGtTcT' )== 1, "score perfect match ". $test_object->score_seq ('GGtTCT' ) );

ok ( $test_object->score_seq ('CAGGGC' ) < 0.07448582545859, "score worst match " );


## A handy help if you do not know what you should expect
#print "\$exp = ".root->print_perl_var_def($value ).";\n";
