#! /usr/bin/perl
use strict;
use warnings;
use stefans_libs::root;
use Test::More tests => 2;
use stefans_libs::flexible_data_structures::data_table;

use FindBin;
my $plugin_path = "$FindBin::Bin";

my ( $value, @values, $exp, $mailtype, $mailuser, $p, $t, $A, $n, $mempercpu, $N, );

my $exec = $plugin_path . "/../bin/GlobalSlurmOptions.pl";
ok( -f $exec, 'the script has been found' );
my $outpath = "$plugin_path/data/output/GlobalSlurmOptions";
if ( -d $outpath ) {
	system("rm -Rf $outpath");
}

$p = 'lu';
$t = '02:00:00';
$A = 'notImportantHere';
$n = 4;
$N = 1;


my $cmd =
    "perl -I $plugin_path/../lib  $exec "
#. " -mail-type, " . $mailtype
#. " -mail-user, " . $mailuser
. " -p " . $p 
. " -t " . $t 
. " -A " . $A 
. " -n " . $n 
#. " -mem-per-cpu " . $mempercpu 
. " -N " . $N 
#. " -debug"
;
my $start = time;
print $cmd."\n";
system( $cmd );
my $duration = time - $start;
print "Execution time: $duration s\n";


#print "\$exp = ".root->print_perl_var_def($value ).";\n";