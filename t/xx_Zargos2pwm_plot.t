#! /usr/bin/perl
use strict;
use warnings;
use stefans_libs::root;
use Test::More tests => 24;
use stefans_libs::flexible_data_structures::data_table;

use FindBin;
my $plugin_path = "$FindBin::Bin";

my ( $value, @values, $exp, $outpath, @options, $infile, );

$infile = "$plugin_path/data/test_zargos_out.txt";

ok ( -f $infile, "test file exists ($infile)" );

$outpath = "$plugin_path/data/output/Zargos2pic";

my $exec = $plugin_path . "/../bin/Zargos2pwm_plot.pl";
ok( -f $exec, 'the script has been found' );
my $outpath = "$plugin_path/data/output/Zargos2pwm_plot";
if ( -d $outpath ) {
	system("rm -Rf $outpath");
}

@options= qw(sampleName data1);

my $cmd =
    "perl -I $plugin_path/../lib  $exec "
. " -outpath " . $outpath 
. " -options " . join(' ', @options )
. " -infile " . $infile 
. " -debug";

system( $cmd );

ok ( -d $outpath, "outpath created ($outpath)" );

if ( -f "$outpath/1_1041.R" ) {
	open ( IN , "$outpath/1_1041.R") or die "could not open the outfile $outpath/1_1041.R\n$!";
	@values = map { chomp; $_ } <IN>;
	#print "\$exp = ".root->print_perl_var_def( \@values ).";\n";
	$exp = [ 
	'library(seqLogo)', 
	'proportion <- function(x){', 
	'   rs <- sum(x);', 
	'   return(x / rs);', 
	'}', 
	'A <- c(213, 84, 204, 108, 195, 286)', 
	'C <- c(66, 328, 36, 56, 461, 52)', 
	'G <- c(539, 541, 9, 32, 25, 283)', 
	'T <- c(223, 88, 792, 845, 360, 420)', 
	'df <- data.frame( A,C,G,T )', 
	'pwm <- apply(df, 1, proportion)', 
	'pwm <- makePWM(pwm)', 
	"pdf( file = '$outpath/1_1041.pdf', width=6, height=6)", 
	'seqLogo(pwm)', 
	'dev.off()' 
	];
	
	is_deeply( \@values,$exp , "R script as expected" );
	close ( IN );
}
else {
	  BAIL_OUT( "No outfile was created as expected!" );
}

foreach ( qw(10_7.pdf  1_1041.pdf  2_1041.R    3_1027.R   4_581.R    5_124.R   6_74.R    7_42.R    8_28.R    9_13.R
10_7.R    1_1041.R    2_1041.pdf                     3_1027.pdf  4_581.pdf  5_124.pdf  6_74.pdf  7_42.pdf  8_28.pdf  9_13.pdf
) ) {
	ok (-f "$outpath/$_", "outfile $_");
}



#print "\$exp = ".root->print_perl_var_def($value ).";\n";