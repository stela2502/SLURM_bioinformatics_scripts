#! /usr/bin/perl -w

#  Copyright (C) 2016-01-26 Stefan Lang

#  This program is free software; you can redistribute it
#  and/or modify it under the terms of the GNU General Public License
#  as published by the Free Software Foundation;
#  either version 3 of the License, or (at your option) any later version.

#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
#  See the GNU General Public License for more details.

#  You should have received a copy of the GNU General Public License
#  along with this program; if not, see <http://www.gnu.org/licenses/>.

=head1 convert_R_S3_R_S4.pl

Defining the same generic function in multiple packages is a problem as loading both packages deletes the ones package functionallity.
Therefore this script checks every unescaped setGeneric() function call with an escaped one.

To get further help use 'convert_R_S4_addCheck2setGeneric.pl -help' at the comman line.

=cut

use Getopt::Long;
use strict;
use warnings;

use stefans_libs::root;

use FindBin;
my $plugin_path = "$FindBin::Bin";

my $VERSION = 'v1.0';

my ( $help, $debug, $database, @R_source, @outfile, $R_source, $outfile );

Getopt::Long::GetOptions(
	"-R_source=s{,}" => \@R_source,
	"-outfile=s{,}"  => \@outfile,

	"-help"  => \$help,
	"-debug" => \$debug
);

my $warn  = '';
my $error = '';

unless ( -f $R_source[0] ) {
	$error .= "the cmd line switch -R_source is undefined!\n";
}
unless ( defined $outfile[0] ) {
	@outfile = @R_source;
	$warn .= "the outfile is set to the infile!\n";
}

if ($help) {
	print helpString();
	exit;
}

if ( $error =~ m/\w/ ) {
	print helpString($error);
	exit;
}

sub helpString {
	my $errorMessage = shift;
	$errorMessage = ' ' unless ( defined $errorMessage );
	return "
 $errorMessage
 command line switches for convert_R_S3_R_S4.pl

   -R_source      :The R file you want to check
   -outfile       :the outfile if different from infile

   -help           :print this help
   -debug          :verbose output
   

";
}

my ($task_description);

$task_description .=
  'perl ' . root->perl_include() . ' ' . $plugin_path . '/convert_R_S3_R_S4.pl';
$task_description .= " -R_source $R_source" if ( defined $R_source );
$task_description .= " -outfile $outfile"   if ( defined $outfile );

open( LOG, ">$outfile[0].log" )
  or die "I could not create the logfile '$outfile[0].log'\n$!\n";
print LOG $task_description . "\n";
close(LOG);

## Do whatever you want!

## I am only interested in top level function setGeneric calls.
#Hence I need to check if there is any top level function call and I need to keep a brackets count.
my $br;

sub update_klammer {
	my $line = shift;
	return if ( $line =~ m/^#/ );
	$br->{'('} += length( $line =~ m/(\()/g );
	$br->{'('} -= length( $line =~ m/(\))/g );
	$br->{'{'} += length( $line =~ m/(\{)/g );
	$br->{'{'} -= length( $line =~ m/(\})/g );
	$br->{'['} += length( $line =~ m/(\[)/g );
	$br->{'['} -= length( $line =~ m/(\])/g );

	#print "line = $line";
	#print "\$br = ".root->print_perl_var_def( $br )."}\n";
}

sub is_setGeneric {
	my ($line) = @_;
	if ( $line =~ m/^\s*setGeneric\(\s*['"]([\w\.]+)["']\s*,/ ) {
		return $1;
	}
	return 0;
}

for ( my $fn = 0 ; $fn < @R_source ; $fn++ ) {
	$R_source = $R_source[$fn];
	$outfile  = $outfile[$fn];

	$br = {
		'{' => 0,
		'[' => 0,
		'(' => 0
	};
	open( IN, "<$R_source" ) or die $!;
	my @file = <IN>;
	close(IN);

	my $fm = root->filemap($R_source);

	#print "\$fm = ".root->print_perl_var_def( $fm).";\n";
	my $functions;

	open( OUT, ">$outfile" ) or die $!;

	my $store = 0;

	my ( $function_store, $function_name, $modified );
	$function_store = '';
	for ( my $i = 0 ; $i < @file ; $i++ ) {
		&update_klammer( $file[$i] );
		if ($store) {
			if ( $br->{'('} == 0 ) {
				## this line finishes the function.
				$function_store .=
				    $file[$i]
				  . "}else {\n"
				  . "\tprint (\"Onload warn generic function '$function_name' already defined - no overloading here!\")\n"
				  . "}\n";
				print OUT $function_store;
				$function_store = '';
				$store          = 0;
				$modified       = 1;
				next;
			}
		}
		elsif ( &is_setGeneric( $file[$i] ) ) {
			$function_name = &is_setGeneric( $file[$i] );
			$file[$i]      = "if ( ! isGeneric('$function_name') ){ $file[$i]";
			$store         = 1;
		}
		if ($store) {
			$function_store .= $file[$i];
		}
		else {
			print OUT $file[$i];
		}
	}

	close(OUT);

	if ($modified) {
		print "changed $outfile\n";
		#print "The setGeneric function has been changed!\n";
	}
	else {
		print "unchanged $outfile\n";
	}

}
