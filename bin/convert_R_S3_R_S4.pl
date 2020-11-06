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

Help in converting a R S3 class into a R S4 class by converting the function definitions into generic methods.

To get further help use 'convert_R_S3_R_S4.pl -help' at the comman line.

=cut

use Getopt::Long;
use strict;
use warnings;

use File::Spec;
use stefans_libs::root;

use FindBin;
my $plugin_path = "$FindBin::Bin";

my $VERSION = 'v1.0';

my ( $help, $debug, $database, $R_source,$classN, $outfile );

Getopt::Long::GetOptions(
	"-R_source=s" => \$R_source,
	"-outfile=s"  => \$outfile,
	"-className=s" => \$classN,

	"-help"  => \$help,
	"-debug" => \$debug
);

my $warn  = '';
my $error = '';

unless ( defined $R_source ) {
	$error .= "the cmd line switch -R_source is undefined!\n";
}
unless ( defined $outfile ) {
	$error .= "the cmd line switch -outfile is undefined!\n";
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

   -R_source    :the old function based lib file
   -outfile     :the new S4isizde lib file
   -className   :the class name for this file (in case the file name is not the class name)

   -help           :print this help
   -debug          :verbose output
   

";
}

my ($task_description);

$task_description .=
  'perl ' . root->perl_include() . ' ' . $plugin_path . '/convert_R_S3_R_S4.pl';
$task_description .= " -R_source $R_source" if ( defined $R_source );
$task_description .= " -outfile $outfile"   if ( defined $outfile );
$task_description .= " -className $classN" if ( defined $classN);

open( LOG, ">$outfile.log" ) or die $!;
print LOG $task_description . "\n";
close(LOG);

## Do whatever you want!

open( IN, "<$R_source" ) or die $!;
my @file = <IN>;
close(IN);

my $fm = root->filemap($R_source);
print "\$fm = " . root->print_perl_var_def($fm) . ";\n";
my $functions;
open( OUT, ">$outfile" ) or die $!;
## populate functions
my ( $funN, $funArgs, $add, $bracket, $save, $man );
$add     = 1;
$bracket = 0;
$classN ||= $fm->{'filename_core'};
for ( my $i = 0 ; $i < @file ; $i++ ) {

	if ( $bracket == 0 ) {
		if ( $file[$i] =~ m/setClass/) {
			$man = &populate_man($classN , $classN, $man, '');
			print OUT man2str( $man );
			$man = undef;
			print OUT $file[$i];
		}
		if ( $file[$i] =~ m/^#' ?/ ) {
			print "123: I have the old man string $file[$i]";
			if ( $file[$i] =~ m/^#'\s*\@([\d\w]+)\s+(.*)/ ) {
				$man->{$1} ||= [];
				push( @{ $man->{$1} }, $2 );
			}
			else {
				$file[$i] =~ m/^#'\s*(.*)/;
				$man->{'description'} ||= [];
				push( @{ $man->{'description'} }, $1 );

 #warn "The collected descritpion: ". join(' ',  @{$man->{'description'}})."\n";
			}
		}
		if (   $file[$i] =~ m/([\w_\.\d]*)\s*=\s*function\s*(\(.*\)?)\s*{*/
			|| $file[$i] =~ m/([\w_\.\d]*)\s*<-\s*function\s*(\(.*\)?)\s*{*/ )
		{
			$funN    = $1;
			$funArgs = $2;
			$funN =~ s/\.$fm->{'filename_core'}$//;
			unless ( $file[$i] =~ m/\)\s*{/ ) {
				$add = 1;
				while ( !$file[ $i + $add ] =~ m/(.*\))\s*{/ ) {
					chomp( $file[ $i + $add ] );
					$funArgs .= $file[ $i + $add ];
					$file[ $i + $add ] = '';
					$add++;
				}
				$funArgs .= $1;
				$file[ $i + $add ] =~ s/$1//;
			}
			$funArgs =~ s/\)\s*{/\)/;

			# now I have function name and args in two variables.
			$bracket =
			  1;    ## I need to check whether the } does close the function!
			print "# before populate_man\n\t\$man = " . root->print_perl_var_def( $man ) . ";\n";
			
			$man = &populate_man( $funN, $classN, $man, $funArgs);
			print OUT man2str( $man );
			$man = undef;
			print OUT "if ( ! isGeneric('$funN') ){setGeneric('$funN', ## Name
	function $funArgs { 
		standardGeneric('$funN')
	}
) }

setMethod('$funN', signature = c ('$classN'),
	definition = function $funArgs {\n";
		}
	}
	else {
		if ( $bracket > 0 ) {
			unless ( $file[$i] =~ m/^\s*#/) {
			map { $bracket++ } $file[$i] =~ m/[\[\({]/g;
			map { $bracket-- } $file[$i] =~ m/[\]\)}]/g;
			}
			if ( $bracket == 0 ) {
				chomp( $file[$i] );
				print OUT $file[$i] . " )\n";
			}
			else {
				print OUT $file[$i];
			}
		}
		else {
			print OUT $file[$i];
		}

	}
}
close(OUT);
print "Done\n";

sub man2str {
	my ($man) = @_;
	my $str = '';
	my $n = { 'description' => 1, 'export' => 1 };

	foreach (
		map {
			unless ( $n->{$_} ) { $n->{$_} = 1; $_; }
		} 'name',
		'aliases',
		'rdname',
		'docType'
	  )
	{
		next unless ( defined $_ );
		next unless ( defined $man->{$_} );
		foreach my $v ( @{ $man->{$_} } ) {
			$str .= "#' \@$_ $v\n";
		}
	}

	# first the description
	$man->{'description'} = join( " ", @{ $man->{'description'} } );

	unless ( $man->{'description'} =~ m/\w/ ) {
		$str .= "#' \@description \n";
	}
	else {
		$man->{'description'} =~ s/\s\s+/ /g;
		my $tmp = '';
		foreach ( split( " ", $man->{'description'} ) ) {
			$tmp .= " $_";
			if ( length($tmp) > 80 ) {
				$str .= "#' \@description $tmp\n";
				$tmp = '';
			}
		}
		if ( $tmp =~ m/\w/ ) {
			$str .= "#' \@description $tmp\n";
		}
	}
	foreach (
		map {
			unless ( $n->{$_} ) { $n->{$_} = 1; $_; }
		} 'param',
		'return',
		'example',
		keys %$man
	  )
	{
		next unless ( defined $_ );
		next unless ( defined $man->{$_} );
		foreach my $v ( @{ $man->{$_} } ) {
			$str .= "#' \@$_ $v\n";
		}
	}
	foreach ('export') {
		next unless ( defined $_ );
		next unless ( defined $man->{$_} );
		foreach my $v ( @{ $man->{$_} } ) {
			$str .= "#' \@$_ $v\n";
		}
	}
	return $str;
}

sub populate_man {
	my ( $funN, $classN, $man, $funArgs ) = @_;
	my $skel = {
		'name'        => [$funN],
		'title'       => [ 'description of function ' . $funN ],
		'aliases'     => ["$funN,$classN-method"],
		'docType'     => ['methods'],
		'rdname'      => ["$funN-methods"],
		'description' => [''],
		'param'       => [],
		'export'      => [''],
	};
	foreach ( keys %$skel ) {
		$man->{$_} = $skel->{$_} unless ( defined $man->{$_} );
	}
	$funArgs =~ s/^\s*\(//;
	$funArgs =~ s/\s*\)\s*$//;
	print $funArgs."\n";
	my ($d, @tmp);
	
	my $definedParam = { map { @tmp = split(/\s+/, $_); "$tmp[0]" => 1 } @{$man->{'param'}} };
	
	foreach (split( ",", $funArgs ) ) {
		next if ( $definedParam->{$_} );
		$_ =~ s/=(.*)$//;
		$d = $1;
		$_ =~ s/\s//g;
		if ( $d ) {
			push ( @{$man->{'param'}} , "$_  TEXT MISSING default=$d");
		}else {
			push ( @{$man->{'param'}} , "$_  TEXT MISSING");
		}
	}
	if ( scalar(@{$man->{'param'}}) == 0) {
		@{$man->{'param'}} = [''];
	}
	return $man;
}

sub next_function_name {
	my ( $array, $id ) = @_;
	my $r;
	for ( my $i = $id ; $i < @$array ; $i++ ) {
		$r = &function_name( @$array[$i] );
		return $r if ( defined $r );
	}
	Carp::confess("No function stat found after line $id\n");
}

sub function_name {
	my ($string) = @_;
	if ( $string =~ m/setMethod\(\s*["']([\w\d\.\_]+)["']/ ) {
		return $1;
	}
	if ( $string =~ m/setGeneric\(\s*["']([\w\d\.\_]+)["']/ ) {
		return $1;
	}
	return undef;
}
