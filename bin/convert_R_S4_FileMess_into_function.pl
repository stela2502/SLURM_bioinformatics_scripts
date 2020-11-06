#! /usr/bin/perl -w

#  Copyright (C) 2016-04-15 Stefan Lang

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

=head1 convert_R_S4_FileMess_into_function.pl

this script converts a set of R functions from probabaly multiple classes into a function name based set of files.

To get further help use 'convert_R_S4_FileMess_into_function.pl -help' at the comman line.

=cut

use Getopt::Long;
use strict;
use warnings;

use stefans_libs::root;

use FindBin;
my $plugin_path = "$FindBin::Bin";

my $VERSION = 'v1.0';


my ( $help, $debug, $database, @infiles, $use_fname, $outpath);

Getopt::Long::GetOptions(
	 "-infiles=s{,}"    => \@infiles,
	 "-outpath=s"    => \$outpath,
	 "-use_fname" => \$use_fname,

	 "-help"             => \$help,
	 "-debug"            => \$debug
);

my $warn = '';
my $error = '';

unless ( defined $infiles[0]) {
	$error .= "the cmd line switch -infiles is undefined!\n";
}
unless ( defined $outpath) {
	$error .= "the cmd line switch -outpath is undefined!\n";
}


if ( $help ){
	print helpString( ) ;
	exit;
}

if ( $error =~ m/\w/ ){
	print helpString($error ) ;
	exit;
}

sub helpString {
	my $errorMessage = shift;
	$errorMessage = ' ' unless ( defined $errorMessage); 
 	return "
 $errorMessage
 command line switches for convert_R_S4_FileMess_into_function.pl

   -infiles       :a list of R files that you want to splice
   -outpath       :the outpath for the new files
   
   -use_fname     :should I include the original filename in the outfile

   -help           :print this help
   -debug          :verbose output
   

"; 
}


my ( $task_description);

$task_description .= 'perl '.$plugin_path .'/convert_R_S4_FileMess_into_function.pl';
$task_description .= ' -infiles '.join( ' ', @infiles ) if ( defined $infiles[0]);
$task_description .= " -outpath $outpath" if (defined $outpath);


mkdir( $outpath ) unless ( -d $outpath );
open ( LOG , ">$outpath/".$$."_convert_R_S4_FileMess_into_function.pl.log") or die $!;
print LOG $task_description."\n";
close ( LOG );


## Do whatever you want!

my ($functions, $fm); ## here I will store all function definitions untill I create the outfiles

foreach my $file ( @infiles ) {
	open ( IN, "<$file" ) or Carp::confess ( "I could not open the file '$file'\n$!\n");
	$fm = root->filemap($file);
	my $this_function ='';
	my ( $is_func, $last_comment );
	$is_func = $last_comment = 0;
	while ( <IN> ) {
		## all functions have a leading documentation and therefore I can search for the #' string
		if ( $_ =~ m/^#'/) {
			if ( $last_comment==0 && $this_function =~ m/\w/ ) { ## I need to save the last bit
				$is_func = &funcname($this_function);
				print ( "$this_function\n has the funciton name $is_func\n" );
				$functions->{$is_func} ||= Func->new($fm->{'filename_base'});
				$functions->{$is_func}->R($this_function);
				$this_function = '';
				$last_comment= 1;
				$is_func = 1;
			}
			else {
				$is_func = 1;
				$last_comment= 1;
			}
		}
		else {
			$last_comment= 0;
		}
		if ( $is_func ) {
			$this_function .= $_;
		}
	}
	if ( length($this_function) > 0 ) {
		$is_func = &funcname($this_function);
		$functions->{$is_func} ||= Func->new($fm->{'filename_base'});
		$functions->{$is_func}->R($this_function);
		$this_function = '';
	}
	close ( IN );
}

## and now all function definitions need to be saved as files
my $outfile;
foreach my $fname ( keys %$functions ) {
	if ( $use_fname ){
		$outfile = "$outpath/$functions->{$fname}->{'file'}"."_$fname.R";
	}
	else {
		$outfile = "$outpath/$fname.R";
	}
	unless ( -f $outfile ) {
		open ( OUT, ">$outfile" ) or die "I could not create the outfile '$outfile'\n$!\n";
		print OUT $functions->{$fname}->R();
		close ( OUT );
		print "created file '$outfile'\n";
	}
}

print "Done\n";

sub funcname {
	my ( $string ) = @_;
	my $funcname = '01_class';
	if ( $string =~ m/setMethod\s*\(\s*['"]([\w_\.]+)["']/ ){
		$funcname = $1;
	}
	return $funcname
}


package Func;

sub new{
	my ( $class, $file ) = @_;

	my ($self);
	$self = {
		'file' => $file,
		'r' => ''
	};

	bless( $self, $class ) if ( $class eq "Func" );

	return $self;
}

sub R {
	my ($self, $add ) = @_;
	$add ||= '';
	$self->{'R'} .= $add;
	return $self->{'R'};
}

