#! /usr/bin/perl -w

#  Copyright (C) 2016-05-31 Stefan Lang

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

=head1 bowtie2_run_aurora.pl

This script takes a  list of fastq or fastq.gz files, created a folder named bowtie2 below these files and outputs the processed data there. Use -debug to not start the SBIN script.

To get further help use 'bowtie2_run_aurora.pl -help' at the comman line.

=cut

use Getopt::Long;
use strict;
use warnings;

use FindBin;
my $plugin_path = "$FindBin::Bin";

my $VERSION = 'v1.0';


my ( $help, $debug, $database, @files, @options, $genome);

Getopt::Long::GetOptions(
	 "-files=s{,}"    => \@files,
	 "-options=s{,}"    => \@options,
	 "-genome=s"    => \$genome,

	 "-help"             => \$help,
	 "-debug"            => \$debug
);

my $warn = '';
my $error = '';

unless ( defined $files[0]) {
	$error .= "the cmd line switch -files is undefined!\n";
}
unless ( defined $options[0]) {
	$error .= "the cmd line switch -options is undefined!\n";
}
unless ( defined $genome) {
	$error .= "the cmd line switch -genome is undefined!\n";
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
 command line switches for bowtie2_run_aurora.pl

   -files       :<please add some info!> you can specify more entries to that
   -options       :<please add some info!> you can specify more entries to that
   -genome       :<please add some info!>

   -help           :print this help
   -debug          :verbose output
   

"; 
}


my ( $task_description);

$task_description .= 'perl '.root->perl_include().' '.$plugin_path .'/bowtie2_run_aurora.pl';
$task_description .= ' -files '.join( ' ', @files ) if ( defined $files[0]);
$task_description .= ' -options '.join( ' ', @options ) if ( defined $options[0]);
$task_description .= " -genome $genome" if (defined $genome);



## Do whatever you want!

