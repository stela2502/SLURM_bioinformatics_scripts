#! /usr/bin/perl -w

=head1 LICENCE

  Copyright (C) 2016-10-12 Stefan Lang

  This program is free software; you can redistribute it 
  and/or modify it under the terms of the GNU General Public License 
  as published by the Free Software Foundation; 
  either version 3 of the License, or (at your option) any later version.

  This program is distributed in the hope that it will be useful, 
  but WITHOUT ANY WARRANTY; without even the implied warranty of 
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
  See the GNU General Public License for more details.

  You should have received a copy of the GNU General Public License 
  along with this program; if not, see <http://www.gnu.org/licenses/>.

=head1  SYNOPSIS

    filter_fastq_file.pl
       -infile       :<please add some info!>
       -sequence       :<please add some info!>
       -outfile       :<please add some info!>


       -help           :print this help
       -debug          :verbose output
   
=head1 DESCRIPTION

  Reads one fastq file and removes one read sequence.

  To get further help use 'filter_fastq_file.pl -help' at the comman line.

=cut

use Getopt::Long;
use Pod::Usage;

use strict;
use warnings;

use FindBin;
my $plugin_path = "$FindBin::Bin";

my $VERSION = 'v1.0';

use stefans_libs::FastqFile;

my ( $help, $debug, $database, $infile, $sequence, $outfile);

Getopt::Long::GetOptions(
	 "-infile=s"    => \$infile,
	 "-sequence=s"    => \$sequence,
	 "-outfile=s"    => \$outfile,

	 "-help"             => \$help,
	 "-debug"            => \$debug
);

my $warn = '';
my $error = '';

unless ( defined $infile) {
	$error .= "the cmd line switch -infile is undefined!\n";
}
unless ( defined $sequence) {
	$error .= "the cmd line switch -sequence is undefined!\n";
}
unless ( defined $outfile) {
	$error .= "the cmd line switch -outfile is undefined!\n";
}


if ( $help ){
	print helpString( ) ;
	exit;
}

if ( $error =~ m/\w/ ){
	helpString($error ) ;
	exit;
}

sub helpString {
	my $errorMessage = shift;
	$errorMessage = ' ' unless ( defined $errorMessage); 
	print "$errorMessage.\n";
	pod2usage(q(-verbose) => 1);
}



my ( $task_description);

$task_description .= 'perl '.$plugin_path .'/filter_fastq_file.pl';
$task_description .= " -infile '$infile'" if (defined $infile);
$task_description .= " -sequence '$sequence'" if (defined $sequence);
$task_description .= " -outfile '$outfile'" if (defined $outfile);



open ( LOG , ">$outfile.log") or die $!;
print LOG $task_description."\n";



## Do whatever you want!

my $OBJ = stefans_libs::FastqFile -> new();
$OBJ -> exclude_read ($infile, $sequence, $outfile );
print LOG "$OBJ->{'OK'} OK fastq entries, $OBJ->{'filtered'} filtered.\n";
close ( LOG );
