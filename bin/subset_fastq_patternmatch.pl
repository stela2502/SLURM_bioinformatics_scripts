#! /usr/bin/perl -w

=head1 LICENCE

  Copyright (C) 2016-10-17 Stefan Lang

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

    subset_fastq_patternmatch.pl
       -infile     :the fastq file to subset
       -outfile    :the new fastq file
       -pattern    :a perl pattern to search for
       -where      :0 -> header; 1 -> sequence; 3-> quality scores


       -help           :print this help
       -debug          :verbose output
   
=head1 DESCRIPTION

  Use a perl pattern match to subselect fastq entries.

  To get further help use 'subset_fastq_patternmatch.pl -help' at the comman line.

=cut

use Getopt::Long;
use Pod::Usage;

use strict;
use warnings;

use stefans_libs::FastqFile;
use stefans_libs::BAMfile;

use FindBin;
my $plugin_path = "$FindBin::Bin";

my $VERSION = 'v1.0';


my ( $help, $debug, $database, $infile, $outfile, $pattern, $where);

Getopt::Long::GetOptions(
	 "-infile=s"    => \$infile,
	 "-outfile=s"    => \$outfile,
	 "-pattern=s"    => \$pattern,
	 "-where=s"    => \$where,

	 "-help"             => \$help,
	 "-debug"            => \$debug
);

my $warn = '';
my $error = '';

unless ( defined $infile) {
	$error .= "the cmd line switch -infile is undefined!\n";
}
unless ( defined $outfile) {
	$error .= "the cmd line switch -outfile is undefined!\n";
}
unless ( defined $pattern) {
	$error .= "the cmd line switch -pattern is undefined!\n";
}
unless ( defined $where) {
	$error .= "the cmd line switch -where is undefined!\n";
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

$task_description .= 'perl '.$plugin_path .'/subset_fastq_patternmatch.pl';
$task_description .= " -infile '$infile'" if (defined $infile);
$task_description .= " -outfile '$outfile'" if (defined $outfile);
$task_description .= " -pattern '$pattern'" if (defined $pattern);
$task_description .= " -where '$where'" if (defined $where);



open ( LOG , ">$outfile.log") or die $!;
print LOG $task_description."\n";
close ( LOG );

my $OBJ;
## Do whatever you want!
if ( $infile =~ m/.fa?s?t?q$/ ){
	$OBJ = stefans_libs::FastqFile -> new();
}elsif( $infile =~ m/.[sb]am$/){
	$OBJ = stefans_libs::BAMfile -> new();
}
$OBJ -> select_4_str ( $infile, $pattern, $where, $outfile );
