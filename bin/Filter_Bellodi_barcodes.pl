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

    Filter_Bellodi_barcodes.pl
       -infile         :the fastq file to process (gz also possible)
       -outfile        :the base outfile - please make sure the path exists
       -adapter        :the adapter sequence to trom from the end of the read (min 4 bp to trim)
       -sample_barcode :a space separated list of sample barcodes to analze
       -match_start    :relax the matching algorithm to have <match_start> possible N's before the pattern match (defualt = 0)

       -help           :print this help
       -debug          :verbose output
   
=head1 DESCRIPTION

  The strategy was a one read barcode with this setup NNN_Barcode_NN_cDNA_Adapter(Which may or may not be present) - and that is what the tool does support. Flexible parts: Barcode and Adapter.

  To get further help use 'Filter_Bellodi_barcodes.pl -help' at the comman line.

=cut

use Getopt::Long;
use Pod::Usage;

use strict;
use warnings;

use FindBin;
my $plugin_path = "$FindBin::Bin";

my $VERSION = 'v1.0';

use stefans_libs::FastqFile;

my ( $help, $debug, $database, $infile, $outfile, $adapter, $matchStart, @sample_barcode);

Getopt::Long::GetOptions(
	 "-infile=s"    => \$infile,
	 "-outfile=s"    => \$outfile,
	 "-adapter=s"    => \$adapter,
     "-sample_barcode=s{,}"    => \@sample_barcode,
     "-match_start=s"  => \$matchStart,  

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
unless ( defined $adapter) {
	$error .= "the cmd line switch -adapter is undefined!\n";
}
unless ( defined $sample_barcode[0]) {
	$error .= "the cmd line switch -sample_barcode is undefined!\n";
}
unless ( defined $matchStart ){
	$matchStart = 0;
	warn "match_start set to 0";
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

$task_description .= 'perl '.$plugin_path .'/Filter_Bellodi_barcodes.pl';
$task_description .= " -infile '$infile'" if (defined $infile);
$task_description .= " -outfile '$outfile'" if (defined $outfile);
$task_description .= " -adapter '$adapter'" if (defined $adapter);
$task_description .= ' -sample_barcode "'.join( '" "', @sample_barcode ).'"' if ( defined $sample_barcode[0]);
$task_description .= " -match_start $matchStart";


open ( LOG , ">$outfile.log") or die $!;
print LOG $task_description."\n";



## Do whatever you want!

my $OBJ = stefans_libs::FastqFile -> new();
foreach my $barcode ( @sample_barcode ) {
	$OBJ -> extract_barcode ( $infile,  $adapter, $barcode, $matchStart, $outfile.$barcode.".fastq" );
	print LOG "$barcode: $OBJ->{'OK'} OK fastq entries, $OBJ->{'filtered'} filtered.\n";
	
}

close ( LOG );