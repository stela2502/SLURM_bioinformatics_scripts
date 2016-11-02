#! /usr/bin/perl -w

=head1 LICENCE

  Copyright (C) 2016-10-18 Stefan Lang

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

    BAM_select_unique_hits.pl
       -infile       :the source bam/sam file
       -outfile      :the traget BAM file
       -max_hits     :the max number of hits (default 1)

       -help           :print this help
       -debug          :verbose output
   
=head1 DESCRIPTION

  Selects uniquely mapped hits based on the hisat NH:i:\d SAM field.

  To get further help use 'BAM_select_unique_hits.pl -help' at the comman line.

=cut

use Getopt::Long;
use Pod::Usage;
use stefans_libs::BAMfile;

use strict;
use warnings;

use FindBin;
my $plugin_path = "$FindBin::Bin";

my $VERSION = 'v1.0';


my ( $help, $debug, $database, $infile, $outfile, $max_hits);

Getopt::Long::GetOptions(
	 "-infile=s"    => \$infile,
	 "-outfile=s"    => \$outfile,
	 "-max_hits=s"   => \$max_hits,
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
unless ( defined $max_hits) {
	$max_hits=1;
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

$task_description .= 'perl '.$plugin_path .'/BAM_select_unique_hits.pl';
$task_description .= " -infile '$infile'" if (defined $infile);
$task_description .= " -outfile '$outfile'" if (defined $outfile);
$task_description .= " -max_hits $max_hits";


open ( LOG , ">$outfile.log") or die $!;
print LOG $task_description."\n";
close ( LOG );


## Do whatever you want!

my $OBJ = stefans_libs::BAMfile -> new();

$OBJ -> select_single_match ( $infile, $max_hits, $outfile );

