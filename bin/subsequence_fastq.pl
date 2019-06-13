#! /usr/bin/perl -w

use stefans_libs::FastqFile;
use stefans_libs::root;
=head1 LICENCE

  Copyright (C) 2017-02-20 Stefan Lang

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

    subsequence_fastq.pl
       -infile      :the input fastq file
       -outfile     :the outfile
       -asFasta     :output fastq not fastq
       
       -start       :start of every sequence (start with 0)
       -length      :length of the fragment (undef == to the end)


       -help           :print this help
       -debug          :verbose output
   
=head1 DESCRIPTION

  A quick way to get a subsequence from ever fastq entry in a file.

  To get further help use 'subsequence_fastq.pl -help' at the comman line.

=cut

use Getopt::Long;
use Pod::Usage;

use strict;
use warnings;

use FindBin;
my $plugin_path = "$FindBin::Bin";

my $VERSION = 'v1.0';

my ( $help, $debug, $database, $infile, $outfile, $asFasta, $start, $length );

Getopt::Long::GetOptions(
	"-infile=s"  => \$infile,
	"-outfile=s" => \$outfile,
	"-asFasta"   => \$asFasta,
	"-start=s"   => \$start,
	"-length=s"  => \$length,

	"-help"  => \$help,
	"-debug" => \$debug
);

my $warn  = '';
my $error = '';

unless ( -f $infile ) {
	$error .= "the cmd line switch -infile is undefined!\n";
}
unless ( defined $outfile ) {
	$error .= "the cmd line switch -outfile is undefined!\n";
}

# fasFasta - no checks necessary
unless ( defined $start ) {
	$error .= "the cmd line switch -start is undefined!\n";
}
unless ( defined $length ) {
	$warn .= "the cmd line switch -length is undefined!\n";
}

if ($help) {
	print helpString();
	exit;
}

if ( $error =~ m/\w/ ) {
	helpString($error);
	exit;
}

sub helpString {
	my $errorMessage = shift;
	$errorMessage = ' ' unless ( defined $errorMessage );
	print "$errorMessage.\n";
	pod2usage( q(-verbose) => 1 );
}

my ($task_description);

$task_description .= 'perl ' . $plugin_path . '/subsequence_fastq.pl';
$task_description .= " -infile '$infile'" if ( defined $infile );
$task_description .= " -outfile '$outfile'" if ( defined $outfile );
$task_description .= " -asFasta " if ($asFasta);
$task_description .= " -start '$start'" if ( defined $start );
$task_description .= " -length '$length'" if ( defined $length );

my $fm = root->filemap($outfile);
mkdir( $fm->{'path'} ) unless ( -d $fm->{'path'} );

open( LOG, ">$outfile.log" ) or die $!;
print LOG $task_description . "\n";
close(LOG);

## Do whatever you want!
my $worker = stefans_libs::FastqFile->new();
open( my $OUT, ">$outfile" )
  or die "I could not create the outfile '$outfile'\n$!\n";

my $restrict;

if ( defined $length ) {
	$restrict = sub {
		my ( $worker, $entry ) = @_;
		$entry->sequence( substr($entry->sequence(),$start,$length) );
		if ( $asFasta ) {
			print $OUT ">".$entry->name()."\n".$entry->sequence()."\n";
		}else {
			$entry->quality( substr($entry->quality(),$start,$length) );
			$entry->write($OUT);
		}
		$entry->clear();
	};
}
else {
	$restrict = sub {
		my ( $worker, $entry ) = @_;
		$entry->sequence( substr($entry->sequence(),$start) );
		if ( $asFasta ) {
			print $OUT ">".$entry->name()."\n".$entry->sequence()."\n";
		}else {
			$entry->quality( substr($entry->quality(),$start) );
			$entry->write($OUT);
		}
		$entry->clear();
	};
}

$worker->filter_multiple_files( $restrict, $infile );

print "Done\n";
