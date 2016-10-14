#! /usr/bin/perl -w

=head1 LICENCE

  Copyright (C) 2016-07-05 Stefan Lang

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

    polish_bed_like_files.pl

       -bed_file :the bed file you want to check
       -drop_mt  :chould I drop the MT chromosome?

       -help           :print this help
       -debug          :verbose output
   
=head1 DESCRIPTION

  This script ensures that a given file contains only bona fide chromosomes by adding a chr in front of a otherwise chr - less entry and removing all none chromosome entries

  To get further help use 'polish_bed_like_files.pl -help' at the comman line.

=cut

use Getopt::Long;
use Pod::Usage;
use File::Copy;

use strict;
use warnings;

use FindBin;
my $plugin_path = "$FindBin::Bin";

my $VERSION = 'v1.0';


my ( $help, $debug, $database, $bed_file, $drop_mt);

Getopt::Long::GetOptions(
	 "-bed_file=s"    => \$bed_file,
       "-drop_mt"    => \$drop_mt,

	 "-help"             => \$help,
	 "-debug"            => \$debug
);

my $warn = '';
my $error = '';

unless ( -f $bed_file) {
	$error .= "the cmd line switch -bed_file is undefined!\n";
}
# drop_mt - no checks necessary


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

$task_description .= 'perl '.$plugin_path .'/polish_bed_like_files.pl';
$task_description .= " -bed_file '$bed_file'" if (defined $bed_file);
$task_description .= " -drop_mt " if ($drop_mt);


## Do whatever you want!
 
open ( IN , "<$bed_file" ) or die "I could not open the bed_file $bed_file\n$!\n";
open ( OUT, ">$bed_file.tmp") or die "I could not create the tmp file '$bed_file.tmp'\n$!\n";

my ( @line );
my $OK = { map { 'chr'.$_ => 1 } 1..22,'X', 'Y' };
$OK -> {'chrMT'} = 1 unless ( $drop_mt );

while ( <IN> ) {
	@line= split("\t", $_);
	if ( @line >=2 ) {
		$line[0] = "chr".$line[0] unless ( $line[0] =~ m/^chr/ );
		next unless ( $OK->{$line[0]} );
		print OUT join("\t",@line );
	}else {
		print OUT $_;
	}
}

close ( IN );
close ( OUT );

unless ( $debug ) {
	move( "$bed_file.tmp", $bed_file );
}

