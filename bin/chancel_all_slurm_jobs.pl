#! /usr/bin/perl -w

=head1 LICENCE

  Copyright (C) 2017-12-12 Stefan Lang

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

=head1 CREATED BY
   
   binCreate.pl from git@github.com:stela2502/Stefans_Lib_Esentials.git commit 2b3dac48125abf3cf3d1c692eb96a10f20faf457
   

=head1  SYNOPSIS

    chancel_all_slurm_jobs.pl
       -uname       :your username (if you want to stop all)
       -start       :the start pid to stop (if you want to stop only a subset)
       -stop        :the end pid to stop (if you want to stop only a subset)


       -help           :print this help
       -debug          :verbose output
   
=head1 DESCRIPTION

  Chancels all slurm jobs for a given user

  To get further help use 'chancel_all_slurm_jobs.pl -help' at the comman line.

=cut

use Getopt::Long;
use Pod::Usage;

use strict;
use warnings;

use FindBin;
my $plugin_path = "$FindBin::Bin";

my $VERSION = 'v1.0';

my ( $help, $debug, $database, $uname, $start, $stop );

Getopt::Long::GetOptions(
	"-uname=s" => \$uname,
	"-start=s" => \$start,
	"-stop=s"  => \$stop,

	"-help"  => \$help,
	"-debug" => \$debug
);

my $warn  = '';
my $error = '';

if ( !defined $uname ) {
	$error .= "the cmd line switch -uname is undefined!\n";
	if ( defined $start and defined $stop ) {
		$error = "";
	}
}

# start - no checks necessary
# stop - no checks necessary

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

$task_description .= 'perl ' . $plugin_path . '/chancel_all_slurm_jobs.pl';
$task_description .= " -uname '$uname'" if ( defined $uname );
$task_description .= " -start " if ( defined $start );
$task_description .= " -stop " if ( defined $stop );

## Do whatever you want!

open( IN, " squeue -u $uname |" );
my (@tmp, $min, $max,  $ip );
$min = 1e+9; 
$max = 0;
while (<IN>) {
	$_ =~ s/^\s*//;
	chomp;
	if ($_ =~ m/^(\d+) / ){
		$ip = $1; 
		push ( @tmp, $ip);
		$max = $ip if ( $max < $ip );
		$min = $ip if ( $min > $ip );
	}
}

close ( IN );

$start ||=$min;
$stop ||= $max;

foreach $ip ( @tmp ) {
	if ( $ip >= $start and $ip <= $stop ){
		print "scancel $ip\n" ;
		system ( "scancel $ip") unless ( $debug );
	}
}
