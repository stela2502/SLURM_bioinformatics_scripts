#! /usr/bin/perl -w

=head1 LICENCE

  Copyright (C) 2018-03-05 Stefan Lang

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
   
   binCreate.pl from git@github.com:stela2502/Stefans_Lib_Esentials.git commit 70d6934e575eaa13cb8390e5c050eeb0d38bba81
   

=head1  SYNOPSIS

    killThemAll.pl
       -uname          :the username to kill obs from (default stefanl)
     -iterations       :ieach isteration contains 5sec wait - how long should I kill everyting? (default 100)


       -help           :print this help
       -debug          :verbose output
   
=head1 DESCRIPTION

  kill my (stefanl) slurm jobs - all of them every 5 seconds.

  To get further help use 'killThemAll.pl -help' at the comman line.

=cut

use Getopt::Long;
use Pod::Usage;

use strict;
use warnings;

use FindBin;
my $plugin_path = "$FindBin::Bin";

my $VERSION = 'v1.0';


my ( $help, $debug, $database, $iterations, $uname);

Getopt::Long::GetOptions(
	 "-iterations=s"    => \$iterations,
         "-uname=s" => \$uname,
	 "-help"             => \$help,
	 "-debug"            => \$debug
);

my $warn = '';
my $error = '';

unless ( defined $iterations) {
	$iterations = 100;
	$warn .= "Set the iterations to 100 x 5sec\n";
}
unless ( defined $uname) {
	$uname = 'stefanl';
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

$task_description .= 'perl '.$plugin_path .'/killThemAll.pl';
$task_description .= " -iterations '$iterations'" if (defined $iterations);




## Do whatever you want!
my $cmd;
$cmd = "chancel_all_slurm_jobs.pl  -uname $uname";
print "I will run this $iterations times:\n$cmd\n";

unless ( $debug ) {
	for ( my $i = 0; $i < $iterations; $i ++ ){
		system( $cmd );
		sleep(5);
	}
}
print "Finished\n";
