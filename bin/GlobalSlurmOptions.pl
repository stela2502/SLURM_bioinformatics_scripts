#! /usr/bin/perl -w

=head1 LICENCE

  Copyright (C) 2017-06-26 Stefan Lang

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
   
   binCreate.pl from  commit 
   

=head1  SYNOPSIS

    GlobalSlurmOptions.pl
       -n       :the default slurm number of cores per processor
       -N       :the default slurm number of processors per script
       -t       :the default slurm max time per process
       
       -A       :the default project to run SLURM scripts under (optional)
       -p       :the default partition to use (optional)
       -mail-user,       :the default mail user to use (optional)
       -mail-type,       :the default mail type to use (optional)
       -mem-per-cpu      :the default mem-per-cpu (optional and probaly not working!)


       -help           :print this help
       -debug          :verbose output
   
=head1 DESCRIPTION

  This script uses the stefans_libs::flexible_data_structures::optionsFile to store all required/supported SLURM options in the file ~/.SLURM_options.txt.

  To get further help use 'GlobalSlurmOptions.pl -help' at the comman line.

=cut

use Getopt::Long;
use Pod::Usage;

use strict;
use warnings;

use stefans_libs::SLURM;

use FindBin;
my $plugin_path = "$FindBin::Bin";

my $VERSION = 'v1.0';


my ( $help, $debug, $database, $n, $N, $t, $A, $p, $mailuser, $mailtype, $mempercpu);
#
#Getopt::Long::GetOptions(
#	 "-n=s"    => \$n,
#	 "-N=s"    => \$N,
#	 "-t=s"    => \$t,
#	 "-A=s"    => \$A,
#	 "-p=s"    => \$p,
#	 "-mail-user=s"    => \$mailuser,
#	 "-mail-type=s"    => \$mailtype,
#	 "-mem-per-cpu=s"    => \$mempercpu,
#
#	 "-help"             => \$help,
#	 "-debug"            => \$debug
#);
#
my $warn = '';
my $error = '';
#
#unless ( defined $n) {
#	$error .= "the cmd line switch -n is undefined!\n";
#}
#unless ( defined $N) {
#	$error .= "the cmd line switch -N is undefined!\n";
#}
#unless ( defined $t) {
#	$error .= "the cmd line switch -t is undefined!\n";
#}


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

$task_description .= 'perl '.$plugin_path .'/GlobalSlurmOptions.pl';
#$task_description .= " -n '$n'" if (defined $n);
#$task_description .= " -N '$N'" if (defined $N);
#$task_description .= " -t '$t'" if (defined $t);
#$task_description .= " -A '$A'" if (defined $A);
#$task_description .= " -p '$p'" if (defined $p);
#$task_description .= " -mail-user, '$mailuser,'" if (defined $mailuser);
#$task_description .= " -mail-type, '$mailtype,'" if (defined $mailtype);
#$task_description .= " -mem-per-cpu '$mempercpu'" if (defined $mempercpu);


## Do whatever you want!

my $obj = stefans_libs::SLURM->new()->{options};


$obj->load(undef,1);# load and do die if there are errors in the config file!
## but this will create a usable file ;-)

print "Config is theoretically OK\n"."Check ".$obj->{'default_file'}."\n";
