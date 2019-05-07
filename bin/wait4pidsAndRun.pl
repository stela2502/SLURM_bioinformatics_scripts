#! /usr/bin/perl -w

=head1 LICENCE

  Copyright (C) 2017-02-08 Stefan Lang

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

    wait4pidsAndRun.pl
       -pids     :the slurm pids to wait for
       -cmd      :the command to run
       -options  :format option space value
       				sleep <time in seconds> default 5


       -help           :print this help
       -debug          :verbose output
   
=head1 DESCRIPTION

  Wait untill a set of slurm jobs have finished and run the cmd.

  To get further help use 'wait4pidsAndRun.pl -help' at the comman line.

=cut

use Getopt::Long;
use Pod::Usage;

use strict;
use warnings;

use stefans_libs::SLURM;


use FindBin;
my $plugin_path = "$FindBin::Bin";

my $VERSION = 'v1.0';


my ( $help, $debug, $database, @pids, $cmd, $options, @options);

Getopt::Long::GetOptions(
       "-pids=s{,}"    => \@pids,
	 "-cmd=s"    => \$cmd,
       "-options=s{,}"    => \@options,

	 "-help"             => \$help,
	 "-debug"            => \$debug
);

my $warn = '';
my $error = '';

unless ( defined $pids[0]) {
	$error .= "the cmd line switch -pids is undefined!\n";
}
unless ( defined $cmd) {
	$error .= "the cmd line switch -cmd is undefined!\n";
}
unless ( defined $options[0]) {
	$warn .= "the cmd line switch -options is undefined!\n";
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

### initialize default options:

$options->{'sleep'} ||= 5;

###


my ( $task_description);

$task_description .= 'perl '.$plugin_path .'/wait4pidsAndRun.pl';
$task_description .= ' -pids "'.join( '" "', @pids ).'"' if ( defined $pids[0]);
$task_description .= " -cmd '$cmd'" if (defined $cmd);
$task_description .= ' -options "'.join( '" "', @options ).'"' if ( defined $options[0]);


for ( my $i = 0 ; $i < @options ; $i += 2 ) {
	$options[ $i + 1 ] =~ s/\n/ /g;
	$options->{ $options[$i] } = $options[ $i + 1 ];
}
###### default options ########
#$options->{'something'} ||= 'default value';
##############################

## Do whatever you want!

my $SLURM = stefans_libs::SLURM->new();

while ( $SLURM->pids_finished(@pids) == 0 ) {
	warn "sleeping for $options->{'sleep'} sec\n" if ($debug);
	sleep($options->{'sleep'});
}
print "staring the downstream process:\n$cmd\n";

system ( $cmd ) unless ($debug );
