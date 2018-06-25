#! /usr/bin/perl -w

=head1 LICENCE

  Copyright (C) 2018-06-14 Stefan Lang

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
   
   binCreate.pl from git@gitlab:stefanlang/Stefans_Lib_Esentials.git commit e7545a27e262074f058adf34953b79b184a19248
   

=head1  SYNOPSIS

    runScripts.pl
       -scripts     :the scripts you want to run
       -max_jobs    :the max number of scripts running at a given timepoint (default 10)
       -outfiles    :one outfile for each script (script is not run if outfile exists)


       -help           :print this help
       -debug          :verbose output
   
=head1 DESCRIPTION

  Run a set of already prepared scripts and get additional control over the number of scripts you have in the slurm list while processing the scripts.

  To get further help use 'runScripts.pl -help' at the comman line.

=cut

use Getopt::Long;
use Pod::Usage;

use stefans_libs::SLURM;


use strict;
use warnings;

use FindBin;
my $plugin_path = "$FindBin::Bin";

my $VERSION = 'v1.0';


my ( $help, $debug, $database, @scripts, $max_jobs, @outfiles);

Getopt::Long::GetOptions(
       "-scripts=s{,}"    => \@scripts,
	   "-max_jobs=s"      => \$max_jobs,
       "-outfiles=s{,}"   => \@outfiles,

	 "-help"             => \$help,
	 "-debug"            => \$debug
);

my $warn = '';
my $error = '';

unless ( defined $scripts[0]) {
	$error .= "the cmd line switch -scripts is undefined!\n";
}
unless ( defined $max_jobs) {
	$max_jobs = 10;
}
unless ( defined $outfiles[0]) {
	$error .= "the cmd line switch -outfiles is undefined!\n";
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

$task_description .= 'perl '.$plugin_path .'/runScripts.pl';
$task_description .= ' -scripts "'.join( '" "', @scripts ).'"' if ( defined $scripts[0]);
$task_description .= " -max_jobs '$max_jobs'" if (defined $max_jobs);
$task_description .= ' -outfiles "'.join( '" "', @outfiles ).'"' if ( defined $outfiles[0]);


my $fm = root->filemap( $outfiles[0]);
unless ( -d $fm->{'path'} ) {
	system( "mkdir -p $fm->{'path'}" );
}

open ( LOG , ">$outfiles[0].log") or die $!;
print LOG $task_description."\n";
close ( LOG );


## Do whatever you want!

my $SLURM =stefans_libs::SLURM->new( { 'max_jobs' => $max_jobs });

$SLURM->{'debug'} = $debug;

for (my $i = 0; 4 < @scripts; $i ++ ) {
	$SLURM->runScript( $scripts[$i], $outfiles[$i] );
}

print "All scripts started!\n";






