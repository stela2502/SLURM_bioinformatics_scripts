#! /usr/bin/perl -w

=head1 LICENCE

  Copyright (C) 2016-06-16 Stefan Lang

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

    runCommand.pl
       -cmd      :The command you want to run
       -max_jobs :The maximum number of jobs to the SLURM system for this user (default 40)
       -options  :format: key_1 value_1 key_2 value_2 ... key_n value_n
       		n  :amount of cores per node
       		N  :amount of nodes
       		A  :the accout to use (lu2016-2-7)
       -outfile :the outfile that will be created during the run (to check whether it should be run or not.)

       -loadModules : a list of modules to load in the script
       
       -I_have_loaded_all_modules :add this switch (no value) and mean it ;-)

       -help           :print this help
       -debug          :verbose output
   
=head1 DESCRIPTION

  Gets one command line and runs it on the SLURM backend.
  
  If the string "$SNIC_TMP"" is found in the cmd all files 
  from this folder are copied back to the output folder.

  To get further help use 'runCommand.pl -help' at the comman line.

=cut

use Getopt::Long;
use Pod::Usage;

use strict;
use warnings;

use stefans_libs::SLURM;

use FindBin;
my $plugin_path = "$FindBin::Bin";

my $VERSION = 'v1.0';


my ( $help, $debug, $database, $cmd, $options, @options, $I_have_loaded_all_modules, $outfile, $max_jobs, @loadModules);

Getopt::Long::GetOptions(
	 "-cmd=s"    => \$cmd,
       "-options=s{,}"    => \@options,
	 "-outfile=s"    => \$outfile,
	 "-loadModules=s{,}" => \@loadModules,
	 "-I_have_loaded_all_modules" => \$I_have_loaded_all_modules,
	 "-max_jobs=s" => \$max_jobs,

	 "-help"             => \$help,
	 "-debug"            => \$debug
);

my $warn = '';
my $error = '';

unless ( defined $cmd) {
	$error .= "the cmd line switch -cmd is undefined!\n";
}
unless ( defined $options[0]) {
	$warn .= "the cmd line switch -options is undefined!\n";
}
unless ( defined $outfile) {
	$error .= "the cmd line switch -outfile is undefined!\n";
}
if ( defined $loadModules[0]){
	## Cool process that later.	
}else {
unless ( $I_have_loaded_all_modules ){
	$error .= "Please make sure you have loaded all required modules and try again! (-I_have_loaded_all_modules missing)\n"
}
}
unless ( $max_jobs ){
	$max_jobs = 40;
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

$task_description .= 'perl '.$plugin_path .'/runCommand.pl';
$task_description .= " -cmd '$cmd'" if (defined $cmd);
$task_description .= ' -options "'.join( '" "', @options ).'"' if ( defined $options[0]);
$task_description .= " -outfile '$outfile'" if (defined $outfile);
$task_description .= " -max_jobs $max_jobs";

my $fm = root->filemap( $outfile);
unless ( -d $fm->{'path'} ) {
	system( "mkdir -p $fm->{'path'}" );
}
for ( my $i = 0 ; $i < @options ; $i += 2 ) {
	$options[ $i + 1 ] =~ s/\n/ /g;
	$options->{ $options[$i] } = $options[ $i + 1 ];
}
open ( LOG , ">$outfile.log") or die $!;
print LOG $task_description."\n";
close ( LOG );


## Do whatever you want!
$options->{'A'} ||= 'lu2016-2-7';
$options->{'n'} ||= 10;
$options->{'N'} ||= 1;
$options->{'t'} ||= '02:00:00';
$options->{'proc'} ||= $options->{'n'}*$options->{'N'};
$options->{'p'} ||= $options->{'proc'};
$options->{'debug'} = $debug;
$options->{'max_jobs'} ||= $max_jobs;

my $SLURM =stefans_libs::SLURM->new( $options );

if ( defined $loadModules[0]) {
	## create a copy - might be better...
	$SLURM->{'SLURM_modules'} = [@loadModules];
}
$SLURM->{'debug'} = $debug;

if ( $cmd =~ m/\$SNIC_TMP/ ) {
	## Ohoh - the output will go into the node specific $SNIC_TMP folder and has to be copied into the right outpath afterwards!
	$cmd .= "\ncp -R \$SNIC_TMP/* $fm->{'path'}\n";
}

$SLURM->run( $cmd, $fm );

print "Done\n";