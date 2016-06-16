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
       -cmd       :<please add some info!>
       -options     :<please add some info!> you can specify more entries to that
                         format: key_1 value_1 key_2 value_2 ... key_n value_n
       -outfile       :<please add some info!>


       -help           :print this help
       -debug          :verbose output
   
=head1 DESCRIPTION

  Gets one command line and runs it on the SLURM backend.

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


my ( $help, $debug, $database, $cmd, $options, @options, $outfile);

Getopt::Long::GetOptions(
	 "-cmd=s"    => \$cmd,
       "-options=s{,}"    => \@options,
	 "-outfile=s"    => \$outfile,

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

$options->{'n'} ||= 10;
$options->{'N'} ||= 1;
$options->{'t'} ||= '02:00:00';
$options->{'proc'} ||= $options->{'n'}*$options->{'N'};
$options->{'p'} ||= $options->{'proc'};
$options->{'debug'} = $debug;

my $SLURM =stefans_libs::SLURM->new( $options );

$SLURM->run( $cmd, $fm );

print "Done\n";