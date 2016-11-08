#! /usr/bin/perl -w

=head1 LICENCE

  Copyright (C) 2016-10-28 Stefan Lang

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

    iCLIP_postprocessing.pl
       -infiles     :the mapped bam files
       -outpath     :the outpath for the peaks and binding site models
       -genome      :the path to the fasta genome files
       -options     :Here I only need your SNIC account e.g.
                       -A 'snic2016-4-13'

       -help           :print this help
       -debug          :verbose output
   
=head1 DESCRIPTION

  A tool to postprocess iCLIP data using Piranja and Zargos on a slurm instance.

  To get further help use 'iCLIP_postprocessing.pl -help' at the comman line.

=cut

use Getopt::Long;
use Pod::Usage;

use stefans_libs::SLURM;

use strict;
use warnings;

use FindBin;
my $plugin_path = "$FindBin::Bin";

my $VERSION = 'v1.0';


my ( $help, $debug, $database, @infiles, $outpath, $options, @options, $genome);

Getopt::Long::GetOptions(
       "-infiles=s{,}"    => \@infiles,
	 "-outpath=s"    => \$outpath,
       "-options=s{,}"    => \@options,
       "-genome=s" => \$genome,

	 "-help"             => \$help,
	 "-debug"            => \$debug
);

my $warn = '';
my $error = '';

unless ( -f $infiles[0]) {
	$error .= "the cmd line switch -infiles is undefined or the first file does not exist!\n";
}
unless ( -d $genome) {
	$error .= "the cmd line switch -genome is undefined or the path does not exist!\n";
}
unless ( defined $outpath) {
	$error .= "the cmd line switch -outpath is undefined!\n";
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

#$options->{'n'} ||= 10;

###


my ( $task_description);

$task_description .= 'perl '.$plugin_path .'/iCLIP_postprocessing.pl';
$task_description .= ' -infiles "'.join( '" "', @infiles ).'"' if ( defined $infiles[0]);
$task_description .= " -outpath '$outpath'" if (defined $outpath);
$task_description .= ' -options "'.join( '" "', @options ).'"' if ( defined $options[0]);


for ( my $i = 0 ; $i < @options ; $i += 2 ) {
	$options[ $i + 1 ] =~ s/\n/ /g;
	$options->{ $options[$i] } = $options[ $i + 1 ];
}
###### default options ########
#$options->{'something'} ||= 'default value';
##############################
mkdir( $outpath ) unless ( -d $outpath );
open ( LOG , ">$outpath/".$$."_iCLIP_postprocessing.pl.log") or die $!;
print LOG $task_description."\n";
close ( LOG );


## Do whatever you want!

## single processor scripts - take forever to run
$options->{'n'} ||= 1;
$options->{'N'} ||= 1;
$options->{'t'} ||= '48:00:00';
$options->{'proc'} ||= $options->{'n'}*$options->{'N'};
$options->{'p'} ||= $options->{'proc'};

# from what I have understood you need to run first Pirnaha and then Zafia, 
# both from the http://smithlabresearch.org/software/


my $SLURM = stefans_libs::SLURM->new($options);
$SLURM->{'debug'} = 1 if ($debug);

foreach (qw(n N t mem A)) {
	delete( $options->{$_} );
}

my $fm = root->filemap( $infiles[0] );
open( SC, ">$fm->{'path'}/InitializeSLURMenv.sh" )
  or die "I could not create the SLURM init script\n$!\n";

$options->{'modules'} ||= 'GCCcore/4.9.3 GCC/4.9.3-2.25 GSL/2.1 ucsc-tools/R2016a GSL/2.1 BEDTools/2.25.0';

print SC "module load $options->{'modules'}\n";

close(SC);
unless ($debug) {
	system("bash $fm->{'path'}/InitializeSLURMenv.sh");
}

my $cmd;
foreach $fm ( map { root->filemap($_) } @infiles ) {
	#bamToBed
	$cmd = $SLURM-> check_4_outfile ( 
		"bamToBed -i $fm->{'total'} > $fm->{'path'}/$fm->{'filename_core'}.bed", 
		"$fm->{'path'}/$fm->{'filename_core'}.bed" 
	)."\n";
	#Piranha peak caller
	$cmd .= $SLURM-> check_4_outfile ( 
		"Piranha -s -p 0.01  $fm->{'path'}/$fm->{'filename_core'}.bed -dist Poisson > $outpath/$fm->{'filename_core'}.peaks",
		"$outpath/$fm->{'filename_core'}.peaks"
	)."\n";
	#thermo estimation of RNA structure close to the peak
	$cmd .= $SLURM-> check_4_outfile ( 
		"thermo -c $genome -o $outpath/$fm->{'filename_core'}.thermo $outpath/$fm->{'filename_core'}.peaks", 
		"$outpath/$fm->{'filename_core'}.thermo" 
	)."\n";
	#zagros identification of putative binding sites
	$cmd .= $SLURM-> check_4_outfile ( 
		"zagros -c $genome -t $outpath/$fm->{'filename_core'}.thermo $outpath/$fm->{'filename_core'}.peaks > $outpath/$fm->{'filename_core'}.zagros", 
		"$outpath/$fm->{'filename_core'}.zagros" 
	)."\n";
	
	$SLURM->run( $cmd, $fm);
}



