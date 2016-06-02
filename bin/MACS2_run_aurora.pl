#! /usr/bin/perl -w

use Getopt::Long;
use Pod::Usage;
use stefans_libs::root;
use stefans_libs::SLURM;
use stefans_libs::flexible_data_structures::data_table;

=head1 LICENCE

  Copyright (C) 2016-06-01 Stefan Lang

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

    MACS2_run_aurora.pl
       -files     :the sorted.bam files or sorted.rmdup.bam files you have
       -cfiles    :the control sorted.bam files or sorted.rmdup.bam files you have
       -options   :additional options to the MACS2 in the format
                   key_1 value_1 key_2 value_2 ... key_n value_n
                   e.g.: g hs f BAM q 0.01
       -broad     :whether to call broad peaks or not
       -bigWig    :produce bigwig files?
       -normdup   :NOT drop pcr duplicates?


       -help   :print this help
       -debug  :verbose output
   
=head1 DESCRIPTION

  This does take a mass of sorted.bam files, checks that they are either rmduped or does that and calls MACS2 afterwards.

  To get further help use 'MACS2_run_aurora.pl -help' at the comman line.

=cut


use strict;
use warnings;

use FindBin;
my $plugin_path = "$FindBin::Bin";

my $VERSION = 'v1.0';


my ( $help, $debug, $database, @files, @cfiles,$options, @options, $broad, $bigWig, $normdup);

Getopt::Long::GetOptions(
     "-files=s{,}"    => \@files,
     "-cfiles=s{,}"    => \@cfiles,
     "-options=s{,}"    => \@options,
	 "-broad"    => \$broad,
	 "-bigWig"   => \$bigWig,
	 "-normdup"    => \$normdup,

	 "-help"             => \$help,
	 "-debug"            => \$debug
);

my $warn = '';
my $error = '';

unless ( -f $files[0]) {
	$error .= "the cmd line switch -files is undefined!\n";
}
unless ( -f $cfiles[0]) {
	my $d = data_table->new ( {'filename' => $files[0] } );
	if ( defined $d->Header_Position( 'files') ) {
		@files = @{$d->GetAsArray('files')};
		unless ( -f $files[0]) {
			$error .= "the cmd line switch -files is undefined!\n";
		}
	}
	if ( defined $d->Header_Position( 'cfiles') ){
		@cfiles = @{$d->GetAsArray('cfiles')};
		unless ( -f $cfiles[0]) {
			$error .= "the cmd line switch -cfiles is undefined!\n";
		}
	}else {
		$error .= "the cmd line switch -cfiles is undefined!\n";
	}
}
unless ( defined $options[0]) {
	$error .= "the cmd line switch -options is undefined!\n";
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

$task_description .= 'perl '.$plugin_path .'/MACS2_run_aurora.pl';
$task_description .= ' -files "'.join( '" "', @files ).'"' if ( defined $files[0]);
$task_description .= ' -options "'.join( '" "', @options ).'"' if ( defined $options[0]);
$task_description .= " -broad" if ($broad);
$task_description .= " -bigWig" if ($bigWig);
$task_description .= " -normdup" if ($normdup);

for ( my $i = 0 ; $i < @options ; $i += 2 ) {
	$options[ $i + 1 ] =~ s/\n/ /g;
	$options->{ $options[$i] } = $options[ $i + 1 ];
}

### set SLURM default options
## n N t
$options ->{'n'} ||=1;
$options ->{'N'} ||=1;
$options ->{'t'} ||='02:00:00';


## Do whatever you want!
my ( $SLURM, $cmd, $tmp, $outfile, $file, $cfile, $fm, $cfm );
$SLURM = stefans_libs::SLURM->new( $options );

## kick all SLURM options that should not be used for the MACS2
foreach (qw(n N t mem)) {
	delete( $options->{$_} );
}

$SLURM->{'debug'} = 1 if ( $debug );
stefans_libs::SLURM::check( $options, qw( f g q ) );

for ( my $i = 0; $i <@files; $i ++ ) {
	$file = $files[$i];
	Carp::confess ( "bam file '$file' does not exist\n$!\n") unless ( -f $file);
	$cfile= $cfiles[$i];
	Carp::confess ( "control bam file '$file' does not exist\n$!\n") unless ( -f $cfile);
	( $tmp, $fm ) = &rmdup( root->filemap( $file ) );
	$cmd = $tmp;
	print "The file rmdup command: $tmp\n";
	( $tmp, $cfm ) = &rmdup( root->filemap( $cfile ) );
	$cmd .= $tmp;
	print "The cfile rmdup command: $tmp\n";
	$fm->{'path'}.="/MACS2_out";
	mkdir ( $fm->{'path'} ) unless ( -d $fm->{'path'} );
	$outfile = "$fm->{'path'}/$fm->{'filename_core'}";
	$cmd .= $SLURM->check_4_outfile( "MACS2 -t $fm->{'total'} -c $cfm->{'total'} ", $outfile."_peaks.bed");
	foreach my $key ( keys %$options ) {
		$cmd .= " -$key $options->{$key}";
	}
	if ( $bigWig ){
		$cmd .= " -B";
	}
	if ( $broad ){
		$cmd .= " --broad";
	}
	$cmd .= " -n $outfile\n";
	$SLURM -> run ( $cmd, $fm );
}


sub rmdup {
	my $fm = shift;
	return ('' , $fm) if ( $normdup );
	return ('' , $fm) unless ( $fm->{'filename_ext'} eq "bam" );
	my $outfile = "$fm->{'path'}/$fm->{'filename_core'}.rmdup.bam";
	$cmd = $SLURM->check_4_outfile( "samtools rmdup $fm->{'total'} $outfile\n", $outfile);
	$fm = root->filemap( $outfile );
	return ( $cmd, $fm );
}