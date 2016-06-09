#! /usr/bin/perl -w

=head1 LICENCE

  Copyright (C) 2016-06-09 Stefan Lang

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

    hisat_run_aurora.pl
       -files     :<please add some info!> you can specify more entries to that
       -options     :<please add some info!> you can specify more entries to that
                         format: key_1 value_1 key_2 value_2 ... key_n value_n
                         
       -genome       :<please add some info!>
       -coverage       :<please add some info!>
       -paired       :<please add some info!>
       -bigwigTracks       :<please add some info!>
       -sra       :<please add some info!>


       -help           :print this help
       -debug          :verbose output
   
=head1 DESCRIPTION

  wrapper script to run hisat on NGS expression data fastq and sra files.

  To get further help use 'hisat_run_aurora.pl -help' at the comman line.

=cut

use Getopt::Long;
use Pod::Usage;

use strict;
use warnings;

use stefans_libs::root;
use stefans_libs::SLURM;
use stefans_libs::scripts::BAM;

use FindBin;
my $plugin_path = "$FindBin::Bin";

use Cwd;
my $dir = getcwd . "/";

my $VERSION = 'v1.0';


my ( $help, $debug, $database, @files, $options, @options, $genome, $coverage, $paired, $bigwigTracks, $sra);

Getopt::Long::GetOptions(
       "-files=s{,}"    => \@files,
       "-options=s{,}"    => \@options,
	 "-genome=s"    => \$genome,
	 "-coverage=s"    => \$coverage,
	 "-paired"    => \$paired,
	 "-bigwigTracks=s"    => \$bigwigTracks,
	 "-sra"    => \$sra,

	 "-help"             => \$help,
	 "-debug"            => \$debug
);

my $warn = '';
my $error = '';

unless ( defined $files[0]) {
	$error .= "the cmd line switch -files is undefined!\n";
}
unless ( defined $options[0]) {
	$error .= "the cmd line switch -options is undefined!\n";
}
unless ( defined $genome) {
	$error .= "the cmd line switch -genome is undefined!\n";
}
unless ( -f $coverage) {
	$error .= "the cmd line switch -coverage is undefined ($coverage)!\n";
}
unless ( defined $bigwigTracks) {
	$error .= "the cmd line switch -bigwigTracks is undefined!\n";
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

$task_description .= 'perl '.$plugin_path .'/hisat_run_aurora.pl';
$task_description .= ' -files "'.join( '" "', @files ).'"' if ( defined $files[0]);
$task_description .= ' -options "'.join( '" "', @options ).'"' if ( defined $options[0]);
$task_description .= " -genome '$genome'" if (defined $genome);
$task_description .= " -coverage '$coverage'" if (defined $coverage);
$task_description .= " -paired" if ($paired);
$task_description .= " -bigwigTracks '$bigwigTracks'" if (defined $bigwigTracks);
$task_description .= " -sra" if ($sra);


for ( my $i = 0 ; $i < @options ; $i += 2 ) {
	$options[ $i + 1 ] =~ s/\n/ /g;
	$options->{ $options[$i] } = $options[ $i + 1 ];
}

## Do whatever you want!

my ( $cmd, $fm, @big_wig_urls, $tmp, $this_outfile );

$options->{'n'} ||= 10;
$options->{'N'} ||= 1;
$options->{'t'} ||= '02:00:00';

my $SLURM = stefans_libs::SLURM->new($options);
$SLURM->{'debug'} = 1 if ($debug);

## kick all SLURM options that should not be used for the hisat
foreach (qw(n N t mem)) {
	delete( $options->{$_} );
}

my $BAM = stefans_libs::scripts::BAM->new($options);

@files = map { $_ =~ s/^\.\///; $_ } @files;

while ( scalar(@files) ){
	$fm  = root->filemap( $files[0] ); ## the files will be depleted by the create_call function!
	$cmd = &chk_cmd(&create_call());
	$cmd .= &chk_cmd($BAM->convert_sam_2_sorted_bam($this_outfile));
	$cmd .= &chk_cmd($BAM->convert_sorted_bam_2_bedGraph($this_outfile,$coverage));
	$cmd .= &chk_cmd($BAM->convert_bedGraph_2_bigwig($this_outfile,$coverage));
	$SLURM->run( $cmd, $fm, $this_outfile);
}

if ( @{$BAM->{'big_wig_urls'}} > 0 ) {
	open( OUT, ">$bigwigTracks" )
	  or die "I could not create the bigwig outfile '$bigwigTracks'\n$!\n";
	print OUT join( "\n", @{$BAM->{'big_wig_urls'}} );
	close(OUT);
	print join( "\n", @{$BAM->{'big_wig_urls'}} ) . "\n\n";
	open( LOG, ">$bigwigTracks.log" )
	  or die "I could not open the log file '$bigwigTracks.log'\n$!\n";
	print LOG $task_description . "\n";
	close(LOG);
}



print "Done\n";

sub create_call {
	return &create_paired_call() if ( $paired );
	return &create_sra_call() if ( $sra );
    my $file = shift(@files);
    my $fm   = root->filemap($file);
    mkdir ( "$fm->{'path'}/HISAT_OUT/" ) unless ( -d "$fm->{'path'}/HISAT_OUT/");
    $fm->{'path'} = "" if ( $fm->{'path'} eq "./" );
    unless ( $fm->{'path'} =~ m/^\// ) { $fm->{'path'} = $dir . $fm->{'path'}; }
    my $s =
"hisat -x $genome -U $fm->{'total'} --threads $options->{proc} > $fm->{'path'}/HISAT_OUT/$fm->{'filename_core'}_hisat.sam\n";
    return $s, "$fm->{'path'}/HISAT_OUT/$fm->{'filename_core'}_hisat.sam";
}

sub create_sra_call {
    my $file = shift( @files);
    my $fm   = root->filemap($file);
    mkdir ( "$fm->{'path'}/HISAT_OUT/" ) unless ( -d "$fm->{'path'}/HISAT_OUT/");
    $fm->{'path'} = "" if ( $fm->{'path'} eq "./" );
    unless ( $fm->{'path'} =~ m/^\// ) { $fm->{'path'} = $dir . $fm->{'path'}; }
    my $s =
"hisat -x $genome --sra-acc $fm->{'total'} --threads $options->{proc} > $fm->{'path'}/HISAT_OUT/$fm->{'filename_core'}_hisat.sam\n";
    return $s, "$fm->{'path'}/HISAT_OUT/$fm->{'filename_core'}_hisat.sam";
}

sub create_paired_call {
    my $file = shift(@files);
    my $pair = shift(@files);
    my $fm   = root->filemap($file);
    my $fm2  = root->filemap($pair);
    $fm->{'path'} = "" if ( $fm->{'path'} eq "./" );
    $fm2->{'path'} = "" if ( $fm2->{'path'} eq "./" );
    $fm->{'path'} .= '/' unless ( $fm->{'path'} =~ m/\/$/ );
    $fm2->{'path'} .= '/' unless ( $fm2->{'path'} =~ m/\/$/ );
    mkdir ( "$fm->{'path'}/HISAT_OUT/" ) unless ( -d "$fm->{'path'}/HISAT_OUT/");
    unless ( $fm->{'path'} =~ m/^\// ) { $fm->{'path'} = $dir . $fm->{'path'}; }
    my $s =
"hisat -x $genome -1 $fm->{'total'} -2 $fm2->{'total'} --threads $options->{proc} > $fm->{'path'}/HISAT_OUT/$fm->{'filename_core'}_hisat.sam\n";
    return $s, "$fm->{'path'}/HISAT_OUT/$fm->{'filename_core'}_hisat.sam";
}

sub chk_cmd {
	my ( $cmd, $ofile ) = @_;
	$this_outfile = $ofile;
	return join( "\n",
		map { $SLURM->check_4_outfile( $_, $ofile ) } split( /\n/, $cmd ) )
	  . "\n";
}

