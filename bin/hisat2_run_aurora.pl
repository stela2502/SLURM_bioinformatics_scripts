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

    hisat2_run_aurora.pl
       -files     :the input fastq or sra files
       -outpath   :the outpath for the mapped files 
                   default to subpath of infile named HISAT2_OUT 
       -sra       :the files are in sra format, not fastq
       -max_jobs  :Define how many different jobs might be sublitted at a time.
                   default 40
       -local     :run the alignements one by one on the local comuter
                   do not use the cluster (not recommended!)
       -options   :format: key_1 value_1 key_2 value_2 ... key_n value_n
       		n    :number of cores per node (default = 10 )
            N    :number of nodes (default =1)
            t    :time untill the script is stopped (default =02:00:0 (2h))
            proc :hisat 2 number of threads (default = n*N )
                   
       -genome       :the hisat2 genome information path
       -coverage     :the chromome length file
       -paired       :analyse paried fasta files 
                      every first file == read 1 every second == read2
       -bigwigTracks :the bigwig tracks file you can upload to the USCS genome browser


       -help           :print this help
       -debug          :verbose output
   
=head1 DESCRIPTION

  wrapper script to run hisat2 on NGS expression data fastq and sra files.

  To get further help use 'hisat2_run_aurora.pl -help' at the comman line.

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


my ( $help, $debug, $database, @files, $options, @options, $genome, $coverage, $paired, $bigwigTracks, $sra, $outpath, $max_jobs);

Getopt::Long::GetOptions(
       "-files=s{,}"    => \@files,
       "-options=s{,}"    => \@options,
	 "-genome=s"    => \$genome,
	 "-coverage=s"    => \$coverage,
	 "-outpath=s" => \$outpath,
	 "-paired"    => \$paired,
	 "-bigwigTracks=s"    => \$bigwigTracks,
	 "-sra"    => \$sra,
	 "-max_jobs=s" => \$max_jobs,

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
unless ( $max_jobs ) {
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


for ( my $i = 0 ; $i < @options ; $i += 2 ) {
	$options[ $i + 1 ] =~ s/\n/ /g;
	$options->{ $options[$i] } = $options[ $i + 1 ];
}
### initialize default options:

$options->{'n'} ||= 10;
$options->{'N'} ||= 1;
$options->{'t'} ||= '02:00:00';
$options->{'proc'} ||= $options->{'n'}*$options->{'N'};
$options->{'p'} ||= $options->{'proc'};

###
 
my ( $task_description);

$task_description .= 'perl '.$plugin_path .'/hisat2_run_aurora.pl';
$task_description .= ' -files "'.join( '" "', @files ).'"' if ( defined $files[0]);
$task_description .= ' -options "'.join( '" "', @options ).'"' if ( defined $options[0]);
$task_description .= " -genome '$genome'" if (defined $genome);
$task_description .= " -coverage '$coverage'" if (defined $coverage);
$task_description .= " -paired" if ($paired);
$task_description .= " -bigwigTracks '$bigwigTracks'" if (defined $bigwigTracks);
$task_description .= " -sra" if ($sra);
$task_description .= " -max_jobs $max_jobs";
$task_description .= " -debug" if ($debug);


## Do whatever you want!

my ( $cmd, $fm, @big_wig_urls, $tmp, $this_outfile );

my $SLURM = stefans_libs::SLURM->new($options);
$SLURM->{'debug'} = 1 if ($debug);

## kick all SLURM options that should not be used for the hisat run
foreach (qw(n N t mem)) {
	delete( $options->{$_} );
}

$fm = root->filemap( $files[0] );
open( SC, ">$fm->{'path'}/InitializeSLURMenv.sh" )
  or die "I could not create the SLURM init script\n$!\n";

foreach (
	'icc/2016.1.150-GCC-4.9.3-2.25 impi/5.1.2.150 HISAT2/2.0.4',
	stefans_libs::scripts::BAM->SLURUM_load(),
  )
{
	print SC "module load $_\n";
}
close(SC);
unless ($debug) {
	system("bash $fm->{'path'}/InitializeSLURMenv.sh");
}

my $BAM = stefans_libs::scripts::BAM->new($options);

@files = map { $_ =~ s/^\.\///; $_ } @files;
my $submitted = 0;

while ( scalar(@files) ){
	$fm  = root->filemap( $files[0] ); ## the files will be depleted by the create_call function!
	$cmd = &chk_cmd(&create_call());
	$cmd .= &chk_cmd($BAM->convert_sam_2_sorted_bam($this_outfile));
	$cmd .= &chk_cmd($BAM->convert_sorted_bam_2_bedGraph($this_outfile,$coverage));
#	$cmd .= &chk_cmd($BAM->convert_bedGraph_2_bigwig($this_outfile,$coverage)); # brken in the scripts
	$tmp = $SLURM->run( $cmd, $fm, $this_outfile);
	$submitted ++ if ( $tmp == 1 );
	if ( $submitted >= $max_jobs ) {
		$submitted -= 50;
		&wait_for_last_finished( $this_outfile );
	}
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


sub wait_for_last_finished {
	my ( $fname ) = @_;
	my $wait;
	while( $wait = &in_pipeline() > 4 ) {
		warn "waiting for $wait processes\n";
		sleep( 50);
	}
}

sub in_pipeline {
	open ( IN ,"whoami |" ) or die $!;
	my @IN = <IN>;
	close ( IN );
	my $name = $IN[0];
	chomp($name);
	open(IN, "squeue -u $name |" ) or die $!;
	@IN = <IN>;
	close ( IN );
	## pending
	return scalar(grep ( /\sPD\s/, @IN));
	## all
	return scalar(@IN)-1;
}


sub create_call {
	return &create_paired_call() if ( $paired );
	return &create_sra_call() if ( $sra );
    my $file = shift(@files);
    my $fm   = root->filemap($file);
	my $p = $outpath;
	$p ||= "$fm->{'path'}/HISAT2_OUT/";
    mkdir ( "$p/" ) unless ( -d "$p/");
    $fm->{'path'} = "" if ( $fm->{'path'} eq "./" );
    unless ( $fm->{'path'} =~ m/^\// ) { $fm->{'path'} = $dir . $fm->{'path'}; }
    my $s =
"hisat2 -x $genome -U $fm->{'total'} --threads $options->{proc} --add-chrname > $p/$fm->{'filename_core'}_hisat.sam\n";
    return $s, "$p/$fm->{'filename_core'}_hisat.sam", "$p/$fm->{'filename_core'}_hisat.sorted.bam";
}

sub create_sra_call {
    my $file = shift( @files);
    my $fm   = root->filemap($file);
	my $p = $outpath;
	$p ||= "$fm->{'path'}/HISAT2_OUT/";
    mkdir ( "$p/" ) unless ( -d "$p/");
    $fm->{'path'} = "" if ( $fm->{'path'} eq "./" );
    unless ( $fm->{'path'} =~ m/^\// ) { $fm->{'path'} = $dir . $fm->{'path'}; }
    my $s =
"hisat2 -x $genome --sra-acc $fm->{'total'} --threads $options->{proc} --add-chrname > $p/$fm->{'filename_core'}_hisat.sam\n";
    return $s, "$p/$fm->{'filename_core'}_hisat.sam", "$p/$fm->{'filename_core'}_hisat.sorted.bam";
}

sub create_paired_call {
    my $file = shift(@files);
    my $pair = shift(@files);
    my $fm   = root->filemap($file);
    my $fm2  = root->filemap($pair);
	my $p = $outpath;
	$p ||= "$fm->{'path'}/HISAT2_OUT/";
    $fm->{'path'} = "" if ( $fm->{'path'} eq "./" );
    $fm2->{'path'} = "" if ( $fm2->{'path'} eq "./" );
    $fm->{'path'} .= '/' unless ( $fm->{'path'} =~ m/\/$/ );
    $fm2->{'path'} .= '/' unless ( $fm2->{'path'} =~ m/\/$/ );
    mkdir ( "$p/" ) unless ( -d "$p/");
    unless ( $fm->{'path'} =~ m/^\// ) { $fm->{'path'} = $dir . $fm->{'path'}; }
    my $s =
"hisat2 -x $genome -1 $fm->{'total'} -2 $fm2->{'total'} --threads $options->{proc} --add-chrname > $p/$fm->{'filename_core'}_hisat.sam\n";
    return $s, "$p/$fm->{'filename_core'}_hisat.sam", "$p/$fm->{'filename_core'}_hisat.sorted.bam";
}

sub chk_cmd {
	my ( $cmd, @ofiles ) = @_;
	$this_outfile = $ofiles[0];
	return join( "\n",
		map { $SLURM->check_4_outfile( $_, @ofiles ) } split( /\n/, $cmd ) )
	  . "\n";
}

