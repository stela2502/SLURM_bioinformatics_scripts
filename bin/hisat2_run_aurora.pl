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
     partitition :explicitely call for a partitition in the SLURM environment (slurm p option)
                  should not be necessary!
                  
       -fast_tmp :the nodes should have a fast local tmp folder that could be used for 
                  the intermediate files (default '$SNIC_TMP')

       -mapper_options :specific options for the hisat2 as one string
                        e.g. -mapper_options '--max-seeds 500'

       -genome       :the hisat2 genome information path
       -coverage     :the chromome length file
       -paired       :analyse paried fasta files
                      every first file == read 1 every second == read2

       -dropDuplicates :use picard to drop duplicates

       -bigwigTracks :the bigwig tracks file you can upload to the USCS genome browser

ls

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

my (
	$help,     $debug,   $database,       @files,
	$options,  @options, $dropDuplicates, $genome,
	$coverage, $paired,  $fast_tmp,       $bigwigTracks,
	$sra,      $outpath, $max_jobs,       $mapper_options,
	$local,
);

Getopt::Long::GetOptions(
	"-files=s{,}"       => \@files,
	"-options=s{,}"     => \@options,
	"-genome=s"         => \$genome,
	"-coverage=s"       => \$coverage,
	"-outpath=s"        => \$outpath,
	"-paired"           => \$paired,
	"-bigwigTracks=s"   => \$bigwigTracks,
	"-sra"              => \$sra,
	"-max_jobs=s"       => \$max_jobs,
	"-dropDuplicates"   => \$dropDuplicates,
	"-mapper_options=s" => \$mapper_options,
	"-fast_tmp=s"       => \$fast_tmp,
	"-local"          => \$local,

	"-help"  => \$help,
	"-debug" => \$debug
);

my $warn  = '';
my $error = '';

unless ( defined $files[0] ) {
	$error .= "the cmd line switch -files is undefined!\n";
}
unless ( defined $options[0] ) {
	$error .= "the cmd line switch -options is undefined!\n";
}
unless ( defined $genome ) {
	$error .= "the cmd line switch -genome is undefined!\n";
}
unless ( -f $coverage ) {
	$error .= "the cmd line switch -coverage is undefined ($coverage)!\n";
}
unless ( defined $bigwigTracks ) {
	$error .= "the cmd line switch -bigwigTracks is undefined!\n";
}
unless ( defined $mapper_options ) {
	$mapper_options = '';
}
unless ($max_jobs) {
	$max_jobs = 40;
}
unless ($fast_tmp) {
	$fast_tmp = '$SNIC_TMP';
}

if ($help) {
	print helpString();
	exit;
}

if ( $error =~ m/\w/ ) {
	helpString($error);
	exit;
}

sub helpString {
	my $errorMessage = shift;
	$errorMessage = ' ' unless ( defined $errorMessage );
	print "$errorMessage.\n";
	pod2usage( q(-verbose) => 1 );
}

for ( my $i = 0 ; $i < @options ; $i += 2 ) {
	$options[ $i + 1 ] =~ s/\n/ /g;
	$options->{ $options[$i] } = $options[ $i + 1 ];
}
### initialize default options:

$options->{'n'}    ||= 10;
$options->{'N'}    ||= 1;
$options->{'t'}    ||= '02:00:00';
$options->{'proc'} ||= $options->{'n'} * $options->{'N'};
$options->{'p'}    ||= $options->{'proc'};

###

my ( $task_description, $mainOutpath );

$task_description .= 'perl ' . $plugin_path . '/hisat2_run_aurora.pl';
$task_description .= ' -files "' . join( '" "', @files ) . '"'
  if ( defined $files[0] );
$task_description .= ' -local ' if ($local); 
$task_description .= ' -options "' . join( '" "', @options ) . '"'
  if ( defined $options[0] );
$task_description .= " -genome '$genome'"     if ( defined $genome );
$task_description .= " -coverage '$coverage'" if ( defined $coverage );
$task_description .= " -paired"               if ($paired);
$task_description .= " -bigwigTracks '$bigwigTracks'"
  if ( defined $bigwigTracks );
$task_description .= " -sra"            if ($sra);
$task_description .= " -dropDuplicates" if ($dropDuplicates);
$task_description .= " -max_jobs $max_jobs";
$task_description .= " -mapper_options '$mapper_options'"
  if ( $mapper_options =~ m/\w/ );
$task_description .= " -fast_tmp '$fast_tmp'" if ( $fast_tmp =~ m/\w/ );
$task_description .= " -debug" if ($debug);

## Do whatever you want!
my $fm = root->filemap($bigwigTracks);

unless ( -d $fm->{'path'} ) {
	mkdir( $fm->{'path'} );
}
open( LOG, ">$bigwigTracks.log" )
  or die "I could not create the log file '$bigwigTracks.log'\n$!\n";
print LOG $task_description;
close(LOG);

my ( @cmd, @big_wig_urls, $tmp, $this_outfile );

my $SLURM = stefans_libs::SLURM->new($options);
$SLURM->{'debug'} = 1 if ($debug);
$SLURM->{'local'} = $local;

## kick all SLURM options that should not be used for the hisat run
foreach (qw(n N t mem)) {
	delete( $options->{$_} );
}

$fm = root->filemap( $files[0] );

unless ( $local ){
$SLURM->{'SLURM_modules'} = [
	'GCC/4.9.3-2.25', 'OpenMPI/1.10.2',
	'icc/2016.1.150-GCC-4.9.3-2.25', 'impi/5.1.2.150', 'SAMtools/1.3.1',
	'HISAT2/2.0.4',
	'BEDTools/2.25.0', 'Java/1.8.0_72', 'picard/2.8.2', 'ucsc-tools/R2016a',

	#	stefans_libs::scripts::BAM->SLURUM_load(),
];
}

$tmp                    = $SLURM->define_Subscript();
$tmp->{'purge'}         = 1;
unless ( $local ){
$tmp->{'SLURM_modules'} = [
	'GCC/4.9.3-2.25', 'OpenMPI/1.10.2',
	'icc/2016.1.150-GCC-4.9.3-2.25', 'impi/5.1.2.150', 'BEDTools/2.25.0',
	'Java/1.8.0_72', 'picard/2.8.2', 'ucsc-tools/R2016a',
];
}
my $BAM = stefans_libs::scripts::BAM->new($options);

@files = map { $_ =~ s/^\.\///; $_ } @files;
my $submitted = 0;

while ( scalar(@files) ) {
	$fm = root->filemap( $files[0] )
	  ;    ## the files will be depleted by the create_call function!
	$cmd[0] = &chk_cmd( &create_call() );
	$cmd[0] .=
	  &chk_cmd( $BAM->convert_sam_2_sorted_bam( $this_outfile, $mainOutpath ) );
	$cmd[1] .= &chk_cmd( create_picard_call($this_outfile) ) if ($dropDuplicates);
	$cmd[1] .=
	  &chk_cmd(
		$BAM->convert_sorted_bam_2_bedGraph( $this_outfile, $coverage ) );

	$cmd[1] .=
	  &chk_cmd( $BAM->convert_bedGraph_2_bigwig( $this_outfile, $coverage ) );
	$tmp = $SLURM->run( \@cmd, $fm, $this_outfile );
	$submitted++ if ( $tmp == 1 );
	if ( $submitted >= $max_jobs ) {
		$submitted -= 50;
		$SLURM->wait_for_last_finished($this_outfile);
	}
}

if ( @{ $BAM->{'big_wig_urls'} } > 0 ) {
	open( OUT, ">$bigwigTracks" )
	  or die "I could not create the bigwig outfile '$bigwigTracks'\n$!\n";
	print OUT join( "\n", @{ $BAM->{'big_wig_urls'} } );
	close(OUT);
	print join( "\n", @{ $BAM->{'big_wig_urls'} } ) . "\n\n";
	open( LOG, ">$bigwigTracks.log" )
	  or die "I could not open the log file '$bigwigTracks.log'\n$!\n";
	print LOG $task_description . "\n";
	close(LOG);
}

print "Done\n";

sub create_picard_call {
	my ($file) = @_;
	my $fm     = root->filemap($file);
	my $p      = $outpath;
	$p ||= "$fm->{'path'}/";
	my $cmd = "picard";
	if ( $local) {
		$cmd = "picard-tools";
	}
	my $s =
	    "$cmd MarkDuplicates INPUT="
	  . $fm->{'total'}
	  . " OUTPUT="
	  . "$p$fm->{'filename_core'}_picard_deduplicated.bam"
	  . " REMOVE_DUPLICATES=TRUE M="
	  . "$p$fm->{'filename_core'}_picard_metrix.txt";
	return $s, "$p$fm->{'filename_core'}_picard_deduplicated.bam";
}

sub create_call {
	return &create_paired_call() if ($paired);
	return &create_sra_call()    if ($sra);
	my $file = shift(@files);
	my $fm   = root->filemap($file);
	my $p    = $outpath;
	$p ||= "$fm->{'path'}/HISAT2_OUT";
	$mainOutpath = $p;
	mkdir("$p/") unless ( -d "$p/" );
	$fm->{'path'} = "" if ( $fm->{'path'} eq "./" );
	unless ( $fm->{'path'} =~ m/^\// ) { $fm->{'path'} = $dir . $fm->{'path'}; }
	my $s =
"hisat2 $mapper_options -x $genome -U $fm->{'total'} --threads $options->{proc} --add-chrname > $fast_tmp/$fm->{'filename_core'}_hisat.sam\n";
	return $s, "$fast_tmp/$fm->{'filename_core'}_hisat.sam",
	  "$p/$fm->{'filename_core'}_hisat.sorted.bam";
}

sub create_sra_call {
	my $file = shift(@files);
	my $fm   = root->filemap($file);
	my $p    = $outpath;
	$p ||= "$fm->{'path'}/HISAT2_OUT";
	$mainOutpath = $p;
	mkdir("$p/") unless ( -d "$p/" );
	$fm->{'path'} = "" if ( $fm->{'path'} eq "./" );
	unless ( $fm->{'path'} =~ m/^\// ) { $fm->{'path'} = $dir . $fm->{'path'}; }
	my $s =
"hisat2 $mapper_options -x $genome --sra-acc $fm->{'total'} --threads $options->{proc} --add-chrname > $fast_tmp/$fm->{'filename_core'}_hisat.sam\n";
	return $s, "$fast_tmp/$fm->{'filename_core'}_hisat.sam",
	  "$p/$fm->{'filename_core'}_hisat.sorted.bam";
}

sub create_paired_call {
	my $file = shift(@files);
	my $pair = shift(@files);
	my $fm   = root->filemap($file);
	my $fm2  = root->filemap($pair);
	my $p    = $outpath;
	$p ||= "$fm->{'path'}/HISAT2_OUT";
	$mainOutpath = $p;
	$fm->{'path'}  = "" if ( $fm->{'path'} eq "./" );
	$fm2->{'path'} = "" if ( $fm2->{'path'} eq "./" );
	$fm->{'path'}  .= '/' unless ( $fm->{'path'} =~ m/\/$/ );
	$fm2->{'path'} .= '/' unless ( $fm2->{'path'} =~ m/\/$/ );
	mkdir("$p/") unless ( -d "$p/" );
	unless ( $fm->{'path'} =~ m/^\// ) { $fm->{'path'} = $dir . $fm->{'path'}; }
	my $s =
"hisat2 $mapper_options -x $genome -1 $fm->{'total'} -2 $fm2->{'total'} --threads $options->{proc} --add-chrname > $fast_tmp/$fm->{'filename_core'}_hisat.sam\n";
	return $s, "$fast_tmp/$fm->{'filename_core'}_hisat.sam",
	  "$p/$fm->{'filename_core'}_hisat.sorted.bam";
}

sub chk_cmd {
	my ( $cmd, @ofiles ) = @_;
	$this_outfile = $ofiles[0];
	return join( "\n",
		map { $SLURM->check_4_outfile( $_, @ofiles ) } split( /\n/, $cmd ) )
	  . "\n";
}
