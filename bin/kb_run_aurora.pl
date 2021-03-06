#!/sw/easybuild/software/Perl/5.26.1-GCCcore-6.4.0/bin/perl -w

=head1 LICENCE

  Copyright (C) 2019-10-24 Stefan Lang

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
   
   binCreate.pl from git@github.com:stela2502/Stefans_Lib_Esentials.git commit c35cfea822cac3435c5821897ec3976372a89673
   

=head1  SYNOPSIS

    kb_run_aurora.pl
       -files          :a list of fastq files for one sample
       -sname          :the sample  name for this set of fastqs   
       -outpath        :the outpath for this data
       -x              :the 10x library version like 10xv1, 10xv2 or 10xv3
       -options        :SLURM options you can specify more entries to thatn    
       		
       		n    :number of cores per node (default = 10 )
            N    :number of nodes (default =1)
            t    :time untill the script is stopped (default =02:00:0 (2h))
            proc :hisat 2 number of threads (default = n*N )
     partitition :explicitely call for a partitition in the SLURM environment (slurm p option)
                  should not be necessary!
                         
       -fast_tmp       :the nodes should have a fast local tmp folder that could be used for 
                         the intermediate files (default '$SNIC_TMP')
                         currently unused
       
       -mapper_options :currently unused
       -index          :the kallisto index files
       -local          :run the analysis mouse/Genecode_M25/on the frontend (DO NOT USE!)

       -help           :print this help
       -debug          :verbose output
   
=head1 DESCRIPTION

  Start a kb wrapper call using kallisto to create a h5ad file.

  To get further help use 'kallisto_run_aurora.pl -help' at the comman line.

=cut

use Getopt::Long;
use Pod::Usage;

use strict;
use warnings;

use stefans_libs::root;
use stefans_libs::SLURM;

use Cwd 'abs_path';

use FindBin;
my $plugin_path = "$FindBin::Bin";

my $VERSION = 'v1.0';


my ( $help, $debug, $database, @files, $x, $outpath, $options, @options, $fast_tmp, $mapper_options, $local, $organism, $sname, $tmp);

Getopt::Long::GetOptions(
     "-files=s{,}"    => \@files,
	 "-outpath=s"    => \$outpath,
	 "-x=s" => \$x,
     "-options=s{,}"    => \@options,
	 "-fast_tmp=s"    => \$fast_tmp,
	 "-mapper_options=s"    => \$mapper_options,
	 "-organism=s" => \$organism,
	 "-sname=s"     => \$sname,
   	 "-local"            => \$local,
	 "-help"             => \$help,
	 "-debug"            => \$debug
);

my $indexPath = "/projects/fs1/common/genome/lunarc/10Xindexes/kallisto";

my $warn = '';
my $error = '';

unless ( defined $files[0]) {
	$error .= "the cmd line switch -files is undefined!\n";
}
else {
	foreach my $f ( @files ) {
		if ( ! -f $f ) {
			$error .="file not found '$f'\n";
		}
	}
}

unless ( defined $outpath) {
	$warn .="The outpath has been set to the sname\n";
	$outpath = $sname;
#	$error .= "the cmd line switch -outpath is undefined!\n";
}

unless ( defined $options[0]) {
	$error .= "the cmd line switch -options is undefined!\n";
}
unless ( defined $fast_tmp) {
	$fast_tmp = '$SNIC_TMP';
	#$error .= "the cmd line switch -fast_tmp is undefined!\n";
}
unless ( defined $x ) {
	$error .= "I need the 10x lib version '-x'!\n";
}
unless ( defined $mapper_options) {
	$mapper_options = '';
	#$error .= "the cmd line switch -mapper_options is undefined!\n";
}
unless ( defined $organism) {
	$error .= "the cmd line switch -organism is undefined (either 'human' or 'mouse')!\n";
}
unless ( defined $sname) {
	$error .= "the cmd line switch -sname is undefined!\n";
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
$options->{'n'}    ||= 20;
$options->{'N'}    ||= 1;
$options->{'t'}    ||= '20:00:00';
$options->{'partitition'} ||= 'dell';
$options->{'A'} ||= 'lsens2018-3-3';
###


my ( $task_description);


$task_description .= 'perl '.$plugin_path .'/kb_run_aurora.pl';
$task_description .= ' -files "'.join( '" "', @files ).'"' if ( defined $files[0]);
$task_description .= " -outpath '$outpath'" if (defined $outpath);
$task_description .= ' -options "'.join( '" "', @options ).'"' if ( defined $options[0]);
$task_description .= " -fast_tmp '$fast_tmp'" if (defined $fast_tmp);
$task_description .= " -mapper_options '$mapper_options'" if (defined $mapper_options);
$task_description .= " -organism '$organism'" if (defined $organism);
$task_description .= ' -local ' if ($local);
$task_description .= " -sname '$sname'" if (defined $sname);


for ( my $i = 0 ; $i < @options ; $i += 2 ) {
	$options[ $i + 1 ] =~ s/\n/ /g;
	$options->{ $options[$i] } = $options[ $i + 1 ];
}
###### default options ########
#$options->{'something'} ||= 'default value';
##############################
mkdir( $outpath ) unless ( -d $outpath );
open ( LOG , ">$outpath/".$$."_kb_run_aurora.pl.log") or die $!;
print LOG $task_description."\n";
close ( LOG );

my ( $i, $g, $c1, $c2 );
if ( $organism eq 'human') {
	die ( "please add the human option into the script!");
}
elsif ( $organism eq 'mouse') {
	## set the values
	#  -i INDEX              Path to kallisto index
	#  -g T2G                Path to transcript-to-gene mapping
	#  -x TECHNOLOGY         Single-cell technology used (`kb --list` to view)
	#  -c1 T2C               Path to cDNA transcripts-to-capture
	#  -c2 T2C               Path to intron transcripts-to-captured
	#/projects/fs1/common/genome/lunarc/10Xindexes/kallisto/mouse/Genecode_M25/
	$indexPath = "$indexPath/mouse/GRCm38.p6";
	$i  = "$indexPath/ensembl.GRCm38.p6.100.idx";
	$g  = "$indexPath/transcripts_to_genes.txt";
	$c1 = "$indexPath/cDNA_transcripts_2_capture.txt";
	$c2 = "$indexPath/introns_to_capture.txt";
	
}
else {
	die ("Sorry, but the organism '$organism' is not supported at the moment - fix the script?")
}

## Do whatever you want!

#!/bin/sh
#SBATCH -n 20
#SBATCH -N 1
#SBATCH -t 12:00:00
#SBATCH -A lu2018-2-44
#SBATCH -p lu
#SBATCH  --reservation=lu2018-2-44
#SBATCH -J Org32d_H9_409B2
#SBATCH -o Org32d_H9_409B2.%j.out
#SBATCH -e Org32d_H9_409B2.%j.err
#module purge
#module load GCC/8.2.0-2.31.1 OpenMPI/3.1.3 kallisto/0.46.0 bustools/0.39.3 Anaconda3/2019.07
#
#loompy fromfq Org32d_H9_409B2.loom Org32d_H9_409B2 /lunarc/nobackup/users/ssoneji/LCA/human_GRCh38_gencode.v31.600 Org32d_H9_409B2.meta.tab Org32d_H9_409B2_R1.fq.gz Org32d_H9_409B2_R2.fq.gz 

@files = map { File::Spec->rel2abs($_) } @files;

$outpath = File::Spec->rel2abs($outpath);

my $SLURM = stefans_libs::SLURM->new( $options, 0 );
$SLURM->{'debug'} = 1 if ($debug);
$SLURM->{'local'} = 1 if ($local);
$SLURM->{'max_jobs'} = 10;

$SLURM->{'SLURM_modules'} = [
	'Singularity/default'	
	];
$SLURM->{'purge'} = 1;

## create the metadata.tab file
## example on their web page:
#name    technology  targetnumcells
#1kPBMC  10xv3       1000
open( META, ">$outpath/$sname.metadata.tab") or die $!;
print META join("\t", qw(name    technology  targetnumcells) )."\n";
print META join("\t", $sname, $x, 10000 )."\n";
close(META);

## enable to use the latest singularity image for the job:
my $singularity = '/projects/fs1/common/software/SingSingCell/1.0/SingleCells_v1.0.sif';
if ( ! -f $singularity ){
	die "the expected sigularity image '$singularity' des not exist\n";
}
my $job = { 
	'fm' => "$outpath/$sname.sh", 
	'cmd' => "singularity exec -B/projects:/projects $singularity ". #"cd $outpath &&".
		"kb count -o $outpath/$sname -t $options->{'n'} --workflow lamanno --h5ad -i $i -g $g -x $x -c1 $c1 -c2 $c2 ".join(" ", sort @files ),
	'outfile' => "$outpath/$sname.h5ad"
};

$tmp = $SLURM->run( $job->{'cmd'}, $job->{'fm'}, $job->{'outfile'} );

if ( $tmp ) {
	print "Sucessfully submitted the job\n";
}


