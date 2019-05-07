#! /usr/bin/perl -w

=head1 LICENCE

  Copyright (C) 2016-12-07 Stefan Lang

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

    quantify_single_cell_NGS.pl
       -bams     :a list of bam files
       -gtf      :the gtf file to use
       -outfile  :the final StefansExpressionSet file
       -tmp_path :the temp path (default tmp in the outfiles's path)
       -amount   :the max amount of bam files quantified with one process (default =500)
       -IDXstats :the outfile from IDXstats_4_bams.pl script
       
       -help     :print this help
       -debug    :verbose output
       
       ## options for the Rsubread::featureCounts call
       
       -gtf_feature_type 
                 :the feature to quantify (default = exon)
       -gtf_attr_type    
                 :the attributeType (default = gene_id)
       -paired   :is the data paried
       -stranded :is the RNA seq been strand specific 
                  default 0; 1 stranded or 2 reversely stranded
                  
       -n        :how many processes to read the data (for the R process)
       
       ## options for the slurm batch process
       
       -slurm A lsens2017-3-2 p dell N 1 t 02:00:00
   
=head1 DESCRIPTION

  Having more than 600 samples in a NGS data set broke the Rsubread and therefore this tool got produced. 
  It will read any list of bam files in chunks of 500 samples per set and use a separate R instance each time.

  To get further help use 'quantify_single_cell_NGS.pl -help' at the comman line.

=cut

use Getopt::Long;
use Pod::Usage;

use stefans_libs::root;
use strict;
use warnings;

use FindBin;
my $plugin_path = "$FindBin::Bin";

my $VERSION = 'v1.0';

my (
	$help,     $debug,         $database, @bams,
	$gtf,      $tmp_path,      $amount,   $gtf_feature_type,
	$stranded, $gtf_attr_type, $paired,   $n,
	@slurm,    $IDXstats,      $outfile,  $ls_input
);

Getopt::Long::GetOptions(
	"-bams=s{,}"          => \@bams,
	"-gtf=s"              => \$gtf,
	"-gtf_feature_type=s" => \$gtf_feature_type,
	"-gtf_attr_type=s"    => \$gtf_attr_type,
	"-paired"             => \$paired,
	"-stranded=s"         => \$stranded,
	"-amount=s"           => \$amount,
	"-n=s"                => \$n,
	"-outfile=s"          => \$outfile,
	"-tmp_path=s"         => \$tmp_path,
	"-slurm=s{,}"         => \@slurm,
	"-IDXstats=s"         => \$IDXstats,

	"-help"  => \$help,
	"-debug" => \$debug
);

my $warn  = '';
my $error = '';

unless ( -f $bams[0] ) {
	eval {
		open( IN, "ls $bams[0] |" );
		$ls_input = $bams[0];
		@bams = map { chomp; $_ } <IN>;
		close(IN);
	};
}
unless ( -f $bams[0] ) {
	$error .= "the cmd line switch -bams is undefined!\n";
}
unless ( -f $gtf ) {
	$error .= "the cmd line switch -gtf is undefined!\n";
}
unless ( defined $gtf_feature_type ) {
	$gtf_feature_type = 'exon';
}
unless ( defined $amount ) {
	$amount = 500;
}
unless ( defined $gtf_attr_type ) {
	$gtf_attr_type = 'gene_id';
}

# paired - no checks necessary
unless ( defined $n ) {
	$n = 2;
}
unless ( defined $tmp_path ) {
	$tmp_path = 'tmp';
}
unless ( defined $outfile ) {
	$error .= "the cmd line switch -outfile is undefined!\n";
}
unless ($stranded) {
	$stranded = 0;
}
elsif ( !( $stranded == 1 or $stranded == 2 ) ) {
	$error .=
	  "-stranded may only have the values 0, 1 or 2 - not '$stranded'\n";
}

map {
	$warn .= "file $_ must not be given as relative path!\n"
	  if ( $_ =~ m/^\./ )
} @bams;

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

my ($task_description);

$task_description .= 'perl ' . $plugin_path . '/quantify_single_cell_NGS.pl';
$task_description .= ' -bams "' . join( '" "', @bams ) . '"'
  if ( defined $bams[0] );
$task_description .= " -gtf '$gtf'" if ( defined $gtf );
$task_description .= " -gtf_feature_type '$gtf_feature_type'"
  if ( defined $gtf_feature_type );
$task_description .= " -gtf_attr_type '$gtf_attr_type'"
  if ( defined $gtf_attr_type );
$task_description .= " -paired "   if ($paired);
$task_description .= " -stranded " if ($stranded);
$task_description .= " -n '$n'"    if ( defined $n );
$task_description .= " -tmp_path '$tmp_path'" unless ( $tmp_path eq "tmp" );
$task_description .= " -outfile '$outfile'" if ( defined $outfile );
$task_description .= " -amount $amount";

my $fm = root->filemap($outfile);

unless ( -d $fm->{'path'} ) {
	mkdir( $fm->{'path'} );
}

if ( $tmp_path eq "tmp" ) {
	$tmp_path = "$fm->{'path'}/tmp/";
}
unless ( -d $tmp_path ) {
	mkdir($tmp_path);
}

open( LOG, ">$outfile.log" ) or die $!;
print LOG $task_description . "\n";
close(LOG);

## Do whatever you want!
my ($max);
my $a = 1;
if ($paired) {
	$paired = "T";
}
else {
	$paired = "F";
}
use stefans_libs::SLURM;
my $slurm;
if (@slurm) {
	for ( my $i = 0 ; $i < @slurm ; $i += 2 ) {
		$slurm->{ $slurm[$i] } = $slurm[ $i + 1 ];
	}
}
$slurm->{'n'}     = $n;
$slurm->{'N'}     = 1 unless ( $slurm->{'N'} );
$slurm->{'debug'} = $debug;
$slurm->{'A'} ||= 'lsens2017-3-2';
$slurm->{'p'}             = 'dell';

$slurm = stefans_libs::SLURM->new($slurm);
$slurm->{'debug'} = $debug;

$slurm->{'SLURM_modules'} = [
	'ifort/2017.4.196-GCC-6.4.0-2.28', 'impi/2017.3.196',
	'R/3.4.3-X11-20171023'
];
$slurm->{'purge'} = 1;
unless( @slurm ) {
	$slurm->{'local'} = 1;
}

for ( my $i = 0 ; $i < @bams ; $i += $amount ) {
	open( OUT, ">$tmp_path/bams.$a.txt" )
	  or die "I could not create the bam file '$tmp_path/bams.$a.txt'\n$!\n";
	$max = $i + $amount - 1;
	$max = $#bams if ( $max > $#bams );
	print OUT
	  join( "\n", map { root->filemap($_)->{'total'}; } @bams[ $i .. $max ] )
	  . "\n";
	close(OUT);
	open( SCR, ">$tmp_path/quantify.$a.R" )
	  or die "I could not create the script file '$tmp_path/script.$a.R'\n$!\n";
	print SCR join(
		"\n",
		"#library( StefansExpressionSet)", &read_bams_R(),
		"library(Rsubread)",
"dat$a <- read.bams( '$tmp_path/bams.$a.txt', '$gtf', nthreads =  $n, GTF.featureType = '$gtf_feature_type',
					GTF.attrType = '$gtf_attr_type',isPairedEnd = $paired , strandSpecific = $stranded )",
		"save(dat$a, file='$tmp_path/Robject$a.RData')",
		"system('touch $tmp_path/Robject$a.RData.finished')"
	);
	close(SCR);
	if ( !-f "$tmp_path/Robject$a.RData" ) {
		
		$slurm->run(
			"R CMD BATCH $tmp_path/quantify.$a.R",
			root->filemap("$tmp_path/Robject$a.RData")
		);
		
	}
	else {
		print "outfile $tmp_path/Robject$a.RData does exist\n";
	}
	$a++;
}
$a--;

open( SCR, ">$tmp_path/sumup.R" )
  or die "I could not create the script file '$tmp_path/sumup.R'\n$!\n";
my @files = map { "$tmp_path/Robject$_.RData" } 1 .. $a;
print SCR join(
	"\n",
	(
		map {
"while ( ! file.exists( '$tmp_path/Robject$_.RData.finished' ) ) { Sys.sleep(10) }"
		} 1 .. $a
	),
	"load('$tmp_path/Robject1.RData')\n",
	"dat <- dat1",
	(
		map {
			my $id = $_;
			join( "\n",
				"load('$tmp_path/Robject$id.RData')",
				"dat\$counts <- cbind(dat\$counts, dat$id\$counts )",
				"dat\$stat <- cbind(dat\$stat, dat$id\$stat[,-1] )",
				"# the stat colnames get lost if dat$id\$stat[,-1] is a vector",
				"if ( class(dat$id\$stat[,-1]) == 'integer' ){",
"  colnames(dat\$stat)[ncol(dat\$stat)] = colnames(dat$id\$stat)[2]",
				"}" )
		} 2 .. $a
	),
	"save(dat, file='$outfile')",
);
close(SCR);
my $tmp = $slurm->{'local'};
$slurm->{'local'} = 1;
$slurm->run(
	"R CMD BATCH $tmp_path/sumum.R",
	root->filemap("$tmp_path/sumup.RData")
);
$slurm->{'local'} = $tmp;

## while all other things are processed lets work on the IDXstats....

if ( $bams[0] =~ m/sorted/ ) {

unless ( -f $IDXstats ) {
	$IDXstats = "$fm->{'path'}/$fm->{'filename_base'}_IDXstats.xls";
	$IDXstats =~ s!//!/!g;
	if ( defined $ls_input ) {
		@bams = ($ls_input);
	}
	my $cmd =
	    "IDXstats_4_bams.pl -bams '"
	  . join( "' '", @bams )
	  . "' -outfile $IDXstats -n $n";
	print
"As I did not get the IDXstats summary file I will create one($IDXstats)\n$cmd\n";
	system($cmd );
}

## convert the IDX stats into a usable R object

my $failedCellsRobj = "$fm->{'path'}/FailedSamples.RData";

unless ( -f $failedCellsRobj ) {
	open( RSCRIPT, ">$fm->{'path'}/idxstats2failedSamples.R" )
	  or die
"I could not create the checkup R script '$fm->{'path'}/idxstats2failedSamples.R'\n$!\n";

	print RSCRIPT "library(stringr)
t <- read.delim('$IDXstats')
FailedSamples <- as.character(
   t[
        which(is.na(apply( t, 1, function(x) { sum(as.numeric(x[grep ('chr[XY1-9][1-9]*$',colnames(t))]))}))==T)
        ,1
   ]
)
names(FailedSamples) <- str_extract( FailedSamples, 'SRR\\d*')
save(FailedSamples,file='$failedCellsRobj' )
";
	close(RSCRIPT);

	system("R CMD BATCH $fm->{'path'}/idxstats2failedSamples.R");
}

## now we can think about summing up the results / hopefully most quantifications have been finished up to now ;-)

}else {
	warn "The idxstats could not be created as the bam files were not sorted.\n";
}
		
print
"Please wait for the sumuo.R process to finish\nThe file '$outfile' will then contain the R counts object that can be further processed in any R script\n";

sub read_bams_R {
	return
	  "read.bams <- function ( bamFiles, annotation, GTF.featureType='exon', 
		GTF.attrType = 'gene_id', isPairedEnd = FALSE, nthreads = 2, as.obj=F, strandSpecific = 0 ) {
	if (file.exists(bamFiles)){
		bamFiles <- readLines(bamFiles)
	}

	counts <- featureCounts(files =bamFiles,annot.ext = annotation ,isGTFAnnotationFile = TRUE,GTF.featureType = GTF.featureType,
		GTF.attrType = GTF.attrType,allowMultiOverlap=T, isPairedEnd =isPairedEnd , nthreads = nthreads, strandSpecific = strandSpecific )
	save(counts, file='count_object.RData' )
	if ( ! as.obj ) {
		counts
	}else {
		counts.tab <- cbind(counts\$annotation,counts\$counts)  # combine the annotation and counts
		counts.tab
	}
}
"
}
