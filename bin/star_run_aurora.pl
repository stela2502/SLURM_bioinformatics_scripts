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

    star_run_aurora.pl
       -files     :the input fastq or sra files
       -outpath   :the outpath for the mapped files
                   default to subpath of infile named HISAT2_OUT
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

       -mapper_options :specific options for the star process as one string
                        e.g. -mapper_options '--max-seeds 500'

       -genome       :the star genome information path
       -coverage     :the chromome length file
       -paired       :analyse paried fasta files
                      every first file == read 1 every second == read2

       -dropDuplicates :use picard to drop duplicates
       
       -justMapping  :Just map and sort the bam files; not create bigwig files
       -noSort       :more restrictive than justMapping as the bams are not even sorted

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
	$help,         $debug,    $database, @files,
	$local,        $options,  @options,  $dropDuplicates,
	$genome,       $coverage, $paired,   $fast_tmp,
	$bigwigTracks, $sra,      $outpath,  $max_jobs,
	$justMapping,  $mapper_options,      $noSort,
);

Getopt::Long::GetOptions(
	"-files=s{,}"       => \@files,
	"-options=s{,}"     => \@options,
	"-genome=s"         => \$genome,
	"-coverage=s"       => \$coverage,
	"-outpath=s"        => \$outpath,
	"-paired"           => \$paired,
	"-bigwigTracks=s"   => \$bigwigTracks,
#	"-sra"              => \$sra,
	"-max_jobs=s"       => \$max_jobs,
	"-dropDuplicates"   => \$dropDuplicates,
	"-mapper_options=s" => \$mapper_options,
	"-fast_tmp=s"       => \$fast_tmp,
	"-justMapping"      => \$justMapping,
	"-local"            => \$local,
	"-noSort" => \$noSort,

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
if (  $sra ) {
	Carp::confes( "sra option does not exists for STAR!");
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
#$options->{'p'}    ||= $options->{'proc'};

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
$task_description .= " -justMapping"          if ($justMapping);
$task_description .= " -noSort"               if ($noSort);

$task_description .= " -debug"                if ($debug);

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

my $SLURM = stefans_libs::SLURM->new( $options, 0 );
$SLURM->{'debug'} = 1 if ($debug);
$SLURM->{'local'} = 1 if ($local);
$SLURM->{'max_jobs'} = $max_jobs;
my $BAM = stefans_libs::scripts::BAM->new($options);
$BAM->{'p'} = $options->{'n'};

## kick all SLURM options that should not be used for the hisat run
foreach (qw(n N t mem)) {
	delete( $options->{$_} );
}

$fm = root->filemap( $files[0] );

$SLURM->{'purge'} = 1;

unless ($local) {
	$SLURM->{'SLURM_modules'} = [
	'GCC/5.4.0-2.26',  'OpenMPI/1.10.3', 'STAR/2.6.0c',
	'SAMtools/1.4', 'BEDTools/2.26.0'
		#'GCC/4.9.3-2.25', 'OpenMPI/1.10.2',
		#'icc/2016.1.150-GCC-4.9.3-2.25', 'impi/5.1.2.150', 'SAMtools/1.3.1',
		#'HISAT2/2.0.4',
		#'BEDTools/2.25.0', 'Java/1.8.0_72', 'ucsc-tools/R2016a',

		#	stefans_libs::scripts::BAM->SLURUM_load(),
	];
}

$tmp = $SLURM->define_Subscript();
$tmp->{'purge'} = 1;
unless ($local) {
	$tmp->{'SLURM_modules'} = [
		'GCC/4.9.3-2.25',                'OpenMPI/1.10.2',
		'icc/2016.1.150-GCC-4.9.3-2.25', 'impi/5.1.2.150',
		'BEDTools/2.25.0',               'Java/1.8.0_72',
		'picard/2.8.2',                  'ucsc-tools/R2016a',
	];
}

@files = map { $_ =~ s/^\.\///; $_ } @files;
my $submitted = 0;
my @jobs;
my @tmp;
my $preprocessed_outfile_to_infile = "";
## create the jobs to sort them by file size later on ;-)

while ( scalar(@files) ) {
	my $fm = root->filemap( $files[0] )
	  ;    ## the files will be depleted by the create_call function!
	my (@cmd, @ofiles);
	
	($cmd[0], $this_outfile ) = &create_call(); ## this is a temporary step
	$cmd[0] = "date\n$cmd[0]date\n"; ## Would be nice to see the execution time in the log...
#	if ( $noSort ){
#		@tmp = $BAM->convert_sam_2_bam( $this_outfile );
#		$tmp[0] = $cmd[0].$tmp[0];
#		$cmd[0] = &chk_cmd( @tmp );
#	}else {
#		@tmp = $BAM->convert_sam_2_sorted_bam($this_outfile );
#		$tmp[0] = $cmd[0].$tmp[0];
#		$cmd[0] = &chk_cmd( @tmp );
#	}
	push ( @ofiles, "$this_outfile" );
		
	if ($dropDuplicates){
		$cmd[1] .= &chk_cmd( create_picard_call($this_outfile) );
		push ( @ofiles, "$this_outfile" );
	}
	
	if (!( $justMapping or $noSort ) ){ ## run this
	$cmd[1] .=
	  &chk_cmd(
		$BAM->convert_sorted_bam_2_bedGraph( $this_outfile, $coverage ) );
	$cmd[1] .=
	  &chk_cmd( $BAM->convert_bedGraph_2_bigwig( $this_outfile, $coverage ) );
	  push ( @ofiles, "$this_outfile" );
	}
	$BAM->{'Outfile_exists'} = 0;
	my @tmp = map {  &chk_cmd_no_cp( $BAM->recover_files( $mainOutpath, $_ ), $_ ) } @ofiles;
	my $null_all = 0;
	if ( $BAM->{'Outfile_exists'} == scalar(@ofiles) ) {
		## OK we can stop here...
		print ( "outfile(s) ". join(", ", @ofiles). " do all exists - no need to run the script\n" );
		$null_all = 1;
	}
	map { $cmd[1] .= $_ } @tmp;
	if ( $null_all ) {
		## nothing has to be run - all outfiles are present!
		$cmd[0] = join( "\n", map { if ( $_ =~m/^#/) { $_} else { "#$_"} } split("\n", $cmd[0] ) );
		$cmd[1] = join( "\n", map { if ( $_ =~m/^#/) { $_} else { "#$_"} } split("\n", $cmd[1] ) );
	}

	push( @jobs, { 'outfile' => "$this_outfile", 'cmd' => [@cmd], 'fm' => $fm, 'run' => !$null_all  } );
}

foreach my $job ( @jobs[ &FileSizeOrder(@jobs) ] ) {
	#print "\$exp = " . root->print_perl_var_def($job) . ";\n";
	my $tmp = $SLURM->run( $job->{'cmd'}, $job->{'fm'}, $job->{'outfile'} );
	$submitted++ if ( $tmp == 1 );
	if ( $submitted >= $max_jobs ) {
		$submitted -= 50;
		$SLURM->wait_for_last_finished( $job->{'outfile'} );
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
exit 0;

sub byFileSize {
	-s $b->{'fm'}->{'total'} <=> -s $a->{'fm'}->{'total'};
}

sub FileSizeOrder {
	my @files = @_;
	my $i     = 0;
	my $order = { map { $_->{'fm'}->{'total'} => $i++ } @files };
	my @ret;
	foreach ( sort byFileSize @files ) {
		push( @ret, $order->{ $_->{'fm'}->{'total'} } );
	}
	return @ret;
}

sub create_picard_call {
	my ($file) = @_;
	my $fm     = root->filemap($file);
	
	#warn "the \$preprocessed_outfile_to_infile would be this: '$preprocessed_outfile_to_infile'\n";
	my $cmd = "picard";
	if ($local) {
		$cmd = "picard-tools";
	}	
	my $s =
	    "$cmd MarkDuplicates INPUT="
	  . $fm->{'total'}
	  . " OUTPUT="
	  . "$fm->{'path'}/$fm->{'filename_core'}_picard_deduplicated.bam"
	  . " REMOVE_DUPLICATES=TRUE M="
	  . "$fm->{'path'}/$fm->{'filename_core'}_picard_metrix.txt";
	return $s, "$fm->{'path'}/$fm->{'filename_core'}_picard_deduplicated.bam";
}

sub create_call {
	return &create_paired_call() if ($paired);
	return &create_sra_call()    if ($sra);
	my $file = shift(@files);
	my $fm   = root->filemap($file);
	my $p    = $outpath;
	$p ||= "$fm->{'path'}/STAR_OUT";
	$mainOutpath = $p;
	mkdir("$p/") unless ( -d "$p/" );
	$fm->{'path'} = "" if ( $fm->{'path'} eq "./" );
	unless ( $fm->{'path'} =~ m/^\// ) { $fm->{'path'} = $dir . $fm->{'path'}; }
	my $unzip = "";
	if ( $fm->{'filename_ext'} eq "gz" ) {
		$unzip = '--readFilesCommand zcat'
	}
	my $encode_opts = "--outFilterType BySJout --outFilterMultimapNmax 20 --alignSJoverhangMin 8"
		.	" --alignSJDBoverhangMin 1 --outFilterMismatchNmax 999 --outFilterMismatchNoverReadLmax 0.04"
		. 	" --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000 --outSAMmultNmax 1";
	my $s =
"STAR $mapper_options --genomeDir $genome --readFilesIn $fm->{'total'} --runThreadN $options->{proc} "
. "$unzip $encode_opts --outFileNamePrefix $fast_tmp/$fm->{'filename_core'} --outSAMtype BAM SortedByCoordinate --twopassMode Basic \n";
	## here I need to know what we finally what to get
	
	return $s, "$fast_tmp/$fm->{'filename_core'}Aligned.sortedByCoord.out.bam", "$fast_tmp/$fm->{'filename_core'}Aligned.sortedByCoord.out.bam"  ;
}

sub create_sra_call {
	my $file = shift(@files);
	my $fm   = root->filemap($file);
	my $p    = $outpath;
	$p ||= "$fm->{'path'}/STAR_OUT";
	$mainOutpath = $p;
	mkdir("$p/") unless ( -d "$p/" );
	$fm->{'path'} = "" if ( $fm->{'path'} eq "./" );
	unless ( $fm->{'path'} =~ m/^\// ) { $fm->{'path'} = $dir . $fm->{'path'}; }
	my $s =
"hisat2 $mapper_options -x $genome --sra-acc $fm->{'total'} --threads $options->{proc} --add-chrname > $fast_tmp/$fm->{'filename_core'}_hisat.sam\n";
	return $s, "$fast_tmp/$fm->{'filename_core'}_hisat.sam" ,"$p/$fm->{'filename_core'}_hisat.sorted.bam", "$p/$fm->{'filename_core'}_hisat.bam" ;
}

sub create_paired_call {
	my $file = shift(@files);
	my $pair = shift(@files);
	my $fm   = root->filemap($file);
	my $fm2  = root->filemap($pair);
	my $p    = $outpath;
	$p ||= "$fm->{'path'}/STAR_OUT";
	$mainOutpath = $p;
	$fm->{'path'}  = "" if ( $fm->{'path'} eq "./" );
	$fm2->{'path'} = "" if ( $fm2->{'path'} eq "./" );
	$fm->{'path'}  .= '/' unless ( $fm->{'path'} =~ m/\/$/ );
	$fm2->{'path'} .= '/' unless ( $fm2->{'path'} =~ m/\/$/ );
	mkdir("$p/") unless ( -d "$p/" );
	unless ( $fm->{'path'} =~ m/^\// ) { $fm->{'path'} = $dir . $fm->{'path'}; }
	my $unzip = "";
	if ( $fm->{'filename_ext'} eq "gz" ) {
		$unzip = '--readFilesCommand zcat'
	}
	my $encode_opts = "--outFilterType BySJout --outFilterMultimapNmax 20 --alignSJoverhangMin 8"
		.	" --alignSJDBoverhangMin 1 --outFilterMismatchNmax 999 --outFilterMismatchNoverReadLmax 0.04"
		. 	" --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000 --outSAMmultNmax 1";
	my $s =
"STAR $mapper_options --genomeDir $genome --readFilesIn $fm->{'total'} $fm2->{'total'} --runThreadN $options->{proc} "
. "$unzip $encode_opts --outFileNamePrefix $fast_tmp/$fm->{'filename_core'} --outSAMtype BAM SortedByCoordinate --twopassMode Basic \n";
	
	return $s, "$fast_tmp/$fm->{'filename_core'}Aligned.sortedByCoord.out.bam", "$fast_tmp/$fm->{'filename_core'}Aligned.sortedByCoord.out.bam"  ;

}

sub chk_cmd {
	my ( $cmd, @ofiles ) = @_;
	$this_outfile = $ofiles[0];
	$preprocessed_outfile_to_infile = undef;
	my @return = map { $SLURM->check_4_outfile( $_,@ofiles, &final_outfile($mainOutpath, @ofiles) ) } split( /\n/, $cmd ) ;
	my $regen = 0;
	map { if ( $_ =~ m/^#/ ){ $regen = 1 } } @return;
	if ( $regen and ref($preprocessed_outfile_to_infile) eq "HASH" ) {
		push ( @return, keys %$preprocessed_outfile_to_infile );
	}
	return join( "\n", @return )
	  . "\n";
}

sub chk_cmd_no_cp {
	my ( $cmd, @ofiles ) = @_;
	$this_outfile = $ofiles[0];
	$preprocessed_outfile_to_infile = undef;
	my @return = map { $SLURM->check_4_outfile( $_,@ofiles, &final_outfile($mainOutpath, @ofiles) ) } split( /\n/, $cmd ) ;
	return join( "\n", @return )
	  . "\n";
}


=head final_outfile( $final_path, @files )

checks if the outfiles are present at the moment and therefore do not need to get created.

=cut

sub final_outfile {
	my ( $final_path, @files) = @_;

	my @ret;
	foreach my $fm (  @files ) {
		$fm = &fm_check($fm);
		my $f  = $fm->{'path'} . "/" . $fm->{'filename'};
		my $final = $final_path."/".$fm->{'filename'};
		unless ( -f $final ) {
			## this needs to be produced
			push( @ret, 1)
		}else { push ( @ret , 0 )}
#		unless ( $final eq $f ){
#			# create a little bash script that does some error checking..
#			my $s = "if [ -f '$final' ];then \n  cp $final $f\n"
#			. "else\n  ( >&2 echo 'file $final was not produced - files in the same path:')\n  for filename in /home/user/*\n  do\n    (>&2 echo \$filename) \n  done;\nfi\n";
#			push ( @ret, $s);
#		}else {
#			push ( @ret, $f );
#		}
	}
	return @ret;
}

sub fm_check {
	my ( $fm ) = @_;
	unless ( ref($fm) eq "HASH") {
 		if ( $fm =~m/^(\$\w*)\/(.*)$/ ) {
 			## oops the tmp file variable!
 			my $tmp = $1;
 			$fm = root->filemap( $2 );
 			$fm->{'path'} = $tmp;
 		}else {
 			$fm = root->filemap( $fm );
 		}
 	}
 	return $fm;
}
