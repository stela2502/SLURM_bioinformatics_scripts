#! /usr/bin/perl -w

use Getopt::Long;
use Pod::Usage;
use stefans_libs::root;
use stefans_libs::SLURM;
use stefans_libs::scripts::BAM;

=head1 LICENCE

  Copyright (C) 2016-05-31 Stefan Lang

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

    bowtie_run_aurora.pl
       -files     :a list of fastq or fastq.gz files   
       -options   :additional options in the format
                   key_1 value_1 key_2 value_2 ... key_n value_n
                   recommended: n 10 N 1 mem 4000M t '02:02:00' p 10
       -genome    :the genome definition as used by bowtie -X
       -coverage  :the genome coverage file to create bigwig tracks
       -paired    :if you have paired data to map use this option and
                   give me the pired files one after the other
                   
       -bowtie_options:
                   All additional bowtie options as one string e.g.
                   '-m 500 -v2 --best --strata'
                   
       -bigwigTracks :If I have a coverage file I can create the 
                      bigwig tracks information and stoire it there

       -help      :print this help
       -debug     :verbose output
    
    required options:
    
    n   :number of cores to request from SLURM (rec. 10)
    N   :number of nodes to request from SLURM (rec. 1)
    mem :minimum memory from SLURM (rec. 4000M for e.g. human genome)
    t   :the time the job should get to finish - do not set too short
         rec 02:00:00 == 2h
  
=head1 DESCRIPTION

  This script takes a  list of fastq or fastq.gz files, created a folder named bowtie below these files and outputs the processed data there. Use -debug to not start the SBIN script.

  To get further help use 'bowtie_run_aurora.pl -help' at the comman line.

=cut

use strict;
use warnings;

use FindBin;
my $plugin_path = "$FindBin::Bin";

my $VERSION = 'v1.0';

my (
	$help,    $debug,  $database, @files,    $options,$bowtie_options,
	@options, $genome, $paired,   $coverage, $bigwigTracks
);

Getopt::Long::GetOptions(
	"-files=s{,}"     => \@files,
	"-options=s{,}"   => \@options,
	"-genome=s"       => \$genome,
	"-coverage=s"     => \$coverage,
	"-paired"         => \$paired,
	"-bigwigTracks=s" => \$bigwigTracks,
	"-bowtie_options=s" => \$bowtie_options,

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
else {
	my $tmp = root->filemap($genome);
	unless ( -d $tmp->{'path'} ) {
		$error .=
		  "the cmd line switch -genome '$tmp->{'path'}'\npath not found!\n";
	}
}
$bowtie_options ||= '';
$coverage ||= '';
unless ( -f $coverage ) {
	$warn .=
"bigwig files can not be created unless you give me the right genome coverage file (-coverage)\n";
}
else {
	unless ( defined $bigwigTracks ) {
		$error .= "the cmd line switch -bigwigTracks  is undefined!\n";
	}
	else {
		my $tmp = root->filemap($bigwigTracks);
		unless ( -d $tmp->{'path'} ) {
			mkdir( $tmp->{'path'} );
		}
	}
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

my ($task_description);

$task_description .= 'perl ' . $plugin_path . '/bowtie_run_aurora.pl';
$task_description .= '       -files "' . join( '" "', @files ) . '"'
  if ( defined $files[0] );
$task_description .= '       -options "' . join( '" "', @options ) . '"'
  if ( defined $options[0] );
$task_description .= " -genome '$genome'"     if ( defined $genome );
$task_description .= " -paired"               if ($paired);
$task_description .= " -coverage '$coverage'" if ( -f $coverage );
$task_description .= " -bigwigTracks '$bigwigTracks'"
  if ( defined $bigwigTracks );
$task_description .= " -bowtie_options '$bowtie_options'" if ( $bowtie_options =~ m/\w/ );

for ( my $i = 0 ; $i < @options ; $i += 2 ) {
	$options[ $i + 1 ] =~ s/\n/ /g;
	$options->{ $options[$i] } = $options[ $i + 1 ];
}
#print "\$exp =  " . root->print_perl_var_def($options) . ";\n";
## Do whatever you want!
my ( $cmd, $fm, $this_cmd, @big_wig_urls, $tmp, $this_outfile );

$options->{'n'} ||= 10;
$options->{'N'} ||= 1;
$options->{'t'} ||= '02:00:00';

my $SLURM = stefans_libs::SLURM->new($options);
$SLURM->{'debug'} = 1 if ($debug);

## kick all SLURM options that should not be used for the bowtie
foreach (qw(n N t mem A)) {
	delete( $options->{$_} );
}

my $BAM = stefans_libs::scripts::BAM->new($options);

##### load the modules that are needed first!
## therefore I need to create a script file and run that using bash
$fm = root->filemap( $files[0] );
$SLURM->{'purge'} = 1;
$SLURM->{'SLURM_modules'} = [
'icc/2015.3.187-GNU-4.9.3-2.25',  'impi/5.0.3.048 Bowtie/1.1.2', 'SAMtools/0.1.20'
];



my $files_modified = 0;
$cmd = "bowtie $bowtie_options";
foreach my $key ( keys %$options ) {
	$cmd .= " -$key $options->{$key}";
}
$cmd .= " --sam $genome ";
for ( my $i = 0 ; $i < @files ; $i++ ) {
	if ( $files[$i] =~ m/\s/ ) {
		$tmp = $files[$i];
		$tmp =~ s/\s+/_/g;
		warn "I rename the file $files[$i] to $tmp\n";
		print "mv '$files[$i]' $tmp\n";
		system("mv '$files[$i]' $tmp");
		$files_modified = 1;
		$files[$i] = $tmp;
	}
}

Carp::confess("Files were moved/renamed - please re-run!\n")
  if ($files_modified);

if ($paired) {
	for ( my $i = 0 ; $i < @files ; $i += 2 ) {
		$this_cmd = $cmd;
		$fm       = root->filemap( $files[$i] );
		$this_outfile =
		  "$fm->{'path'}/bowtie/$fm->{'filename_core'}_bowtie.sam";
		unless ( -d $fm->{'path'} . "/bowtie" ) {
			mkdir( $fm->{'path'} . "/bowtie" );
		}
		$this_cmd .= "-1 '$files[$i]' -2 '" . $files[ $i + 1 ] . "'";
		$this_cmd .= " > $this_outfile";
		$this_cmd = $SLURM->check_4_outfile( $this_cmd,
			"$fm->{'path'}/bowtie/$fm->{'filename_core'}_bowtie.sorted.bam" );
		$this_cmd = $SLURM->check_4_outfile( $this_cmd,
			"$fm->{'path'}/bowtie/$fm->{'filename_core'}_bowtie.sam" );
		$this_cmd .= &chk_cmd($BAM->convert_sam_2_sorted_bam($this_outfile));
		$this_cmd .= "#". &chk_cmd($BAM->convert_sorted_bam_2_bedGraph($this_outfile,$coverage));
		$this_cmd .= "#". &chk_cmd($BAM->convert_bedGraph_2_bigwig($this_outfile,$coverage));
		$SLURM->run( $this_cmd, $fm, $this_outfile );
	}
}
else {
	for ( my $i = 0 ; $i < @files ; $i++ ) {
		$this_cmd = $cmd;
		$fm       = root->filemap( $files[$i] );
		$this_outfile =
		  "$fm->{'path'}/bowtie/$fm->{'filename_core'}_bowtie.sam";
		unless ( -d $fm->{'path'} . "/bowtie" ) {
			mkdir( $fm->{'path'} . "/bowtie" );
		}
		$this_cmd .= "'$files[$i]'";
		$this_cmd .=
		  " > $fm->{'path'}/bowtie/$fm->{'filename_core'}_bowtie.sam";
		$this_cmd = &chk_cmd( $this_cmd,
			"$fm->{'path'}/bowtie/$fm->{'filename_core'}_bowtie.sorted.bam" );
		$this_cmd = $SLURM->check_4_outfile( $this_cmd,
			"$fm->{'path'}/bowtie/$fm->{'filename_core'}_bowtie.sam" );
		$this_cmd .= &chk_cmd(
			$BAM->convert_sam_2_sorted_bam("$fm->{'path'}/bowtie/$fm->{'filename_core'}_bowtie.sam")
		)."\n";
		$this_cmd .= "#".join("\n#", split("\n",&chk_cmd($BAM->convert_sorted_bam_2_bedGraph($this_outfile,$coverage)) ) )."\n";
		$this_cmd .= "#".join("\n#", split("\n", &chk_cmd($BAM->convert_bedGraph_2_bigwig($this_outfile,$coverage))));
		$SLURM->run( $this_cmd, $fm );
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

if ($debug) {
	$fm = root->filemap( $files[0] );
	print
"before you run the sbatch calls make sure you have loaded all required modules.\n"
	  . "I recommend: 'bash $fm->{'path'}/InitializeSLURMenv.sh'\n";
}

sub chk_cmd {
	my ( $cmd, $ofile ) = @_;
	$this_outfile = $ofile;
	return join( "\n",
		map { $SLURM->check_4_outfile( $_, $ofile ) } split( /\n/, $cmd ) )
	  . "\n";
}

