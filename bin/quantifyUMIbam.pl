#! /usr/bin/perl -w

=head1 LICENCE

  Copyright (C) 2016-07-07 Stefan Lang

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

    quantifyUMIbam.pl
       -bam_files     :the aligned bam or sam files
       -organism      :the genome organism (to identify the genome in the database)
       -version       :the genome version (to identify the genome in the database)
       -genome_bed_file  :a previousely created genome information bed file
       
       -outpath       :the outpath for all result files

       -help           :print this help
       -debug          :verbose output
   
=head1 DESCRIPTION

  read a set of bam files and quantify the UMI baed bam files (produced from fastq files splitted with spliceUMI_fastq.pl)

  To get further help use 'quantifyUMIbam.pl -help' at the comman line.

=cut

use Getopt::Long;
use Pod::Usage;
use stefans_libs::database::genomeDB;
use Parallel::ForkManager;

use strict;
use warnings;

use FindBin;
my $plugin_path = "$FindBin::Bin";

my $VERSION = 'v1.0';

my ( $help, $debug, $database, @bam_files, $organism, $version,
	$genome_bed_file, $outpath, $proc );

Getopt::Long::GetOptions(
	"-bam_files=s{,}"    => \@bam_files,
	"-organism=s"        => \$organism,
	"-version=s"         => \$version,
	"-genome_bed_file=s" => \$genome_bed_file,
	"-outpath=s"         => \$outpath,
	"-proc=s" => \$proc,

	"-help"  => \$help,
	"-debug" => \$debug
);

my $warn  = '';
my $error = '';
$genome_bed_file ||= '';
unless ( -f $bam_files[0] ) {
	$error .= "the cmd line switch -bam_files is undefined!\n";
}
unless ( defined $organism and !-f $genome_bed_file ) {
	$error .= "the cmd line switch -organism is undefined!\n";
}
unless ( defined $version and !-f $genome_bed_file ) {
	$error .= "the cmd line switch -version is undefined!\n";
}
unless ( -f $genome_bed_file ) {
	$warn .= "the cmd line switch -genome_bed_file is undefined\n";
}
unless ( defined $outpath ) {
	$error .= "the cmd line switch -outpath is undefined!\n";
}
unless ( defined $proc ) {
	$proc = 3;
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

$task_description .= 'perl ' . $plugin_path . '/quantifyUMIbam.pl';
$task_description .= ' -bam_files "' . join( '" "', @bam_files ) . '"'
  if ( defined $bam_files[0] );
$task_description .= " -organism '$organism'" if ( defined $organism );
$task_description .= " -version '$version'"   if ( defined $version );
$task_description .= " -genome_bed_file '$genome_bed_file'"
  if ( defined $genome_bed_file );
$task_description .= " -outpath '$outpath'" if ( defined $outpath );
$task_description .= " -proc $proc";

mkdir($outpath) unless ( -d $outpath );
open( LOG, ">$outpath/" . $$ . "_quantifyUMIbam.pl.log" ) or die $!;
print LOG $task_description . "\n";
close(LOG);

## Do whatever you want!
## Get the genome bed file
my ($genome_bed);
unless ( -f $genome_bed_file or -f  $outpath . "/genome.bed" ) {
	my $genomeDB = genomeDB->new( variable_table->getDBH() );
	my $interface =
	  $genomeDB->GetDatabaseInterface_for_Organism_and_Version( $organism,
		$version );

	my $run_options = {
		'tag'  => 'mRNA',
		'all_exons' => 1,
	};

	$genome_bed = $interface->get_as_bed_file($run_options);
	$genome_bed = $genome_bed -> sort();
	$genome_bed->write_file( $outpath . "/genome.bed" );
}
elsif ( -f $genome_bed_file ) {
	$genome_bed = stefans_libs_file_readers_bed_file->new();
	$genome_bed->read_file($genome_bed_file);
}else {
	$genome_bed = stefans_libs_file_readers_bed_file->new();
	$genome_bed->read_file( $outpath . "/genome.bed" );
}

my (@line, $fm,@header, @result, $not_matched);


my $pm = Parallel::ForkManager->new($proc);
my $data_table = data_table->new();

FILES:
foreach my $file (@bam_files) {
	$pm->start and next FILES; # do the fork
	
	$fm = root->filemap ( $file );
	my $IN;
	if ( $file =~ m/bam$/ ) {
		open( $IN, "samtools view $file |" )
		  or die "I could not open the bam file '$file'\n$!\n";
	}elsif ( $file =~ m/gz$/) {
		open( $IN, "gunzip -c $file |" )
		  or die "I could not open the gzipped sam file '$file'\n$!\n";
	}else {
		open( $IN, "<$file" )
		  or die "I could not open the sam file $file\n$!\n";
	}
	open ( my $OUT  , ">$outpath/$fm->{'filename_core'}.txt" ) or die "could not create the outfile $!\n";
	print $OUT "GeneSymbol\tuniqe\ttotal\n";
	$not_matched = 0;
	my ($gName, $last_pos, $genes, $lastChr);
	$lastChr = '';
	while (<$IN>) {
		next if ( $_ =~ m/^\@/ );
		@line = split( "\t", $_ );
		unless ( $lastChr eq $line[2] ) {
			#warn  "I dump all gene info: \$genes = ".root->print_perl_var_def($genes ).";\n";
			map { &finalize_gene ($last_pos, $genes, $_, $OUT );} keys %$last_pos ;
			$lastChr = $line[2];
			@result = $genome_bed->efficient_match_chr_position( $line[2], $line[3] );
		}
		if ( defined $result[0] and $result[0] =~ m/^\d+$/ ) {
			if (   @{@{$genome_bed->{'data'}}[$result[0]]}[2] > $line[3] ) {
				#match to THIS entry cool no action required
			}
			elsif ( ! defined @{$genome_bed->{'data'}}[$result[0]+1] ) {
				## make sure we did not have a prolem here
				@result = $genome_bed->efficient_match_chr_position( $line[2], $line[3] );
			}
			elsif ( @{@{$genome_bed->{'data'}}[$result[0]]}[2] < $line[3] and @{@{$genome_bed->{'data'}}[$result[0]+1]}[1] > $line[3]){
				## the read is in between two features
				$not_matched ++;
				next;
			}elsif ( @{@{$genome_bed->{'data'}}[$result[0]+1]}[1] <= $line[3] and  @{@{$genome_bed->{'data'}}[$result[0]+1]}[2] >= $line[3] ) {
				$result[0] ++; ## very quick seleting the next entry
			}else {
				## oops - did we miss one?
				@result = $genome_bed->efficient_match_chr_position( $line[2], $line[3] );
				if ( @result == 0 ) { ## for speed up try to get the closest
					@result = $genome_bed->efficient_match_chr_position( $line[2], $line[3], $line[3], 1e+7 );
					foreach ( @result ) {
						next unless ( defined $_ );
						next if (@{@{$genome_bed->{'data'}}[$_]}[2] < $line[3] );
						@result= ( $_ );
						last;
					}
					$not_matched ++;
					next;
				}
			}
		}else {
			## never matched that...
			next if ( @{@{$genome_bed->{'subset_4_PDL'}->{$line[2]}->{'data'}}[0]}[1] >  $line[3] );
			next if ( @{@{$genome_bed->{'subset_4_PDL'}->{$line[2]}->{'data'}}[($genome_bed->{'subset_4_PDL'}->{$line[2]}->Rows())-1]}[2] < $line[3] );
			@result = $genome_bed->efficient_match_chr_position( $line[2], $line[3] );
		}
		
		if ( @result > 0 and ($line[3] >= @{@{$genome_bed->{'data'}}[$result[0]]}[1] and $line[3] <= @{@{$genome_bed->{'data'}}[$result[0]]}[2] ) ) {
			@header = split(":", $line[0] );
			$gName = @{@{$genome_bed->{'data'}}[$result[0]]}[3];
			$last_pos->{$gName} ||= {};
			$last_pos->{$gName} ->{'pos'} = $line[3];
			unless ( $last_pos->{$gName} ->{$header[$#header]} ){
				$last_pos->{$gName} ->{$header[$#header]} = 1;
				$genes ->{$gName} ||= {'total' => 0, 'unique' => 0};
				$genes ->{$gName}->{'total'} ++;
				$genes ->{$gName}->{'unique'} ++;
			}else {
				$genes ->{$gName}->{'total'} ++;
			}
			foreach ( keys %$last_pos ) {
				if ( $line[3] - $last_pos->{$_}->{'pos'} > 1000000 ) {
					&finalize_gene ($last_pos, $genes, $_, $OUT );
				}
			}
		}
		else {$not_matched ++}
		$lastChr = $line[2];
	}
	close ( $OUT );
	close($IN);
	warn "$not_matched reads not in a transcript for file $file\n";
	$pm->finish;
}

sub finalize_gene{
	my ($last_pos, $genes, $gName, $outfile) = @_;
	print $outfile "$gName\t$genes->{$gName}->{'unique'}\t$genes->{$gName}->{'total'}\n";
	delete($last_pos->{$gName});
	delete( $genes->{$gName});
}

$pm->wait_all_children;

