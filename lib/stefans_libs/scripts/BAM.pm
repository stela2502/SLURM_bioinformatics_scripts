package stefans_libs::scripts::BAM;

#use FindBin;
#use lib "$FindBin::Bin/../lib/";
use strict;
use warnings;

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


=for comment

This document is in Pod format.  To read this, use a Pod formatter,
like 'perldoc perlpod'.

=head1 NAME

stefans_libs::scripts::BAM

=head1 DESCRIPTION

A lib to create sam/bam specific converion scriplets to be used with any clzuster script creator.

=head2 depends on


=cut


=head1 METHODS

=head2 new ( $hash )

new returns a new object reference of the class stefans_libs::scripts::BAM.
All entries of the hash will be copied into the objects hash - be careful t use that right!

=cut

sub new{

	my ( $class, $hash ) = @_;

	my ( $self );

	$self = {
		'p' => 2, # used processors
		'big_wig_urls' => [],
  	};
  	foreach ( keys %{$hash} ) {
  		$self-> {$_} = $hash->{$_};
  	}

  	bless $self, $class  if ( $class eq "stefans_libs::scripts::BAM" );

  	return $self;

}

=head2 SLURUM_load()

Returns the modules that have to be loaded into the SLURM in order to make this work.

=cut

sub SLURUM_load {
	return ( 'icc/2016.1.150-GCC-4.9.3-2.25', 'impi/5.1.2.150', 'HISAT2/2.0.4',
	'SAMtools/1.3.1-HTSlib-1.3.1', 'BEDTools/2.25.0',
	'ucsc-tools/R2016a' );
}

=head2 convert_sam_2_sorted_bam ( <samfile> )

Creates a set of bash commands to convert the samfile to a sorted bam file using samtools.

Returns the script + the final outfile.

=cut

sub convert_sam_2_sorted_bam {
	my ($self, $fm, $final_path) = @_;

	$fm = $self->fm_check($fm);
	
	if ( defined $final_path ){
		Carp::confess ( "Lib change - instead of copying the file to the final path here call ".ref($self)."->recover_files($final_path, (<files you want to get>))");
	}

 	my $ret_file;
	my $f  = $fm->{'path'} . "/" . $fm->{'filename_core'};
	$ret_file = "$f.sorted.bam";

	return join("\n",
		 "samtools view -Sb  $f.sam | samtools sort -\@ "
		  . ( $self->{p} )
		  . " -o $f.sorted.bam -",
		"if  [ -f $f.sorted.bam ]&&[ -s $f.sorted.bam ]; then", "rm -f $f.sam",
		"fi",""
	), $ret_file;
}

=head3 fm_check ($fm )

This function checks if the file is a root->filemap() result and converts it into one if not.

=cut

sub fm_check {
	my ( $self, $fm ) = @_;
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
=head3 recover_files( $final_path, @files )

This function checks whether the files are already in the final path. If they are not the files get moved there.

=cut

sub recover_files {
	my ( $self, $final_path, @files) = @_;
	my $cmd = "";
	foreach my $fm ( @files ) {
		$fm = $self->fm_check($fm);
		my $f  = $fm->{'path'} . "/" . $fm->{'filename'};
		my $final = $final_path."/".$fm->{'filename'};
		unless ( $final eq $f ){
			$cmd .= "mv $f $final\n";
		}
	}
	return $cmd
}

sub convert_sam_2_bam {
	my ($self, $fm, $final_path) = @_;
	
	$fm = $self->fm_check($fm);
 	
 	if ( defined $final_path ){
		Carp::confess ( "Lib change - instead of copying the file to the final path here call ".ref($self)."->recover_files($final_path, (<files you want to get>))");
	}
	
 	my $ret_file;
	my $f  = $fm->{'path'} . "/" . $fm->{'filename_core'};
	
	$ret_file = "$f.bam";
	
	return join("\n",
		 "samtools view -Sb  $f.sam > $f.bam ",
		"if  [ -f $f.bam ]&&[ -s $f.bam ]; then", "rm -f $f.sam",
		"fi",""
	), $ret_file;
}

=head3 convert_sorted_bam_2_bigwig ( <sorted bam file>, <genome file for bedtools genomecov>)

Creates a set of bash commands to convert the sorted bym file to a bigwig file using bedtools genomecov.

=cut

sub convert_sorted_bam_2_bedGraph {
	my ($self, $fm, $coverage )= @_;
	unless ( ref($fm) eq "HASH") {
 		$fm = root->filemap( $fm );
 	}
	my $infile = $fm->{'total'};
	return "## no -coverage option - no bigwig conversion\n"
	  unless ( -f $coverage );
	my $outfile = $fm->{'path'} . "/". $fm->{'filename_core'}.".bedGraph";
	$outfile =~ s/.sorted.bedGraph$/.bedGraph/;
	return "bedtools genomecov -bg -split -ibam $infile -g $coverage | sort -k1,1 -k2,2n > $outfile", $outfile;
}

sub convert_bedGraph_2_bigwig {
	my ($self, $fm, $coverage )= @_;
	unless ( ref($fm) eq "HASH") {
 		$fm = root->filemap( $fm );
 	}
	my $outfile = $fm->{'path'} . "/". $fm->{'filename_core'}.'.bw';
#	#### remove this part when bedGraphToBigWig is working on aurora #########
#	my $w = "bedGraphToBigWig is not working on aurora !?\n#If it is working remove this warning!";
# 	warn $w."\n";
# 	return "#$w", $outfile;
# 	##########################################################################
	my $infile = $fm->{'total'};
	my $cmd ="polish_bed_like_files.pl -bed_file $infile\nbedGraphToBigWig $infile $coverage $outfile";
	push( @{$self->{'big_wig_urls'}},
"track type=bigWig name=\"$fm->{'filename_core'}\" description=\"$fm->{'filename_core'}\""
		  . " bigDataUrl=http://bone.bmc.lu.se/Public/$fm->{'filename_core'}.bw"
	);
	return $cmd, $outfile;
}

1;
