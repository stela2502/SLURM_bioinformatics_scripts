package stefans_libs::FastqFile;

#use FindBin;
#use lib "$FindBin::Bin/../lib/";
use strict;
use warnings;

=head1 LICENCE

  Copyright (C) 2016-10-11 Stefan Lang

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

stefans_libs::FastqFile

=head1 DESCRIPTION

A simple interface to a fastq file

=head2 depends on


=cut


=head1 METHODS

=head2 new ( $hash )

new returns a new object reference of the class stefans_libs::FastqFile.
All entries of the hash will be copied into the objects hash - be careful t use that right!

=cut

sub new{

	my ( $class, $hash ) = @_;

	my ( $self );

	$self = {
		fastq_entry => [],
  	};
  	foreach ( keys %{$hash} ) {
  		$self-> {$_} = $hash->{$_};
  	}

  	bless $self, $class  if ( $class eq "stefans_libs::FastqFile" );

  	return $self;

}

sub open_file{
	my ( $self, $fname ) = @_;
	unless ( -f $fname ) {
		Carp::confess ( "The infile $fname does not exist here!\n" );
	}
	my $file;
	if ( $fname =~ m/gz$/ ) {
		open ( $file , "zcat $fname |") 
	}else {
		open ( $file , "<$fname") 
	}
	return $file;
}

=head2 file_filter ( fname, sub funciton {}  )

the filter function has to process the fastq file all by itself!
It is called like that:

while ( <$file> ) {
	chomp($_);
	&{$filter}( $self, $_ );
}

=cut

sub filter_file {
	my ( $self, $fname, $filter ) = @_;
	my $file = $self->open_file($fname);
	while ( <$file> ) {
		chomp($_);
		&{$filter}( $self, $_ );
	}
	close ( $file);
	return $self;
}

sub fastq_line{
	my ( $self, $line) = @_;
	push ( @{$self->{fastq_entry}}, $line );
	if ( @{$self->{fastq_entry}} == 4 ){
		my $ret =[@{$self->{fastq_entry}}];
		$self->{fastq_entry} = [];
		return $ret;
	}
	return undef;
}

sub exclude_read {
	my ( $self, $fastqfile, $str, $outfile ) = @_;
	my $OUT;
	open ( $OUT, ">$outfile" ) or die "I could not create the outfile $outfile\n $!\n";
	my $function = sub {
		my $tmp = $self->fastq_line($_);
		if ( defined $tmp ) {
			if ( @$tmp[1] eq $str){
				$self->{'filtered'} ++;
				return;
			}
			$self->{'OK'} ++;
			print $OUT join("\n", @$tmp)."\n";
		}
	};
	$self->{'OK'} = $self->{'filtered'} = 0;
	$self->filter_file($fastqfile, $function );
	close ( $OUT );
	print "$self->{'OK'} OK fastq entries, $self->{'filtered'} filtered.\n";
	return $self;
}

=head2 extract_barcode ( $self, $file, $adapter, $sample_barcode, $startingPos, $outfile )

THIS IS AN EXAMPLE IMPLEMENTATION using a fixed read setup:

Read example
NNN_Barcode_NN_cDNA_Adapter(Which may or may not be present)

re-implement this function if your read example differs from the here used one. 

=cut

sub extract_barcode{
	my ( $self, $file, $adapter, $sample_barcode, $startingPos, $outfile ) = @_;
	my $OUT;
	open ( $OUT, ">$outfile" ) or die "I could not create the outfile $outfile\n $!\n";
	$startingPos ||= 1;
	my $pattern = "(...)$sample_barcode(..)(.+)";

	#warn $pattern."\n\n";
	
	sub overlap {
    my ($str1, $str2) = @_;

    # Equalize Lengths
    if (length $str1 < length $str2) {
        $str2 = substr $str2, 0, length($str1);
    } elsif (length $str1 > length $str2) {
        $str1 = substr $str1, length($str1) - length($str2);
    }

    # Reduce until match found
    while ($str1 ne $str2) {
        substr $str1, 0, 1, '';
        chop $str2;
    }

   	 return $str1;
	}
	my $function = sub {
		my $tmp = $self->fastq_line($_);
		my $seq;
		my $overlap;
		if ( defined $tmp ) {
			if ( @$tmp[1] =~ m/$pattern/ and  $-[1] < $startingPos ){
				$self->{'OK'} ++;
				my @tt = split(" ", @$tmp[0]);
				@$tmp[0] =$tt[0].":$1$2 ".$tt[1];
				$seq = $3;
				$overlap = &overlap( $seq, $adapter );
				if ( length($overlap) > 4 ) {
					$seq =~ s/$overlap$// if ( length($overlap) > 0 );
				}
				@$tmp[1] =~ m/$seq/;
				@$tmp[3] = substr( @$tmp[3], $-[0], $+[0]- $-[0] );
				@$tmp[1] = $seq;
				print $OUT join("\n", @$tmp)."\n";
			}else {
				$self->{'filtered'} ++;
			}
		}
	};
	$self->{'OK'} = $self->{'filtered'} = 0;
	$self->filter_file($file, $function );
	close ( $OUT );
	print "$sample_barcode:$self->{'OK'} OK fastq entries, $self->{'filtered'} filtered.\n";
	return $self;
}



1;
