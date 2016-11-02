package stefans_libs::file_readers::transfaq;

#use FindBin;
#use lib "$FindBin::Bin/../lib/";
use strict;
use warnings;
use stefans_libs::file_readers::transfaq::motive;
use stefans_libs::flexible_data_structures::data_table;
use stefans_libs::file_readers::bed_file;

=head1 LICENCE

  Copyright (C) 2016-11-01 Stefan Lang

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

stefans_libs::file_readers::transfaq

=head1 DESCRIPTION

Read files in transfaq format like the iCLIP binding site zafira output.

=head2 depends on


=cut

=head1 METHODS

=head2 new ( $hash )

new returns a new object reference of the class stefans_libs::file_readers::transfaq.
All entries of the hash will be copied into the objects hash - be careful t use that right!

=cut

sub new {

	my ( $class, $hash ) = @_;

	my ($self);

	$self = {
		'motives'  => [],
		'rational' => [],
	};
	foreach ( keys %{$hash} ) {
		$self->{$_} = $hash->{$_};
	}

	bless $self, $class if ( $class eq "stefans_libs::file_readers::transfaq" );
	if ( $self->{'filename'} ) {
		$self->read_file( $self->{'filename'} );
	}
	return $self;

}

=head2 read_file ( $filename)

reads the file and creates motives and rational tables.
The motives can be plotted as sequence logos using the R seqLogo package,
the rationals (sequences that built up the motive) can e.g. matched back to the original peaks.

The file can not be written again!

=cut

sub read_file {
	my ( $self, $filename ) = @_;
	if ( -f $filename ) {
		$self->{'filename'} = $filename;
	}
	elsif ( -f $self->{'filename'} ) {
		$filename = $self->{'filename'};
	}
	open( my $in, "<$filename" )
	  or die "I could not open the file $filename\n$!\n";
	my ( @line, $logo );
	while (<$in>) {
		chomp($_);
		@line = split( /;?\s+/, $_ );
		if ( $_ =~ m/^[P\d]\d\s/ ) {

			#print "Interesting line $_\n";

			if ( $line[0] eq "P0" ) {
				if ( defined $logo ) {
					push( @{ $self->{'motives'} }, $logo );
				}
				$logo = stefans_libs::file_readers::transfaq::motive->new();
			}
			else {
				$logo->AddDataset(
					{
						'A' => $line[1],
						'C' => $line[2],
						'G' => $line[3],
						'T' => $line[4]
					}
				);
			}
		}
		elsif ( $_ =~ m/^BS/ ) {
			if ( defined $logo ) {
				push( @{ $self->{'motives'} }, $logo );
				$logo = undef;
			}
			my $here = @{ $self->{'motives'} } - 1;
			unless ( defined @{ $self->{'rational'} }[$here] ) {
				@{ $self->{'rational'} }[$here] = data_table->new();
				@{ $self->{'rational'} }[$here]->Add_2_Header(
					[
						'seq',         'peak_id', 'start', 'length',
						'orientation', 'value'
					]
				);
			}
			$line[2] =~ s/sequence_//;
			@{ $self->{'rational'} }[$here]->AddDataset(
				{
					'seq'         => $line[1],
					'peak_id'     => $line[2],
					'start'       => $line[3],
					'length'      => $line[4],
					'orientation' => $line[6],
					'value'       => $line[7]
				}
			);
		}
	}
	close($in);
	if ( defined $logo ) {
		push( @{ $self->{'motives'} }, $logo );
		$logo = undef;
	}
	return $self;
}

sub create_logo {
	my ( $self, $outpath, $logo_id ) = @_;

	unless ( defined $logo_id ) {
		map { $self->create_logo( $outpath, $_ ) }
		  0 .. ( @{ $self->{'motives'} } - 1 );
		return $self;
	}

	@{ $self->{'motives'} }[$logo_id]->create_logo( $outpath, $logo_id + 1 );

	return $self;
}

=head2 map_motives_to_peaks ( $bedFile, $motiv_id )

map_motives_to_peaks used one or all motives in this object and maps back the internal 'rational' data to the oriiginal peaks.

=cut

sub map_motives_to_peaks {
	my ( $self, $bedFile, $motiv_id ) = @_;
	my $return = stefans_libs::file_readers::bed_file->new();
	my ( $bgd_site, $peak );
	unless ( defined $motiv_id ) {
		my $tmp;
		map {
			$tmp = $self->map_motives_to_peaks( $bedFile, $_ );
			push( @{ $return->{'data'} }, @{ $tmp->{'data'} } );
		} 0 .. ( @{ $self->{'rational'} } - 1 );
	}else {

	for (
		my $bdg_id = 0 ;
		$bdg_id < @{ $self->{'rational'} }[$motiv_id]->Rows() ;
		$bdg_id++
	  )
	{
		$bgd_site =
		  @{ $self->{'rational'} }[$motiv_id]->get_line_asHash($bdg_id);
		$peak = $bedFile->get_line_asHash( $bgd_site->{'peak_id'} - 1 );

		#'chromosome', 'start', 'end', 'value'
		push(
			@{ $return->{'data'} },
			[
				$peak->{'chromosome'},
				$peak->{'start'} + $bgd_site->{'start'},
				$peak->{'start'} + $bgd_site->{'start'} + $bgd_site->{'length'},
				'motive_'.$motiv_id,
				@{ $self->{'motives'} }[$motiv_id]
				  ->score_seq( $bgd_site->{'seq'} )
			]
		);
	}
	}
	return $return->sort();
}

sub print2file {
	Carp::confess("You can not write a transfaq motive file using this class");
}

1;
