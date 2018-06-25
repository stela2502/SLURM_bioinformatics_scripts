package stefans_libs::FastqFile::FastqEntry;

#use FindBin;
#use lib "$FindBin::Bin/../lib/";
use strict;
use warnings;

=head1 LICENCE

  Copyright (C) 2017-02-17 Stefan Lang

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

stefans_libs::FastqFile::FastqEntry

=head1 DESCRIPTION

A simple fastq entry class

=head2 depends on


=cut

=head1 METHODS

=head2 new ( $hash )

new returns a new object reference of the class stefans_libs::FastqFile::FastqEntry.
All entries of the hash will be copied into the objects hash - be careful t use that right!

=cut

sub new {

	my ( $class, $hash ) = @_;

	my ($self);

	$self = {};
	foreach ( keys %{$hash} ) {
		$self->{$_} = $hash->{$_};
	}

	bless $self, $class if ( $class eq "stefans_libs::FastqFile::FastqEntry" );
	$self->clear();

	return $self;

}

sub clear {
	my ($self) = @_;
	$self->{'data'} = [];
	return $self;
}

sub is_filled {
	my $self = shift;
	return @{ $self->{'data'} } == 4;
}

sub Add {
	my ( $self, $line ) = @_;
	if ( defined $line ) {
		chomp($line);
		unless ( $self->is_filled() ) {
			push( @{ $self->{'data'} }, $line );
		}
		else {
			Carp::confess("You trie to over feed a fastq entry!");
		}
	}
	return $self;
}

sub copy {
	my ($self) = @_;
	my $ret = stefans_libs::FastqFile::FastqEntry->new();
	foreach ( @{ $self->{'data'} } ) {
		$ret->Add($_);
	}
	return $ret;
}

sub write {
	my ( $self, $glob ) = @_;
	@{ $self->{'data'} }[2] = '+';
	if ( defined $glob ) {
		print $glob join( "\n", @{ $self->{'data'} } ) . "\n";
	}
	else {
		print join( "\n", @{ $self->{'data'} } ) . "\n";
	}
	return $self;
}

sub sequence {
	my ( $self, $v ) = @_;
	@{ $self->{'data'} }[1] = $v if ( defined $v );
	return @{ $self->{'data'} }[1];
}

sub name {
	my ( $self, $v ) = @_;
	@{ $self->{'data'} }[0] = $v if ( defined $v );
	return @{ $self->{'data'} }[0];
}

sub quality {
	my ( $self, $v, $int ) = @_;
	$int ||= 0;
	@{ $self->{'data'} }[3] = $v if ( defined $v );
	if ($int) {
		return map { $_ - 33 } unpack( "C*", @{ $self->{'data'} }[3] );
	}
	return @{ $self->{'data'} }[3];
}

=head is_polyA ($ignore, $cutoff)

Checks if more than 95% of the read -$ignore (||=10) first bases are A's.
It also removes bad quality bases to before checking the A contend ($cutoff = 10).
Does not tuch the read itself.
Return 1 or 0.

=cut

sub is_polyA {
	my ( $self, $ignore, $cutoff ) = @_;
	$ignore ||= 10;
	$cutoff ||= 10;
	$self = $self->copy();
	$self->drop_low_quality($cutoff);
	return 1 if ( length( $self->sequence ) == 0 );
	$self->trim( 'start', 10 );
	return 1 if ( length( $self->sequence ) == 0 );
	my $notA = 0;

	map { $notA++ unless ( $_ eq "A" ) } split( "", $self->sequence() );
	if ( $notA / length( $self->sequence() ) <= 0.05 )
	{    ## less than 5 % not A's - likely a polyA only read
		return 1;
	}
	return 0;
}

=head3 trim ($self, $where, $position)

trim up to bp $position if $where = 'start' and from $position to end if $where = 'end'

=cut

sub trim {
	my ( $self, $where, $position ) = @_;
	return $self unless ($position);    # 0 or undefined
	my $length = length( $self->sequence() );

	if ( $where eq "start" ) {
		$self->sequence(
			substr( $self->sequence(), $position, $length - $position ) );
		$self->quality(
			substr( $self->quality(), $position, $length - $position ) );
	}
	elsif ( $where eq "end" ) {
		$self->sequence( substr( $self->sequence(), 0, $length - $position ) );
		$self->quality( substr( $self->quality(), 0, $length - $position ) )
		  ;
	}
	else {
		Carp::confess(
"trim ($self, $where, $position) can not trim the position '$where', only start or end"
		);
	}

	return $self;
}

=head filter_low_quality ($cutoff)

This will filter areas at the beginning and end of the read where the mean quality is below $cutoff (||=10).

=cut

sub filter_low_quality {
	my ( $self, $cutoff ) = @_;
	$cutoff ||= 10;
	my @seq   = split( "", $self->sequence() );
	my @origQ = split( "", $self->quality() );
	my @qual = $self->quality( undef, 1 );
	my ( $sum, $n );
	## trim start
	$sum = $n = 0;
	for ( my $i = 0 ; $i < @seq ; $i++ ) {
		$sum += $qual[$i];
		$n++;
		if ( $sum / $n > $cutoff ) {
			$n--;
			last;
		}
	}
	if ( $n > 0 ) {
		splice( @seq,   0, $n );
		splice( @origQ, 0, $n );
		splice( @qual,  0, $n );
	}

	# trim end
	$sum = $n = 0;
	for ( my $i = @seq - 1 ; $i >= 0 ; $i-- ) {
		$sum += $qual[$i];
		$n++;
		if ( $sum / $n > $cutoff ) {
			$n--;
			last;
		}
	}
	if ( $n > 0 ) {
		my $s = @seq - $n;
		splice( @seq,   $s, $n );
		splice( @origQ, $s, $n );
		splice( @qual,  $s, $n );
	}
	$self->sequence( join( "", @seq ) );
	$self->quality( join( "", @origQ ) );

	#warn "seq sength after filter:".scalar(@seq)."\n";
	return $self;
}

=head drop_low_quality ($cutoff)

will drop all sequence entries with a quality below $cutoff (||=10).
THIS WILL ADD INDELS TO THE SEQUENCE

=cut

sub drop_low_quality {
	my ( $self, $cutoff ) = @_;
	$cutoff ||= 10;
	my @seq   = split( "", $self->sequence() );
	my @origQ = split( "", $self->quality() );

	my @qual = $self->quality( undef, 1 );
	my $c = 0;
	for ( my $i = @qual - 1 ; $i >= 0 ; $i-- ) {
		if ( $qual[$i] <= $cutoff ) {
			$c = 1;
			splice( @seq,   $i, 1 );
			splice( @origQ, $i, 1 );
		}
	}
	if ($c) {
		$self->quality( join( "", @origQ ) );
		$self->sequence( join( "", @seq ) );
	}
	return $self;
}

=head3 distance_to (string, start, length)

This function expects the string to be the same length as the internal string.
You will get the sum of all qualities of the mismatched strings + the total missmatches.

Using start and length you can select a substring in the fastq entry.

=cut

sub distance_to {
	my ( $self, $string, $start, $length ) = @_;
	my ( @values, @qualites, @ref );
	$string = [$string] unless ( ref($string) eq "ARRAY" );
	@qualites = $self->quality( undef, 1 );

	if ( defined $start and defined $length ) {
		@values = ( split( "", substr( $self->sequence(), $start, $length ) ) );
		@qualites = @qualites[ $start .. ( $start + $length - 1 ) ];
	}
	elsif ( defined $start ) {
		@values = ( split( "", substr( $self->sequence(), $start ) ) );
		@qualites = @qualites[ $start .. ( @qualites - 1 ) ];
	}
	else {
		@values = ( split( "", $self->sequence() ) );
	}
	@ref = split( "", @$string[0] );
	print
	  join( "\n", join( "", @ref ), join( "", @values ), join( "", @qualites ) )
	  . "\n"
	  if ( $self->{'debug'} );
	Carp::confess( "The arrays do not have the same dimension: "
		  . join( ", ", map { scalar(@$_) } \@ref, \@values, \@qualites ) )
	  unless ( scalar(@values) == scalar(@ref)
		and scalar(@ref) == scalar(@qualites) );
	return
	  map { $self->_distance( [ split( "", $_ ) ], \@values, \@qualites ); }
	  @$string;

}

sub _distance {
	my ( $self, $str, $ref, $qual ) = @_;
	my $sum = my $mm = 0;
	for ( my $i = 0 ; $i < @$str ; $i++ ) {
		unless ( @$str[$i] eq @$ref[$i] ) {
			$sum += @$qual[$i];
			$mm++;
		}
	}
	return { 'dist' => $sum, 'mm' => $mm };
}

sub Add_UMI_Tag {
	my ( $self, $tag ) = @_;
	my @tmp = split( " ", $self->name() );
	$tmp[0] .= ":$tag";
	$self->name( join( " ", @tmp ) );
	return $self;
}

sub Get_UMI_Tag {
	my ($self) = @_;
	my @tmp = split( " ", $self->name() );
	@tmp = split( ":", $tmp[0] );
	return pop(@tmp);
}

1;
