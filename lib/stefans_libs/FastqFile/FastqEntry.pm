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
	if ( defined $glob ) {
		print $glob join( "\n", @{ $self->{'data'} } ) . "\n";
	}
	else {
		print join( "\n", @{ $self->{'data'} } ) . "\n";
	}
	return $self;
}

sub sequence {
	my ( $self, $v) = @_;
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
	if ( $int ) {
		return  map { $_ - 33 } unpack("C*", @{ $self->{'data'} }[3] );
	}
	return @{ $self->{'data'} }[3];
}

=head3 distance_to (string, start, length)

This function expects the string to be the same length as the internal string.
You will get the sum of all qualities of the mismatched strings + the total missmatches.

Using start and length you can select a substring in the fastq entry.

=cut

sub distance_to {
	my ( $self, $string, $start, $length) = @_;
	my (@values, @qualites, @ref);
	$string = [ $string ] unless ( ref($string) eq "ARRAY" );
	@qualites = $self->quality(undef,1);
	
	if ( defined $start and defined $length) {
		@values = (split( "", substr( $self->sequence(), $start, $length )));
		@qualites =  @qualites[$start..($start+$length-1)] ;
	}elsif ( defined $start ) {
		@values = (split( "", substr( $self->sequence(), $start )));
		@qualites =  @qualites[$start..(@qualites-1)] ;
	}else {
		@values = (split( "",$self->sequence()));
	}
	@ref = split( "",@$string[0]);
	print join("\n",
		join("",@ref ),
		join("", @values),
		join("", @qualites)
	)."\n" if ( $self->{'debug'});
	Carp::confess ( "The arrays do not have the same dimension: ".join(", ", map{scalar(@$_)} \@ref, \@values, \@qualites ) ) 
		unless (scalar(@values) == scalar(@ref) and scalar(@ref)==scalar(@qualites));
	return map { 
		$self->_distance( [split( "",$_)], \@values, \@qualites );
	} @$string;
	
}

sub _distance {
	my ( $self, $str, $ref, $qual ) = @_;
	my $sum = my $mm= 0;
	for( my $i = 0; $i <@$str; $i++) {
		unless ( @$str[$i] eq @$ref[$i]){
			$sum += @$qual[$i];
			$mm ++;
		}
	}
	return { 'dist' =>$sum, 'mm' => $mm} ;
}

sub Add_UMI_Tag {
	my ( $self, $tag ) = @_;
	my @tmp = split(" ", $self->name());
	$tmp[0] .= ":$tag";
	$self->name(join(" ", @tmp ) );
	return $self;
}

sub Get_UMI_Tag {
	my ( $self ) = @_;
	my @tmp = split(" ", $self->name());
	@tmp = split(":", $tmp[0]);
	return pop(@tmp);
}

1;
