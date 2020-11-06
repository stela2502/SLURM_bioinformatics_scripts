package root;

use Carp;

#  Copyright (C) 2008 Stefan Lang

#  This program is free software; you can redistribute it
#  and/or modify it under the terms of the GNU General Public License
#  as published by the Free Software Foundation;
#  either version 3 of the License, or (at your option) any later version.

#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
#  See the GNU General Public License for more details.

#  You should have received a copy of the GNU General Public License
#  along with this program; if not, see <http://www.gnu.org/licenses/>.

use strict;
use Date::Simple;
use File::Spec::Functions;
use Cwd;

=for comment

This document is in Pod format.  To read this, use a Pod formatter,
like "perldoc perlpod".

=head1 NAME

stefans_libs::root

=head1 DESCRIPTION

root is a small collection of methods that are used in many libary parts of my programming work.


=cut

=head1 METHODS

=head2 new

The method new returns a root object. No Variables are needed.

=cut

sub new {
	my ($class) = @_;

	my ($self);
	$self = {};

	bless( $self, $class ) if ( $class eq "root" );

	return $self;
}

=head2 Today

Today somply returns a Data::Simple object.

I<return Date::Simple->new();

L<Date::Simple>
    
=cut

sub Today {
	return Date::Simple->new();
}

=head2 Max

Returns the maximum integer value of the array of variables given.

=cut

sub Max {
	my ( $self, @list ) = @_;
	my $max = 0;
	for ( my $i = 0 ; $i < @list ; $i++ ) {
		$max = $list[$i] if ( $max < $list[$i] && defined $list[$i] );
	}

	#    @list = sort numeric @list;
	#    my $max = $list[@list-1];
	return $max;
}

=head2 Min

Returns the minimum integer value of the array of variables given.

=cut

sub Min {
	my ( $self, @list ) = @_;
	my $min = 1e999;
	for ( my $i = 0 ; $i < @list ; $i++ ) {
		$min = $list[$i] if ( $min > $list[$i] && defined $list[$i] );
	}

	#    @list = sort numeric @list;
	#    my $min = $list[0];
	return $min;
}



sub relative_path {
	my ($self, $fm1, $fm2) = @_;
	## fm1 is the outfile
	## fm2 is the @files entry
	## calc the difference fm1 -> fm2
	my ( @outfile, @infile );
	@outfile = split( "/", $fm1->{'path'});
	@infile = split( "/", $fm2->{'path'} );
	my @ret;
	for ( my $i = 0; $i < @outfile; $i ++ ) {
		unless ( $outfile[$i] eq $infile[$i] ) {
			@ret = ( map { '..' } 1..(@outfile-$i));
			push ( @ret , @infile[$i..$#infile] );
			return join("/",@ret );
		}
	}
	return "./";
}

sub print_perl_var_def {
	my ( $self, $var ) = @_;
	my $return = '';
	if ( ref($var) eq "HASH" ) {
		$return = "{\n";
		foreach my $name ( sort keys %$var ) {
			$return .= "  '$name' => "
			  . $self->print_perl_var_def( $var->{$name} ) . ",\n";
		}
		chop($return);
		chop($return);
		$return .= "\n}";
	}
	elsif ( ref($var) eq "ARRAY" ) {
		$return = "[ ";
		foreach (@$var) {
			$return .= $self->print_perl_var_def($_) . ", ";
		}
		chop($return);
		chop($return);
		$return .= " ]";
	}
	else {
		unless ( defined $var and length($var) > 0 ){
			$return .= 'undef';
		}
		else {
			$var =~ s/'/\\'/g;
			$return .= "'$var'";
		}
	}
	return $return;
}

sub CreatePath {
	my ( $self, $path ) = @_;
	if ( defined $path ) {
		my ( @temp, @path );
		@temp = split( "/", $path );

		for ( my $i = 0 ; $i < @temp ; $i++ ) {
			$path[$i] = $temp[$i];
			mkdir( join( "/", @path ) ) unless ( -d join( "/", @path ) );
		}
		return 1;
	}
	return $path;
}

sub numeric {
	return $a <=> $b;
}

sub print_hashEntries {
	my ( $hash, $maxDepth, $topMessage ) = @_;
	my $p = &get_hashEntries_as_string( $hash, $maxDepth, $topMessage );
	print $p;
	return $p;
}

sub Latex_Label {
	my ( $self, $str ) = @_;
	$str =~ s/_//g;
	$str =~ s/\\//g;
	return $str;
}

sub Print_hashEntries {
	return &get_hashEntries_as_string( @_ );
}
sub get_hashEntries_as_string {
	my ( $hash, $maxDepth, $topMessage ) = @_;

	#	print "output of hashes deactivated in root::print_hashEntries()\n";
	#	return 1;
	my $string = '';
	if ( defined $topMessage ) {
		$string .= "$topMessage\n";
	}
	else {
		$string .= "DEBUG entries of the data structure $hash:\n";

	}

	#warn "production state - no further info\n";
	#return 1;
	if ( $hash =~ m/ARRAY/ ) {
		my $i = 0;
		foreach my $value (@$hash) {
			$string .=
			  printEntry( "List entry $i", $value, 1, $maxDepth, $string );
			$i++;
		}
	}
	elsif ( $hash =~ m/HASH/ ) {
		my $key;
		foreach $key ( sort keys %$hash ) {
			$string .= printEntry( $key, $hash->{$key}, 1, $maxDepth, $string );
		}
	}

	return $string;
}

sub perl_include {
	return
" -I /storage/www/Genexpress/lib/ -I /home/stefan/LibsNewStructure/lib/ ";
}

sub printEntry {
	my ( $key, $value, $i, $maxDepth ) = @_;

	my $max    = 10;
	my $string = '';
	my ( $printableString, $maxStrLength );
	$maxStrLength = 50;

	if ( defined $value ) {
		for ( $a = $i ; $a > 0 ; $a-- ) {
			$string .= "\t";
		}
		$printableString = $value;
		
		if ( length($value) > $maxStrLength ) {
			$printableString = substr( $value, 0, $maxStrLength );
			$printableString = "$printableString ...";
		}
		$printableString =~ s/'/\\'/g;
		$string .= "$key\t$printableString\n";
	}
	else {
		for ( $a = $i ; $a > 0 ; $a-- ) {
			$string .= "\t";
		}
		$printableString = $key;
		if ( length($printableString) > $maxStrLength ) {
			$printableString = substr( $key, 0, $maxStrLength );
			$printableString = "$printableString ...";
		}
		$printableString =~ s/'/\\'/g;
		$string .= "$printableString\n";
	}
	return $string if ( $maxDepth == $i );
	if ( defined $value ) {
		if ( ref($value) eq "ARRAY" ) {
			$max = 20;
			foreach my $value1 (@$value) {
				$string .=
				  printEntry( $value1, undef, $i + 1, $maxDepth, $string )."\n"
				  if ( defined $value1 );
				last if ( $max-- == 0 );
			}
		}
		elsif (  $value =~ m/HASH/ ) {
			$max = 20;
			while ( my ( $key1, $value1 ) = each %$value ) {
				$string .=
				  printEntry( $key1, $value1, $i + 1, $maxDepth, $string );
				last if ( $max-- == 0 );
			}
		}
	}
	if ( defined $key ) {
		if ( ref($key) eq "ARRAY" ) {
			$max = 20;
			foreach my $value1 (@$key) {
				$string .=
				  printEntry( $value1, undef, $i + 1, $maxDepth, $string );
				last if ( $max-- == 0 );
			}
		}
		elsif ( ref($key) eq "HASH" ) {
			$max = 20;
			while ( my ( $key1, $value1 ) = each %$key ) {
				$string .=
				  printEntry( $key1, $value1, $i + 1, $maxDepth, $string );
				last if ( $max-- == 0 );
			}
		}
	}
	return $string;
}

1;
