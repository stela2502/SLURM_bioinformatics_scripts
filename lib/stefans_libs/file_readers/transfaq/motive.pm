package stefans_libs::file_readers::transfaq::motive;

#  Copyright (C) 2016-11-01 Stefan Lang 

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
use warnings;

use stefans_libs::flexible_data_structures::data_table;
use base 'data_table';

=head1 General description

A file interface for the transfay motive

=cut
sub new {

    my ( $class, $debug ) = @_;
    my ($self);
    $self = {
        'debug'           => $debug,
        'arraySorter'     => arraySorter->new(),
        'header_position' => { 
            'A' => 0,
            'C' => 1,
            'G' => 2,
            'T' => 3,
        },
        'default_value'   => [],
        'header'          => [
            'A',
            'C',
            'G',
            'T',       ],
        'data'            => [],
        'index'           => {},
        'last_warning'    => '',
        'subsets'         => {}
    };
    bless $self, $class if ( $class eq "stefans_libs::file_readers::transfaq::motive" );

    return $self;
}


## two function you can use to modify the reading of the data.

sub pre_process_array{
	my ( $self, $data ) = @_;
	##you could remove some header entries, that are not really tagged as such...
	return 1;
}

sub After_Data_read {
	my ($self) = @_;
	return 1;
}



sub Add_2_Header {
    my ( $self, $value ) = @_;
    return undef unless ( defined $value );
    unless ( defined $self->{'header_position'}->{$value} ) {
        Carp::confess( "You try to change the table structure - That is not allowed!\n".
            "If you really want to change the data please use ".
            "the original data_table class to modify the table structure!\n"
        ) ;
    }
    return $self->{'header_position'}->{$value};
}

sub score_seq {
	my ( $self, $seq ) = @_;
	if ( length($seq) != $self->Rows() ) {
		warn "The sequence does not have the right length (me:". $self->Rows() ." != seq:".length($seq)."\n";
		return undef; 
	}
	my $ret = 0;
	my $i = 0;
	foreach ( split( "", uc($seq) ) ) {
		$ret += @{@{$self->{data}}[$i++]}[ $self->Header_Position( $_ ) ];
	}
	return int($ret *1000 / $self->max() + 0.99);
}

sub max { 
	my $self = shift;
	return $self->{'max'} if ( defined $self->{'max'} );
	my $lmax;
	for ( my $i = 0; $i < $self->Rows(); $i ++ ) {
		$lmax = 0;
		map { $lmax = $_ if ( $_ > $lmax) } @{@{$self->{data}}[$i]};
		$self->{'max'} += $lmax;
	}
	return $self->{'max'};
}

sub create_logo {
	my ( $self, $outpath, $logo_id ) = @_;
	
	for ( my $i = 0 ; $i < $self->Lines() ; $i++ ) {
		&div( @{ $self->{'data'} }[$i] );
	}
	my $sum = &sum( @{ $self->{'data'} }[0] );
	
	open( OUTR, ">$outpath/$logo_id" . "_$sum.R" ) or die $!;
	print OUTR "library(seqLogo)\n"
	  . "proportion <- function(x){\n"
	  . "   rs <- sum(x);\n"
	  . "   return(x / rs);\n" . "}\n";
	foreach my $l ( @{ $self->{'header'} } ) {
		print OUTR "$l <- c("
		  . join( ", ", @{$self->GetAsArray($l) } ) . ")\n";
	}
	print OUTR "df <- data.frame( "
	  . join( ",", @{ $self->{'header'} } ) . " )\n"
	  . "pwm <- apply(df, 1, proportion)\n"
	  . "pwm <- makePWM(pwm)\n"
	  . "pdf( file = '$outpath/$logo_id"
	  . "_$sum.pdf', width=6, height=6)\n"
	  . "seqLogo(pwm)\n"
	  . "dev.off()\n";
	close(OUTR);

	system(
		"R CMD BATCH --no-save --no-restore --no-readline -- $outpath/$logo_id"
		  . "_$sum.R" );

	return $self;
}
sub div {
	my @val = @{$_[0]};
	my $sum = &sum($_[0]);
	map { $val[$_] = $val[$_]/$sum } 0..(@val-1);
}
sub sum {
	my @val = @{$_[0]};
	my $ret = 0; 
	map { $ret += $_ } @val;
	return $ret;
}


1;
