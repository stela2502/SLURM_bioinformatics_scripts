package stefans_libs::FastqFile;

#use FindBin;
#use lib "$FindBin::Bin/../lib/";
use strict;
use warnings;

use stefans_libs::FastqFile::FastqEntry;

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

sub new {

	my ( $class, $hash ) = @_;

	my ($self);

	$self = { fastq_entry => [], };
	foreach ( keys %{$hash} ) {
		$self->{$_} = $hash->{$_};
	}

	bless $self, $class if ( $class eq "stefans_libs::FastqFile" );

	return $self;

}

sub open_file {
	my ( $self, $fname ) = @_;
	unless ( -f $fname ) {
		Carp::confess("The infile $fname does not exist here!\n");
	}
	my $file;
	if ( $fname =~ m/gz$/ ) {
		open( $file, "zcat $fname |" );
	}
	else {
		open( $file, "<$fname" );
	}
	return $file;
}

=head2 file_filter ( sub function { my ( $fastqFile, $entry) = @_; } ,fname  )

the filter function has to process the fastq file entry by entry.

=head3 example usage

my $fastqObj = stefans_libs::FastqFile->new();

open ( my $OUT, ">some_file" );
$restrict = sub {
	my ( $worker, $entry ) = @_;
	$entry->sequence( substr($entry->sequence(),$start,$length) );
	if ( $asFasta ) {
		print $OUT ">".$entry->name()."\n".$entry->sequence()."\n";
	}else {
		$entry->quality( substr($entry->quality(),$start,$length) );
		$entry->write($OUT);
	}
	$entry->clear();
};

$fastqObj->filter_file( $restrict, <fastq file> );

close ( $OUT );

=cut

sub filter_file {
	my ( $self, $function, $fname ) = @_;
	my $file = $self->open_file($fname);
	my ($line);
	if ( $self->{'debug'} ) {
		my $i = 0;
	  FILTER:
		while (1) {
			$line = $self->read_file_line($file);
			last unless ($line);
			$line = $self->fastq_line( $line, 0 );
			next unless ($line);
			&{$function}( $self, $line );
			if ( $i++ == 1000 ) {
				warn
"FastqFile is in debug mode - filter stopped after 1000 reads\n";
				last FILTER;
			}
		}
	}
	else {
		while (1) {
			$line = $self->read_file_line($file);
			last unless ($line);
			$line = $self->fastq_line( $line, 0 );
			next unless ($line);
			&{$function}( $self, $line );
		}
	}

	close($file);
	return $self;
}

=head2 filter_multiple_files ($function, @files )

this function iterates over several fastq files (no GLOB just the filename!).
It will apply the function $function on n total fastq entries only:
&{$function}( @fastq_entries )


=head3 example usage

taken from Cropmium_Single_Cell_Perl package SplitToCells.pl script:

my $func = sub{
	my ($fastqfile, @entries) = @_;## f1, f2, i1
	
	my $cellid = "S_". $entries[2]->sequence()."_C_". substr($entries[1]->sequence,0,$options->{'cell_barcode_length'});
	my $UMI_tag = substr($entries[1]->sequence,$options->{'cell_barcode_length'});
	$counter->{$cellid} ++;
	$entries[0]->Add_UMI_Tag($cellid."_".$UMI_tag);
	$entries[0]->write( $OUT );	
	
	for ( my $i = 0 ; $i < @entries ; $i++ ) {
		$entries[$i]->clear();
	}
};

my $worker = stefans_libs::FastqFile->new();


$worker->filter_multiple_files(
	$func, $R2, $R1, $I1
);

=cut

sub filter_multiple_files {
	my ( $self, $function, @files ) = @_;
	my @objects = map { $self->open_file($_) } @files;
	my (@lines);
	if ( $self->{'debug'} ) {
		my $i = 0;
	  FILTER:
		while (1) {
			@lines = map { $self->read_file_line($_) } @objects;
			last unless ( $lines[0] );
			for ( my $i = 0 ; $i < @lines ; $i++ ) {
				$lines[$i] = $self->fastq_line( $lines[$i], $i );
			}
			next unless ( $lines[0] );
			&{$function}( $self, @lines );
			if ( $i++ == 1000 ) {
				warn
"FastqFile is in debug mode - filter stopped after 1000 reads\n";
				last FILTER;
			}
		}
	}
	else {
		while (1) {
			@lines = map { $self->read_file_line($_) } @objects;
			last unless ( $lines[0] );
			for ( my $i = 0 ; $i < @lines ; $i++ ) {
				$lines[$i] = $self->fastq_line( $lines[$i], $i );
			}
			next unless ( $lines[0] );
			&{$function}( $self, @lines );
		}
	}

	map { close($_) } @objects;
	return $self;
}

=head2 read_file_line

## code taken from http://stackoverflow.com/questions/2498937/how-can-i-walk-through-two-files-simultaneously-in-perl

=cut

sub read_file_line {
	my ( $self, $fh ) = @_;
	if ( $fh and my $line = <$fh> ) {
		return $line;
	}
	return undef;
}

sub fastq_line {
	my ( $self, $line, $id ) = @_;
	$id ||= 0;
	$self->{fastq_entry} ||= [];
	@{ $self->{fastq_entry} }[$id] ||=
	  stefans_libs::FastqFile::FastqEntry->new();
	@{ $self->{fastq_entry} }[$id]->Add($line);
	if ( @{ $self->{fastq_entry} }[$id]->is_filled() ) {
		return @{ $self->{fastq_entry} }[$id];
	}
	return undef;
}

sub exclude_read {
	my ( $self, $fastqfile, $str, $outfile ) = @_;
	my $OUT;
	open( $OUT, ">$outfile" )
	  or die "I could not create the outfile $outfile\n $!\n";
	my $function = sub {
		my ( $worker, $entry ) = @_;
		if ( $entry->sequence() eq $str ) {
			$worker->{'filtered'}++;
			$entry->clear();
			return;
		}
		$self->{'OK'}++;
		$entry->write($OUT);
		$entry->clear();
	};
	$self->{'OK'} = $self->{'filtered'} = 0;
	$self->filter_file( $function, $fastqfile );
	close($OUT);
	print "$self->{'OK'} OK fastq entries, $self->{'filtered'} filtered.\n";
	return $self;
}

=head2 extract_barcode ( $self, $file, $adapter, $sample_barcode, $startingPos, $outfile )

THIS IS AN EXAMPLE IMPLEMENTATION using a fixed read setup:

Read example
NNN_Barcode_NN_cDNA_Adapter(Which may or may not be present)

re-implement this function if your read example differs from the here used one. 

=cut

sub extract_barcode {
	my ( $self, $file, $adapter, $sample_barcode, $startingPos, $outfile ) = @_;
	my $OUT;
	open( $OUT, ">$outfile" )
	  or die "I could not create the outfile $outfile\n $!\n";
	$startingPos ||= 1;
	my $pattern = "(...)$sample_barcode(..)(.+)";

	#warn $pattern."\n\n";

	sub overlap {
		my ( $str1, $str2 ) = @_;

		# Equalize Lengths
		if ( length $str1 < length $str2 ) {
			$str2 = substr $str2, 0, length($str1);
		}
		elsif ( length $str1 > length $str2 ) {
			$str1 = substr $str1, length($str1) - length($str2);
		}

		# Reduce until match found
		while ( $str1 ne $str2 ) {
			substr $str1, 0, 1, '';
			chop $str2;
		}

		return $str1;
	}
	my $function = sub {
		my ( $worker, $entry ) = @_;
		my $seq;
		my $overlap;
		my $add;
		if ( $entry->sequence() =~ m/$pattern/ and $-[1] < $startingPos ) {
			$self->{'OK'}++;
			$add     = ":$1$2";
			$seq     = $3;
			$overlap = &overlap( $seq, $adapter );
			if ( length($overlap) > 4 ) {
				$seq =~ s/$overlap$//;
				$add = ":AD" . length($overlap) . $add;
			}
			my @tt = split( " ", $entry->name() );
			$entry->name( $tt[0] . "$add " . $tt[1] );
			$entry->sequence() =~ m/$seq/;
			$entry->quality(
				substr( $entry->quality(), $-[0], $+[0] - $-[0] ) );
			$entry->sequence($seq);
			$entry->write($OUT);
			$entry->clear();
		}
		else {
			$self->{'filtered'}++;
			$entry->clear();
		}

	};
	$self->{'OK'} = $self->{'filtered'} = 0;
	$self->filter_file( $function, $file );
	close($OUT);
	print
"$sample_barcode:$self->{'OK'} OK fastq entries, $self->{'filtered'} filtered.\n";
	return $self;
}

sub select_4_str {
	my ( $self, $file, $str, $where, $outfile ) = @_;
	$where ||= 0;
	my $OUT;
	open( $OUT, ">$outfile" )
	  or die "I could not create the outfile $outfile\n $!\n";
	my $function = sub {
		my ( $worker, $entry ) = @_;
		if ( @{ $entry->{'data'} }[$where] =~ m/$str/ ) {
			$worker->{'OK'}++;
			$entry->write($OUT);
			$entry->clear();
		}
		else {
			$worker->{'filtered'}++;
			$entry->clear();
		}
	};
	$self->{'OK'} = $self->{'filtered'} = 0;
	$self->filter_file( $function, $file );
	close($OUT);
	print "$self->{'OK'} OK fastq entries, $self->{'filtered'} filtered.\n";
	return $self;
}

1;
