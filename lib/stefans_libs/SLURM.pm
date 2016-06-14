package stefans_libs::SLURM;
#  Copyright (C) 2016-05-31 Stefan Lang

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

#use FindBin;
#use lib "$FindBin::Bin/../lib/";
use strict;
use warnings;

use stefans_libs::root;

=for comment

This document is in Pod format.  To read this, use a Pod formatter,
like 'perldoc perlpod'.

=head1 NAME

stefans_libs::SLURM

=head1 DESCRIPTION

The lib helps to create a SLURM batch script

=head2 depends on

=head2 INTERESTING

get q qlogin shell using SLURM:
1 node 10 cores for 2h max
srun -N 1 -n 10 t 02:00:00 --pty /bin/bash


=head1 METHODS

=head2 new

new returns a new object reference of the class stefans_libs::SLURM.

=cut

sub new{

	my ( $class, $options ) = @_;

	my $self = {
	};
	foreach my $key  (keys %$options ) {
		$self->{$key} = $options->{$key};
	}
	$self->{'debug'} ||= 0;
	

  	bless $self, $class  if ( $class eq "stefans_libs::SLURM" );

  	return $self;

}


=head3 script ( $cmd )
Creates a script file string like that:

#! /bin/bash
#SBATCH -n 2
#SBATCH -N 1
#SBATCH -t 01:00:00
#SBATCH -J '.$name'
#SBATCH -o '.$name.'_omp_%j.out'
#SBATCH -e '.$name.'_omp_%j.err'

=cut

sub script {
	my ( $self, $cmd, $name ) = @_;
	&check ( { cmd=>$cmd, name=> $name}, 'cmd', 'name' );
	my @o = qw( n N t);
	$self->check( @o );
	my $ret = '#! /bin/bash'."\n";
	foreach my $option ( @o ) {
		if ( length( $option) == 1 ) {
			$ret .= "#SBATCH -$option $self->{$option}\n";
		}else {
			$ret .= "#SBATCH --$option $self->{$option}\n";
		}
		
	}
	$ret .= join("\n", "#SBATCH -J $name","#SBATCH -o $name"."%j.out","#SBATCH -e $name"."%j.err");
	$ret .= "\n$cmd\n";
	return $ret;
}

=head3 run ( $cmd, $fm )

Creates the run script based on the $fm =stefans_libs::root->filemap($file).
The script will be named $fm->{filename_core}.sh and the jobname will be $fm->{filename_core}.
The script will be created in $fm->{path} and run using SBATCH if the debug value is not set.

=cut

sub run {
	my ( $self, $cmd, $fm, @targets ) = @_;
	my $s = $self->script($cmd, $fm->{'filename_core'} );
	open ( OUT ,">$fm->{path}/$fm->{'filename_core'}.sh" ) or Carp::confess ( "I can not create the script file '$fm->{path}/$fm->{'filename_core'}.sh'\n$!\n");
	print OUT $s;
	close ( OUT );
	print "sbatch $fm->{path}/$fm->{'filename_core'}.sh\n";
	my @ALL = split("\n", $cmd);
	my @OK = grep( ! /^#/, @ALL );
	@OK = grep ( ! /^\s*$/, @OK );
	unless ( $self->{'debug'}) {
		if ( @OK > 0 ) {
			system( "sbatch $fm->{path}/$fm->{'filename_core'}.sh" )
		}else {
			print "Useless to run script - no commands in the file\n";
		}
	}
	return 1;
}
=head3 check_4_outfile( $cmd, $outfile)

Adds a '#' before the command if the outfile exists.

=cut

sub check_4_outfile {
	my ( $self, $cmd, $outfile ) =@_;
	$cmd = "#$cmd" if ( -f $outfile );
	return $cmd;
}

sub check {
	my ( $self, @require ) =@_;
	my $error = '';
	my $type = 'downstream program';
	if ( ref( $self ) eq "stefans_libs::SLURM"){
		$type = 'SLURM';
	}
	map { unless (defined $self->{$_}){ $error .= "MISSING $type option $_\n" } } @require;
	Carp::confess ( $error ) if ( $error =~ m/\w/ );
}

sub get_options_from_script{
	my ( $self, $filename ) = @_;
	my $ignore = { 'J'=> 1, 'o' =>1, 'e' =>1 };
	unless ( -f $filename ) {
		warn "I can not read from a not existinf file '$filename' at stefans_libs::SLURM::get_options_from_script\n";
		return $self;
	}
	else {
		open ( IN, "<$filename" ) or die $!;
		while ( <IN> ) {
			if ( $_ =~m/^#SBATCH \-+(\w+)\s+(.+)\s*\n/ ) {
				$self->{$1} = $2 unless ( $ignore->{$1} );
			}
		}
	}
	return $self;
}

1;
