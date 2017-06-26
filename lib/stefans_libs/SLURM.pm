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

use File::HomeDir;
use File::Spec::Functions;

use stefans_libs::flexible_data_structures::optionsFile;

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

sub new {

	my ( $class, $options ) = @_;

	my $self = {
		'SLURM_modules' => [],
		'sub_SLURMS'    => [],
		'options' => stefans_libs::flexible_data_structures::optionsFile->new({
			'default_file' => File::Spec::Functions::catfile( File::HomeDir->home(), '.SLURM_options.txt' ), 
			'required' => ['n', 'N', 't' ], 
			'optional' => ['A','p', 'mail-user', 'mail-type', 'mem-per-cpu' ], 
		}),
	};

	$self->{'debug'} ||= 0;

	bless $self, $class if ( $class eq "stefans_libs::SLURM" );

	$self->options($options);

	$self->{'max_jobs'} ||= 40;

	open( IN, "whoami |" ) or die $!;
	my @IN = <IN>;
	close(IN);
	my $name = $IN[0];
	chomp($name);
	$self->{'username'} = $name;

	$self->clean_slurm_options($options);
	return $self;
}

sub define_Subscript {
	my ($self) = @_;
	my $sub = ref($self)->new();
	push( @{ $self->{'sub_SLURMS'} }, $sub );
	return $sub;
}

sub options {
	my ( $self, $options ) = @_;
	if ( ref($options) eq "HASH" ) {
		foreach my $key ( keys %$options ) {
			$self->{'options'}->add($key, $options->{$key} );
		}
		unless ( $self->{'options'}->OK()) {
			$self->{'options'}->load();
		}
	}
	
	if ( defined $self->{'options'}->value('mail-user') ) {
		$self->{'options'}->add('mail-type', 'END', 0 );
		my $OK = { map { $_ => 1 } 'BEGIN', 'END', 'FAIL', 'REQUEUE', 'ALL' };
		unless ( $OK->{  $self->{'options'}->value('mail-type') } ) {
			warn
"option 'mail-type' ". $self->{'options'}->value('mail-type') ." is not supported - set to 'END'\n";
			$self->{'options'}->add('mail-type', 'END', 1 );
		}
	}
	return $self;
}

sub wait_for_last_finished {
	my ( $self, $fname ) = @_;
	my $wait;
	while ( $wait = $self->in_pipeline('PD') > 4 ) {
		warn "waiting for $wait processes\n";
		sleep(50);
	}
}

sub pids_finished {
	my ( $self, @pids ) = @_;
	return 1 unless (@pids);
	my $cmd = "squeue -u $self->{'username'} |";
	open( IN, $cmd ) or die $!;
	my @IN = <IN>;
	close(IN);
	my $ret = 1;
	my (@line);
	my $search = { map { $_ => 1 } @pids };

	foreach (@IN) {
		$_ =~ s/^\s*//;
		@line = split( /\s+/, $_ );
		$ret = 0 if ( $search->{ $line[0] } );
	}
	return $ret;
}

sub in_pipeline {
	my $self   = shift;
	my $select = shift;

	#print "squeue -u $self->{'username'}\n";
	open( IN, "squeue -u $self->{'username'} |" ) or die $!;
	my @IN = <IN>;
	close(IN);
	## pending
	## I also need to take care about the different partititions:
	#print "In total ".scalar(@IN)." jobs\n";
	if ( $self->{'partitition'} ) {
		@IN = grep( /\s$self->{'partitition'}\s/, @IN );

		#print scalar(@IN)." jobs in partitcion $self->{'partitition'}\n";
	}
	else {
		@IN = grep( /\ssnic\s/, @IN );
	}
	if ($select) {
		return scalar( grep ( /\s$select\s/, @IN ) );
	}
	## all
	return scalar(@IN) - 1;
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
#SBATCH --mail-user=fred@institute.se
#SBATCH --mail-type=END
#SBATCH -A snic2016-4-13

=cut

sub script {
	my ( $self, $cmd, $name ) = @_;
	&check( { cmd => $cmd, name => $name }, 'cmd', 'name' );
	
	#	$self->{'options'}->{'required'} => ['n', 'N', 't' ],
	#	$self->{'options'}->{'optional'} => ['A','p', 'mail-user', 'mail-type', 'mem-per-cpu' ],

	$self->{'options'}->check();
	
	my $ret = '#! /bin/bash' . "\n";
	
	if ( $self->{'options'}->value('partitition') ) { # a hack if the script already uses p for something else
		$self->{'options'}->add( 'p', $self->{'options'}->value('partitition'), 1);
		$self->{'options'}->drop ( 'partitition' );
	}
	elsif ( $self->{'options'}->value('A') =~ m/^lu/ )
	{      ## add one for each partitition you want to explicitly support
		$self->{'options'}->add( 'p', 'lu',1)
	}
	my $tmp;
	foreach my $option ( @{ $self->{'options'}->{'required'} }, @{ $self->{'options'}->{'optional'} }) {
#		if ($option eq "A" ) {
#			Carp::confess ( "I got the option 'A': ".$self->{'options'}->value($option)."\n");
#		}
		next unless ( defined $self->{'options'}->value($option) );
		$tmp = "-";
		$tmp = '--' if ( length($option) > 1 );
		$ret .= "#SBATCH $tmp$option ".$self->{'options'}->value($option)."\n";

	}
	$ret .= join( "\n",
		"#SBATCH -J $name",
		"#SBATCH -o $name" . "%j.out",
		"#SBATCH -e $name" . "%j.err" );
	if ( @{ $self->{'SLURM_modules'} } ) {
		$ret .= "\n" . $self->load_SLURM_modules();
	}
	if ( ref($cmd) eq "ARRAY" ) {
		$ret .= "\n" . shift(@$cmd);
		for ( my $i = 0 ; $i < @$cmd ; $i++ ) {
			$ret .= @{ $self->{'sub_SLURMS'} }[$i]->subscript( @$cmd[$i] );
		}
	}
	else {
		$ret .= "\n$cmd\n";
	}
	return $ret;
}

=head2 subscript ( $cmd )

Create only the module loads and the cmd parts of the script.
This function is used internally to create scripts that need different sets of modules loaded at different timepoints.

=cut

sub subscript {
	my ( $self, $cmd ) = @_;
	my $ret = '';
	if ( @{ $self->{'SLURM_modules'} } ) {
		$ret .= $self->load_SLURM_modules();
	}
	$ret .= "\n" . $cmd . "\n";
	return $ret;
}

sub load_R_x11 {
	my ($self) = @_;
	my $cmd = $self->load_SLURM_modules( 'ifort/2016.1.150-GCC-4.9.3-2.25',
		'impi/5.1.2.150', 'R/3.2.3-libX11-1.6.3' );
	print "please run\n$cmd" if ( $cmd =~ m/\w/ );
	return $cmd;
}

sub load_SLURM_modules {
	my ( $self, @modules ) = @_;
	my $loaded;
	if ( !$self->{'purge'} ) {
		system("bash -c 'module list 2> /tmp/modulelist.tmp'");
		open( IN, "</tmp/modulelist.tmp" )
		  or die
		  "could not open the tmp module list (/tmp/modulelist.tmp)\n$!\n";

		foreach (<IN>) {
			next if ( $_ =~ m/Currently Loaded Modules/ );
			chomp($_);
			$_ =~ s/\s+\d+\)\s+/;/g;
			map {
				if ( $_ =~ m/\w/ ) { $loaded->{$_} = 1 }
			} split( ";", $_ );
		}
		close(IN);
		unlink("/tmp/modulelist.tmp");
		print "Hope there is some loaded module?: "
		  . root->print_perl_var_def($loaded)
		  if ( $self->{debug} );

	}
	else {
		$loaded = {};
	}

	my $modules_to_load = '';
	if ( @{ $self->{'SLURM_modules'} } ) {
		push( @modules, @{ $self->{'SLURM_modules'} } );
	}
	foreach (@modules) {
		$modules_to_load .= "$_ " unless ( $loaded->{$_} );
	}
	if ( $modules_to_load =~ m/\w/ ) {
		$modules_to_load = "module load $modules_to_load";
	}
	if ( $self->{'purge'} ) {
		return "module purge\n$modules_to_load";
	}
	return $modules_to_load;
}

sub clean_slurm_options {
	my ( $self, $hash ) = @_;
	foreach my $t ( qw( n N t A), 'mail-user', 'mail-type', 'mem-per-cpu' ) {
		delete( $hash->{$t} ) if defined( $hash->{$t} );
	}
}

=head3 run ( $cmd, $fm )

Creates the run script based on the $fm =stefans_libs::root->filemap($file).
The script will be named $fm->{filename_core}.sh and the jobname will be $fm->{filename_core}.
The script will be created in $fm->{path} and run using SBATCH if the debug value is not set.

=cut

sub run {
	my ( $self, $cmd, $fm ) = @_;

	unless ( ref($fm) eq "HASH" ) {
		$fm = root->filemap($fm);
	}
	my $s = $self->script( $cmd, $fm->{'filename_core'} );
	open( OUT, ">$fm->{path}/$fm->{'filename_core'}.sh" )
	  or Carp::confess(
"I can not create the script file '$fm->{path}/$fm->{'filename_core'}.sh'\n$!\n"
	  );
	print OUT $s;
	close(OUT);
	my @ALL = split( "\n", $cmd );
	my @OK = grep( !/^#/, @ALL );
	@OK = grep ( !/^\s*$/, @OK );

	if ( @OK > 0 and !$self->{'debug'} ) {

#	print "test if I am allowed to submitt the job: ($self->{'partitition'},$self->{'max_jobs'}) ".$self-> in_pipeline()." >= $self->{'max_jobs'}?\n";
		if ( $self->in_pipeline() >= $self->{'max_jobs'} ) {
			$self->wait_for_last_finished();
		}
		print "sbatch $fm->{path}/$fm->{'filename_core'}.sh\n";
		system("rm $fm->{'filename_core'}*.err $fm->{'filename_core'}*.out");
		open( PID, "sbatch $fm->{path}/$fm->{'filename_core'}.sh |" );
		my $tmp = join( "", <PID> );
		print $tmp;
		if ( $tmp =~ m/Submitted batch job (\d+)/ ) {
			return $1;
		}
		return 1;
	}
	elsif ( @OK == 0 ) {
		print
"All outfiles present for $fm->{path}/$fm->{'filename_core'}.sh - not run\n";
		return 0;
	}
	else {
		print
		  "sbatch $fm->{path}/$fm->{'filename_core'}.sh # not run (DEBUG)\n";
		return 0;
	}
	return 0;
}

=head3 check_4_outfile( $cmd, @outfiles)

Adds a '#' before the command if any outfile exists.

=cut

sub check_4_outfile {
	my ( $self, $cmd, @outfiles ) = @_;
	Carp::confess("I can not check the outfile for '$cmd' - udefined!\n")
	  unless ( defined $outfiles[0] );
	foreach my $outfile (@outfiles) {
		if ( -f $outfile ) {
			warn "outfile '$outfile' is present - I will not re-create it!\n";
			$cmd = "#$cmd";
			last;
		}
	}
	return $cmd;
}

sub check {
	my ( $self, @require ) = @_;
	my $error = '';
	my $type  = 'downstream program';
	if ( ref($self) eq "stefans_libs::SLURM" ) {
		$type = 'SLURM';
	}
	map {
		unless ( defined $self->{$_} ) { $error .= "MISSING $type option $_\n" }
	} @require;

	Carp::confess($error) if ( $error =~ m/\w/ );
}

sub get_options_from_script {
	my ( $self, $filename ) = @_;
	my $ignore = { 'J' => 1, 'o' => 1, 'e' => 1 };
	unless ( -f $filename ) {
		warn
"I can not read from a not existinf file '$filename' at stefans_libs::SLURM::get_options_from_script\n";
		return $self;
	}
	else {
		open( IN, "<$filename" ) or die $!;
		while (<IN>) {
			if ( $_ =~ m/^#SBATCH \-+(\w+)\s+(.+)\s*\n/ ) {
				$self->{$1} = $2 unless ( $ignore->{$1} );
			}
		}
	}
	return $self;
}

1;
