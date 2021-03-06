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

use Digest::MD5 qw(md5_hex);
use File::HomeDir;
use File::Spec::Functions;
use stefans_libs::database::slurmscripts;

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

	my ( $class, $options, $clean ) = @_;

	$clean = 0 unless ( defined $clean );

	my $debug = 0;
	
	my $self = {
		'local'     => 0,
		'dbfile'    => $options->{dbfile} || "~/workload.sqlite3",
		'shell'         => 'bash',
		'SLURM_modules' => [],
		'sub_SLURMS'    => [],
		'options' => stefans_libs::flexible_data_structures::optionsFile->new(
			{
				'default_file' => File::Spec::Functions::catfile(
					File::HomeDir->home(), '.SLURM_options.txt'
				),
				'required' => [ 'n', 'N', 't' ],
				'optional' =>
				  [ 'A', 'p', 'mail-user', 'mail-type', 'mem-per-cpu', 'w', 'begin' ],
			},
		'slurmscripts' => stefans_libs::database::slurmscripts->new(undef, $debug), ## creates a file ~/.slurmscripts.sqlite3
		),
	};

	$self->{'debug'} ||= 0;

	bless $self, $class if ( $class eq "stefans_libs::SLURM" );

	$self->options($options);

	$self->{'max_jobs'} ||= 80;

	open( IN, "whoami |" ) or die $!;
	my @IN = <IN>;
	close(IN);
	my $name = $IN[0];
	chomp($name);
	$self->{'username'} = $name;
	if ($clean) {
		$self->clean_slurm_options($options);
	}
	$self->{'sleep_sec'} ||= 10;  ## per default sleep 10 sec between any check.
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
			$self->{'options'}->add( $key, $options->{$key} );
		}
		unless ( $self->{'options'}->OK() ) {
			$self->{'options'}->load();
		}
	}

	if ( defined $self->{'options'}->value('mail-user') ) {
		$self->{'options'}->add( 'mail-type', 'END', 0 );
		my $OK = { map { $_ => 1 } 'BEGIN', 'END', 'FAIL', 'REQUEUE', 'ALL' };
		unless ( $OK->{ $self->{'options'}->value('mail-type') } ) {
			warn "option 'mail-type' "
			  . $self->{'options'}->value('mail-type')
			  . " is not supported - set to 'END'\n";
			$self->{'options'}->add( 'mail-type', 'END', 1 );
		}
	}
	return $self;
}

sub wait_for_last_finished {
	my ( $self, $fname ) = @_;
	my $wait;
	while ( $wait = $self->in_pipeline('PD') > 4 ) {
		warn "waiting for $wait processes\n";
		sleep( $self->{'sleep_sec'} );
	}
}

sub pids_finished {
	my ( $self, @pids ) = @_;
	return 1 if ( $self->{'local'} );
	return 1 unless (@pids);
#	my $cmd = "squeue -u $self->{'username'} |";
#	open( IN, $cmd ) or die $!;
#	my @IN = <IN>;
#	close(IN);
#	my $ret = 1;
#	my (@line);
#	my $search = { map { $_ => 1 } @pids };
#
#	foreach (@IN) {
#		$_ =~ s/^\s*//;
#		@line = split( /\s+/, $_ );
#		$ret = 0 if ( $search->{ $line[0] } );
#	}
	## changed to database usage and external help!
	return $self->{'slurmscripts'}->get_data_table_4_search(
		{
			'search_columns' =>
			  [ ref($self) . '.PID' ],
			'where' => [ [ ref($self) . '.Finished', '=', 'my_value' ], [ ref($self).'.PID', '=', 'my_value'] ],
		},
		0, \@pids
	)->GetAsHash('PID');
	#return $ret;
}

sub in_pipeline {
	my $self   = shift;
	my $select = shift;
	return 0 if ( $self->{'local'} );

	#print "squeue -u $self->{'username'}\n";
#	open( IN, "squeue -u $self->{'username'} |" ) or die $!;
#	my @IN = <IN>;
#	close(IN);
#	## pending
#	## I also need to take care about the different partititions:
#	#print "In total ".scalar(@IN)." jobs\n";
#	if ( $self->{'partitition'} ) {
#		@IN = grep( /\s$self->{'partitition'}\s/, @IN );
#
#		#print scalar(@IN)." jobs in partitcion $self->{'partitition'}\n";
#	}
#	if ($select) {
#		return scalar( grep ( /\s$select\s/, @IN ) );
#	}

	#print join("", @IN );
	## all
	return $self->{'slurmscripts'}->waiting();
#	return scalar(@IN) - 1;
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

	if ( $self->{'options'}->value('partitition') )
	{    # a hack if the script already uses p for something else
		$self->{'options'}
		  ->add( 'p', $self->{'options'}->value('partitition'), 1 );
		$self->{'options'}->drop('partitition');
	}
	elsif ( $self->{'options'}->value('A') =~ m/^lu/ )
	{    ## add one for each partitition you want to explicitly support
		$self->{'options'}->add( 'p', 'lu', 1 );
	}
	if ( $self->{'options'}->value('w') ) {    ## that is a special one
		$ret .= "#SBATCH -w, --" . $self->{'options'}->value('w') . "\n";
	}
	if ( $self->{'options'}->value('A') eq "lsens2017-3-2" ){
		$self->{'options'}->value('A', 'lsens2018-3-3');
		warn "Please update your project info to 'lsens2018-3-3'\n";
	}
	my $tmp;
	foreach my $option (
		@{ $self->{'options'}->{'required'} },
		@{ $self->{'options'}->{'optional'} }
	  )
	{
#		if ($option eq "A" ) {
#			Carp::confess ( "I got the option 'A': ".$self->{'options'}->value($option)."\n");
#		}
		next unless ( defined $self->{'options'}->value($option) );
		next if ( $option eq "w" );    #special one processed earlier
		$tmp = "-";
		$tmp = '--' if ( length($option) > 1 );
		$ret .=
		  "#SBATCH $tmp$option " . $self->{'options'}->value($option) . "\n";

	}
	#Tina2_ce
	my $oldN = $name;
	
	$name = substr( $name, 0,4 ). substr( md5_hex($name), 0,4 );
	warn "$oldN reports as $name\n";
	
	$ret .= join( "\n",
		"#SBATCH -J $name",
		"#SBATCH -o $name" . ".%j.out",
		"#SBATCH -e $name" . ".%j.err" );
	if ( @{ $self->{'SLURM_modules'} } ) {
		$ret .= "\n" . $self->load_SLURM_modules();
	}
	if ( ref($cmd) eq "ARRAY" ) {
		$ret .= "\n" . shift(@$cmd);
		chomp($ret);
		$ret .= "\ndate\n";
		for ( my $i = 0 ; $i < @$cmd ; $i++ ) {
			$ret .= @{ $self->{'sub_SLURMS'} }[$i]->subscript( @$cmd[$i] );
		}
	}
	else {
		chomp($ret);
		$ret .= "\n$cmd\ndate\n";
	}
	$ret .= "exit 0\n";
	return $ret;
}

=head2 subscript ( $cmd )

Create only the module loads and the cmd parts of the script.
This function is used internally to create scripts that need 
different sets of modules loaded at different timepoints.

If there is no command other than copying or moving files the modules will not be loaded.

=cut

sub subscript {
	my ( $self, $cmd ) = @_;
	
	unless ( $cmd =~ m/\w/ ) {
		return '';
	}
	my $ret = "## here comes a subscript entry:\n";
	my $OK = { 'cp' => 1 ,'mv' => 1, 'date' => 1, 'module' => 1, 'result' => 0 };
	my @tmp;
	foreach ( split( "\n", $cmd ) ){
		next unless ( $_ =~m/\w/);
		@tmp = split(/\s/, $_ );
		unless ( $OK ->{$tmp[0]} ) {
			$OK ->{'result'} = 1
		}
	}
	if ( @{ $self->{'SLURM_modules'} } and $OK->{'result'} ) {
		$ret .= $self->load_SLURM_modules();
	}
	chomp($ret);
	chomp($cmd);
	$ret .= "\n" . $cmd . "\ndate\n";
	return $ret;
}

sub load_R_x11 {
	my ($self) = @_;
	my $cmd = $self->load_SLURM_modules( 'ifort/2016.1.150-GCC-4.9.3-2.25',
		'impi/5.1.2.150', 'R/3.2.3-libX11-1.6.3' );
	print "please run\n$cmd" if ( $cmd =~ m/\w/ );
	return $cmd;
}

sub restore_SLURM_state{
	my ( $self, $module ) = @_;
	return "module purge\nmodule restore $module\n";
}

sub load_SLURM_modules {
	my ( $self, @modules ) = @_;
	#return '' if ( $self->{'local'} );
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
#		print "Hope there is some loaded module?: "
#		  . root->print_perl_var_def($loaded)
#		  if ( $self->{debug} );

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
	warn "I am removing slurm options from the options hash\n";
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
	my ( $self, $cmd, $fm, $ofile ) = @_;
	if ( defined $ofile ) {
		unless ( -f $ofile ) {
			return $self->run_notest( $cmd, $fm );
		}
		else {
			if ( ref($fm) eq "HASH" ) {
				$fm = $fm->{'total'};
			}
			warn "outfile was present ($ofile), script not run.\n";
		}
	}
	else {
		return $self->run_notest( $cmd, $fm );
	}
	return 1;
}

sub get_files_from_path {
	my ( $self, $path, @matches ) = @_;
	opendir( DIR, $path )
	  or die "I could not read from the directory '$path'\n$!\n";
	my @dat = grep { !/^\./ } readdir(DIR);
	closedir(DIR);

	#print "I get files from path '$path'\n". join("\n", @dat[0..2],'.','.', $dat[$#dat])."\n";
	my @ret;
	foreach my $select (@matches) {

		#print "I select all files matching '$select'\n";
		# . join( "\n", @dat ) . "\n";
		push(@ret, grep { /$select/ } @dat);

		#print "Still in the game:" . join( "\n", @ret ) . "\n\n";
	}
	#print "selected:\n" . join( "\n", @ret ) . "\n\n";
	@ret = map { "$path/$_" } @ret;
	return @ret;
}



sub run_notest {
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
	@OK = grep ( !/^\s*$/,   @OK );
	@OK = grep ( !/^module/, @OK )
	  ;    ## the module load lines should also not make the script run.
	if ( @OK > 0 ) {
		return $self->runScript( "$fm->{path}/$fm->{'filename_core'}.sh" ); 
	}
	else {
		print
"All outfiles present for $fm->{path}/$fm->{'filename_core'}.sh - not run\n";
		return 0;
	}
}

## this functionallity has to be replaced by a sqlite3 database.
## create the scripts and add them to a database.
## some other entity needs to run the scripts as we are likely in a docker or singularity image nowadays.

sub runScript{
	my ( $self, $scriptfile, $outfile ) = @_;
	$outfile ||='';
	
	unless ( -f $scriptfile ) {
		Carp::confess( "I need a script file to run on the cluster\nNot '$scriptfile'\n" );
	}
	if ( -f $outfile ){
		print "Outfile $outfile is present - script not run\n";
		return -1;
	}
	## And here I will check once more whether the script is actually containing any interesting function calls:
	open ( TEST , "<$scriptfile") or die $!;
	my $important = 0;
	while ( <TEST> ) {
		next if ($_ =~ m/^#/);
		chomp();
		next if ( $_ =~m/\s*exit 0\s*/ );
		next if ( $_ =~m/\s*date\s*/ );
		next if ( $_ =~m/^\s*module load / or $_ =~ m/^\s*ml / );
		next if ( $_ =~m/^\s*module purge/ );
		$important ++;
	}
	close ( TEST );
	unless ( $important ) {
		warn "Script did not contain any command apart from module calls data or exit statements\n\t$scriptfile not run\n";
		return -1;
	}
	my $fm = root->filemap($scriptfile);
	unless ( $self->{'debug'} ) {
		if ( $self->{'local'} or $self->{'run_local'}) {
			## add a local stdout and stderr file like slurm does.
			
			system( $self->{'shell'}
				  . " $scriptfile "
				  . "2> $fm->{path}/$fm->{'filename_core'}.local$$.err "
				  . " > $fm->{path}/$fm->{'filename_core'}.local$$.out " );
			return $$;
		}
		else {    ## use slurm pipeline
			
			if ( $self->in_pipeline() >= $self->{'max_jobs'} ) {
				print
"More than $self->{'max_jobs'} jobs in the pipeline - waiting:\n";
				local $| = 1;
				while ( $self->in_pipeline() >= $self->{'max_jobs'} ) {
					print ".";
					sleep( $self->{'sleep_sec'} );
				}
				print "\n";
			}
			print "sbatch $scriptfile\n";
			map { unlink($_) }
			  $self->get_files_from_path( './', "$fm->{'filename_core'}.\\d*.err",
				"$fm->{'filename_core'}.\\d*.out" );
			
			$self->{'slurmscripts'}->AddDataset({'cmd' => "sbatch $fm->{'total'}", 'Finished' => 0 } );
			#### get me here
			print "Added cmd to database\n";
			#open( PID, "sbatch $scriptfile |" );
			#my $tmp = join( "", <PID> );
			#close ( PID );
			#print $tmp;
			#if ( $tmp =~ m/Submitted batch job (\d+)/ ) {
			#	print "Submitted batch job '$1'\n";
			#	return $1;
			#}
		}
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
	Carp::confess("I can not check the outfile for '$cmd' - undefined files array!\n")
	  unless ( defined $outfiles[0] );
	foreach my $outfile (@outfiles) {
		if ( -f $outfile ) {
			warn "outfile '$outfile' is present - I will not re-create it!\n"
			  if ( $self->{'debug'} );
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
		unless ( defined $self->{$_} ) {
			$error .= "MISSING $type option $_\n";
		}
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
				print "get option $1 and value $2\n";
				$self->{'options'}->add( $1, $2 ) unless ( $ignore->{$1} );
			}
		}
	}
	warn( $self->{'options'}->AsString() );
	return $self;
}

1;
