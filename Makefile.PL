#!/usr/bin/env perl
# IMPORTANT: if you delete this file your app will not work as
# expected.  You have been warned.

#use ExtUtils::MakeMaker;
use inc::Module::Install;

my $VERSION = 0.6;

name  'SLURM_bioinformatic_scripts';
all_from 'lib/SLURM_bioinformatic_scripts.pm';
author 'Stefan Lang <Stefan.Lang@med.lu.se>';

#version_from "lib/SLURM_bioinformatic_scripts.pm";

#requires 'Stefans_Lib_Essentials' => 0;
#requires 'stefans_libs-GenomeDB' => 0;
requires 'File::Copy' => 0;
requires 'Parallel::ForkManager' => 0;
requires 'stefans_libs::BAMfile' => 0;

sub install_scripts_in_path {
	my ( $dir ) = @_;
	my $PATH;
	opendir( $PATH, "$dir/" ) or die "I could not open the bin folder\n$!\n";
	my @files =  readdir($PATH);
	closedir ( $PATH );
	warn "processing bin folder $dir\n";
	map { if ( $_ =~ m/^\./ ){
		} elsif ( -f "$dir/$_" and $_ =~ m/.pl$/ ){ 
			&install_script("$dir/$_"); 
		} elsif ( -d "$dir/$_" ) { 
			&install_scripts_in_path("$dir/$_") 
		} } @files;
	
}

&install_scripts_in_path('bin');

auto_install();
WriteAll();


