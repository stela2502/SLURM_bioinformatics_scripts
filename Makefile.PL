#!/usr/bin/env perl
# IMPORTANT: if you delete this file your app will not work as
# expected.  You have been warned.

#use ExtUtils::MakeMaker;
use inc::Module::Install;

name  'SLURM_bioinformatic_scripts';
all_from 'lib/SLURM_bioinformatic_scripts.pm';
author 'Stefan Lang <Stefan.Lang@med.lu.se>';

version_from "lib/SLURM_bioinformatic_scripts.pm";

requires 'Stefans_Lib_Essentials' => 0;
requires 'File::Copy' => 0;
requires 'Parallel::ForkManager' => 0;

opendir( DIR, "bin/" ) or die "I could not open the bin folder\n$!\n";

map { install_script "bin/$_" } grep !/^\./,  grep '*.pl', readdir(DIR);

closedir(DIR);

auto_install();
WriteAll;


