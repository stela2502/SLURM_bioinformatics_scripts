#!/usr/bin/env perl
# IMPORTANT: if you delete this file your app will not work as
# expected.  You have been warned.

use ExtUtils::MakeMaker;



my $VERSION = 0.6;

#name  'SLURM_bioinformatic_scripts';
#all_from 'lib/SLURM_bioinformatic_scripts.pm';
#author 'Stefan Lang <Stefan.Lang@med.lu.se>';
opendir( DIR, "bin/" ) or die "I could not open the bin folder\n$!\n";
my @files = (  grep /\.pl$/, readdir(DIR) );
closedir(DIR);

#print  "These are the sripts I want to install: ". join("\n", @files )."\n";


WriteMakefile(
        NAME         => "SLURM_bioinformatic_scripts",
        VERSION_FROM => "lib/SLURM_bioinformatic_scripts.pm",
        EXE_FILES    => [ map {"bin/$_" } @files, 'iCLIP/map_motive_back_to_peak.pl'  ],
        PREREQ_PM    => {
                'Stefans_Libs_Essentials'               => 0,
		'File::Copy' => 0, 
		'Parallel::ForkManager' => 0,
		'stefans_libs::BAMfile' => 0,
          }

);

