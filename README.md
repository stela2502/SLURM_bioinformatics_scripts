# SLURM_bioinformatics_scripts
SLURM_bioinformatics_scripts is perl libraray that helps to run a bioinformatics core facility.

## INSTALL

As you most likely run the scripts on a server you are not admin of you might consider installing Perl <a href='http://search.cpan.org/~apeiron/local-lib-1.008004/lib/local/lib.pm'>local::lib</a>.

In addition you need git:

You will need my base libs which you can get from <a href='https://github.com/stela2502/Stefans_Lib_Esentials'>my Stefans_Lib_Esentials repository</a>. And my stefans_libs::BAMfile from <a href='https://github.com/stela2502/stefans_libs-BAMfile'>the respective guthub page</a>.

Clone this repository and use the 
perl Makefile.PL
make
make install
combination.

## USAGE

The package comes with a set of scripts. They report there usage on the command line and should be in your path after the install.

bowtie2_run_aurora.pl would return (at the time of writing the readme)

the cmd line switch -files is undefined!
the cmd line switch -options is undefined!
the cmd line switch -genome is undefined!
.
Usage:
        bowtie2_run_aurora.pl
           -files     :a list of fastq or fastq.gz files   
           -options   :additional options in the format
                       key_1 value_1 key_2 value_2 ... key_n value_n
           -genome    :the genome definition as used by bowtie2 -X
           -coverage  :the genome coverage file to create bigwig tracks
           -paired    :if you have paired data to map use this option and
                       give me the pired files one after the other
                   
           -bigwigTracks :If I have a coverage file I can create the 
                          bigwig tracks information and stoire it there

           -help      :print this help
           -debug     :verbose output
    
        required options:
    
        n   :number of cores to request from SLURM (rec. 10)
        N   :number of nodes to request from SLURM (rec. 1)
        mem :minimum memory from SLURM (rec. 4000M for e.g. human genome)
        t   :the time the job should get to finish - do not set too short
             rec 02:00:00 == 2h

