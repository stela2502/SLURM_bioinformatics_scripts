#! /usr/bin/perl -w

=head1 LICENCE

  Copyright (C) 2018-02-06 Stefan Lang

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

=head1 CREATED BY
   
   binCreate.pl from  commit 
   

=head1  SYNOPSIS

    run_cellRanger_10x.pl
       -modules   :a list of slurm modules to load
       -outpath   :the outpath for all files created
       -tmp       :the fast tmp path default to '$SNIC_TMP'
       -organism  :which organism the data is from (mouse or human)
       -gtf       :which annotation version to use (will report options if left empty)
       -indexFile :the index definition file (see https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/mkfastq)
       -options   :format: key_1 value_1 key_2 value_2 ... key_n value_n
       		
       		what :which cell ranger module to call (mkfastq, count or aggr; all if empty)
       		
       		n    :number of cores per node (default = 20 )
            t    :time untill the script is stopped (default =10:00:0 (10h))
            p    :partitition (default 'dell')
                   


       -help           :print this help
       -debug          :verbose output
   
=head1 DESCRIPTION

  Runs the cell ranger pipeline specificly on Lunarc aurora ls2

  To get further help use 'run_cellRanger_10x.pl -help' at the comman line.

=cut

use Getopt::Long;
use Pod::Usage;

use strict;
use warnings;

use stefans_libs::root;
use stefans_libs::SLURM;

use FindBin;
my $plugin_path = "$FindBin::Bin";

my $VERSION = 'v1.0';


my ( $help, $debug, $database, @modules, $outpath,$organism, $tmp, $indexFile, $options, @options, $gtf, $runPath);

Getopt::Long::GetOptions(
       "-modules=s{,}"    => \@modules,
	 "-outpath=s"    => \$outpath,
	 "-runPath=s" => $runPath,
	 "-tmp=s"    => \$tmp,
	 "-indexFile=s"    => \$indexFile,
	 "-organism=s"  => \$organism,
       "-options=s{,}"    => \@options,
     "-gtf=s" => \$gtf,

	 "-help"             => \$help,
	 "-debug"            => \$debug
);

my $warn = '';
my $error = '';

unless ( -d $runPath) {
	$error .= "the cmd line switch -runPath is undefined!\n";
}
unless ( defined $modules[0]) {
	$error .= "the cmd line switch -modules is undefined!\n";
}
unless ( defined $outpath) {
	$error .= "the cmd line switch -outpath is undefined!\n";
}
unless ( defined $tmp) {
	$tmp = '$SNIC_TMP';
	#$error .= "the cmd line switch -tmp is undefined!\n";
}
unless ( defined $indexFile) {
	$error .= "the cmd line switch -indexFile is undefined!\n";
}




if ( $help ){
	print helpString( ) ;
	exit;
}

if ( $error =~ m/\w/ ){
	helpString($error ) ;
	exit;
}

sub helpString {
	my $errorMessage = shift;
	$errorMessage = ' ' unless ( defined $errorMessage); 
	print "$errorMessage.\n";
	pod2usage(q(-verbose) => 1);
}

### initialize default options:

#$options->{'n'} ||= 10;

###


my ( $task_description);

$task_description .= 'perl '.$plugin_path .'/run_cellRanger_10x.pl';
$task_description .= ' -modules "'.join( '" "', @modules ).'"' if ( defined $modules[0]);
$task_description .= ' -runPath "$runPath"';
$task_description .= " -outpath '$outpath'" if (defined $outpath);
$task_description .= " -tmp '$tmp'" if (defined $tmp);
$task_description .= " -indexFile '$indexFile'" if (defined $indexFile);
$task_description .= " -gtf '$gtf'" if (defined $gtf);
$task_description .= " -organism '$organism'" if (defined $organism);
$task_description .= ' -options "'.join( '" "', @options ).'"' if ( defined $options[0]);


for ( my $i = 0 ; $i < @options ; $i += 2 ) {
	$options[ $i + 1 ] =~ s/\n/ /g;
	$options->{ $options[$i] } = $options[ $i + 1 ];
}
###### default options ########
$options->{'p'} ||= 'dell';
$options->{'N'}   = 1;
$options->{'n'} ||= 20;
$options->{'A'} ||= 'Plese set a useful default value';
$options->{'t'} ||= '10:00:00';

$options->{'what'} ||= '';
###############################


mkdir( $outpath ) unless ( -d $outpath );
open ( LOG , ">$outpath/".$$."_run_cellRanger_10x.pl.log") or die $!;
print LOG $task_description."\n";
close ( LOG );


## Do whatever you want!

my $run_out;

if ( $options->{'what'} eq "mkfastq" or $options->{'what'} eq "" ) {

	
}

## set the $run_out to the right position

if ( $options->{'what'} eq 'count' or $options->{'what'} eq "" ) {
	unless ( defined $organism ) {
		die "You can not run 'count' without giving me the organism to work with!\n";
	}
	my $path = "";
	opendir ( DIR, "") or die $!;
	my $opts = { map { $_ => 1} readdir(DIR) };
	closedir(DIR);
	unless ( defined $gtf){
		die "use one of these\n -options gtf ". join( ' or ', sort keys %$opts )."\n";
	}
	unless ( -d File::Spec->catfile( $path, $gtf ) ) {
		die "the path ".File::Spec->catfile( $path, $gtf )." does not exist, please use one gtf file from this list: ".join( ' or ', sort keys %$opts )."\n";
	}
	my $INDEX =  read_indexFile ( $indexFile );
	if ( $debug ) {
		print "The index file:\n".$INDEX->AsString();
	}
}
if ( $options->{'what'} eq 'aggr' or $options->{'what'} eq "" ) {
	
}


sub read_indexFile {
	#here I only need to get the sample specific parts...
	
	open( IN,  "<$indexFile") or die "I could not open the index file: '$indexFile'\n$!\n";
	my $data_table = stefans_libs::flexible_data_structures::data_table->new();
	my $use = 0;
	while ( <IN> ) {
		if ( $_ =~ m/^\[Data\]/ ) {
			$use = 1; 
			next;
		}
		if ( $use ) {
			chomp();
			if ( $data_table->{'last_warning'} eq "" ) {
				$data_table->{'last_warning'} = "All OK";
				$data_table -> Add_2_Header( [split(",", $_)]);
				next;
			}
			$_ =~ s/\s+/_/g;
			$_ =~ s/\-+/_/g;
			push (@{$data_table->{'data'}}, [split(",", $_)] );
		}
	}
	close ( IN );
	return $data_table;
}
