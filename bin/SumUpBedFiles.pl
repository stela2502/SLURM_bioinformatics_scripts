#! /usr/bin/perl -w

=head1 LICENCE

  Copyright (C) 2017-02-08 Stefan Lang

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

=head1  SYNOPSIS

    SumUpBedFiles.pl
       -bed_files   :the list of bed files to merge
       -outfile     :the final merged bed file
       -options     :format key space value
     
              N: Amount of nodes for the sum up process (default 1)
              n: Amount of cores for the sub up process (default 1) 
                 increase this if you encounter memory problems
              t: time before one SLURM process is stopped ('00:04:00')
              A: the SLURM account to use ('lu2016-2-7')

       -help           :print this help
       -debug          :verbose output
   
=head1 DESCRIPTION

  Sum up a huge numer of bed files using a SLURM cluster (here aurora)

  To get further help use 'SumUpBedFiles.pl -help' at the comman line.

=cut

use Getopt::Long;
use Pod::Usage;

use strict;
use warnings;

use stefans_libs::SLURM;

use FindBin;
my $plugin_path = "$FindBin::Bin";

my $VERSION = 'v1.0';


my ( $help, $debug, $database, @bed_files, $outfile, $options, @options);

Getopt::Long::GetOptions(
       "-bed_files=s{,}"    => \@bed_files,
	 "-outfile=s"    => \$outfile,
       "-options=s{,}"    => \@options,

	 "-help"             => \$help,
	 "-debug"            => \$debug
);

my $warn = '';
my $error = '';

unless ( -f $bed_files[0]) {
	$error .= "the cmd line switch -bed_files is undefined!\n";
}
unless ( defined $outfile) {
	$error .= "the cmd line switch -outfile is undefined!\n";
}
unless ( defined $options[0]) {
	$warn .= "the cmd line switch -options is undefined!\n";
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

$task_description .= 'perl '.$plugin_path .'/SumUpBedFiles.pl';
$task_description .= ' -bed_files "'.join( '" "', @bed_files ).'"' if ( defined $bed_files[0]);
$task_description .= " -outfile '$outfile'" if (defined $outfile);
$task_description .= ' -options "'.join( '" "', @options ).'"' if ( defined $options[0]);


for ( my $i = 0 ; $i < @options ; $i += 2 ) {
	$options[ $i + 1 ] =~ s/\n/ /g;
	$options->{ $options[$i] } = $options[ $i + 1 ];
}
###### default options ########
$options->{'N'} ||= 1;
$options->{'n'} ||= 1;
$options->{'t'} ||= '00:04:00';
$options->{'A'} ||= "lu2016-2-7";
###############################

### predefined options #######
$options->{'SLURM_modules'} = [ 'GCC/4.9.3-2.25', 'OpenMPI/1.10.2', 'BEDTools/2.25.0'];
$options->{'debug'} = $debug;
##############################

unless ( -f  $outfile.log ) { # that should only be created by the first process
	open ( LOG , ">$outfile.log") or die $!;
	print LOG $task_description."\n";
	close ( LOG );
}

## Do whatever you want!

my $SLURM = stefans_libs::SLURM->new($options);

#cat fileA  fileB fileC |
#sort -k1,1 -k2,2n| bedtools merge -c 5 -o count  > resultF
my $pid;
if ( @bed_files < 4 ) {
	&merge( $outfile, @bed_files );
	my $fm = root->filemap($outfile );
	print "Outfile '$outfile' should have been created\n";
	if ( $bed_files[0] =~m/merge\d+_\d+\.bed/ ){
		print "you now can savely delete all intermediate files using:\n".
			"rm -f $fm->{'path'}/merge[0-9][0-9]*_[0-9][0-9]*.*\n"
	}
}
else {
	## now I need to merge 2 to 3 files together into one new tmp file
	my @pids;
	my $fm = root->filemap($outfile );
	my $tmp_file;
	my @tmp_files;
	while( @bed_files > 3 ) {
		$tmp_file = $fm->{'path'}."/merge".$$."_".scalar(@pids).".bed";
		push( @pids, &merge( $tmp_file, shift(@bed_files), shift(@bed_files)  ) );
		push ( @tmp_files, $tmp_file );
	}
	$tmp_file = $fm->{'path'}."/merge".$$."_".scalar(@pids).".bed";
	push( @pids, &merge( $tmp_file, @bed_files  ) );
	push ( @tmp_files, $tmp_file );
	## wait for all temp files to be created and
	## restart the sum process on the merged tmp files untill I can create the final outfile using the first option (if( @bed_files < 4 ) ).
	$task_description = 'SumUpBedFiles.pl'
		. ' -bed_files '.join( ' ', @tmp_files )
		. " -outfile $outfile";
	$task_description .= ' -options '.join( ' ', @options ) if ( defined $options[0]);
	$task_description .= ' -debug' if ( $debug );
	my $cmd = "wait4pidsAndRun.pl -pids ".join(" ", @pids). " -cmd '$task_description' -options sleep 10";
	print "I now start the wait and run script\n$cmd\n";
	system( $cmd );
}

sub merge {
	my ( $outfile, @bed_files ) = @_;
	my $cmd = "cat ".join(" ",@bed_files)." | sort -k1,1 -k2,2n | bedtools merge > $outfile";
	print  "$cmd\n" if ( $debug );
	return $SLURM -> run ( $cmd, $outfile );
}
