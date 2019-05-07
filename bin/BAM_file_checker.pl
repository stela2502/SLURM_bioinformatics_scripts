#! /usr/bin/perl -w

=head1 LICENCE

  Copyright (C) 2017-01-30 Stefan Lang

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

    BAM_file_checker.pl
       -files     :a list of bam/sam files to be flagstat-ed
       -outfile   :the report outfile containing all failures
       -n         :how many processes to start for the samttols call (default 4)


       -help           :print this help
       -debug          :verbose output

=head1 DESCRIPTION

  This script runs samtools flagstat and removes all bam files that throw an error.

  To get further help use 'BAM_file_checker.pl -help' at the comman line.

=cut

use Getopt::Long;
use Pod::Usage;

use strict;
use warnings;

use FindBin;
my $plugin_path = "$FindBin::Bin";

use stefans_libs::flexible_data_structures::data_table;
use Parallel::ForkManager;

my $VERSION = 'v1.0';


my ( $help, $debug, $database, @files, $outfile, $n);

Getopt::Long::GetOptions(
       "-files=s{,}"    => \@files,
	 "-outfile=s"    => \$outfile,
	 "-n=s"    => \$n,

	 "-help"             => \$help,
	 "-debug"            => \$debug
);

my $warn = '';
my $error = '';

unless ( defined $files[0]) {
	$error .= "the cmd line switch -files is undefined!\n";
}
unless ( defined $outfile) {
	$error .= "the cmd line switch -outfile is undefined!\n";
}
unless ( defined $n) {
	$n = 4;
	$warn .= "the cmd line switch -n is undefined!\n";
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



my ( $task_description);

$task_description .= 'perl '.$plugin_path .'/BAM_file_checker.pl';
$task_description .= ' -files "'.join( '" "', @files ).'"' if ( defined $files[0]);
$task_description .= " -outfile '$outfile'" if (defined $outfile);
$task_description .= " -n '$n'" if (defined $n);



open ( LOG , ">$outfile.log") or die "I could not create the logfile '$outfile.log'\n$!\n";
print LOG $task_description."\n";
close ( LOG );


## Do whatever you want!

my $pm = Parallel::ForkManager->new($n);

open (OUT, ">$outfile" ) or die "I could not create the outfile '$outfile'\n$!\n";
print OUT "#filename\tsamtools flagstat error\n";

FILES:
foreach my $file ( @files ) {
		my $pid = $pm->start and next FILES;
		unless ( -f "$file.flagstat" ) {
			system ( "samtools flagstat $file 1> $file.flagstat 2> $file.flagstat.err" );
			open ( IN , "< $file.flagstat.err" ) or die "I could not open the outfile '$file.flagstat.err'\n$!\n";
			my $err = join("", <IN> );
			close ( IN );
			if ( $err  =~ m/command not found/ ) {
				map{ unlink($_) if ( -f $_) } "$file.flagstat", "$file.flagstat.err";
				Carp::confess( "This function can not succeed as the analysis command fails with:\n$err\n");
			}
			if ( $err =~ m/\w/ ) {
				$err =~ s/\n/ /g;
				print OUT "$file\t$err\n";
				unless ( $debug ) {
					map{ unlink($_) if ( -f $_) } "$file.flagstat", "$file.flagstat.err", $file;
				}
			}
		}
		$pm->finish; # Terminates the child process

}

$pm->wait_all_children;

close ( OUT );

print "report in file $outfile\n";
