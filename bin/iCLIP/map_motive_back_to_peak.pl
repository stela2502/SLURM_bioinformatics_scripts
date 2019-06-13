#! /usr/bin/perl -w

=head1 LICENCE

  Copyright (C) 2016-11-01 Stefan Lang

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

    map_motive_back_to_peak.pl
       -peaks       :the peak bed file that was used to create the -motives file
       -motives     :the transfac motiv file
       -outfile     :the bedGraph outfile
       
       -separate    :match each motive separately back to the peaks

       -help           :print this help
       -debug          :verbose output
   
=head1 DESCRIPTION

  This script used the peak bed file together with the zagros output to create a binding site bedgraph file where the numer represents the percental overlap with the most likely consensus binging motive.

  To get further help use 'map_motive_back_to_peak.pl -help' at the comman line.

=cut

use Getopt::Long;
use Pod::Usage;

use strict;
use warnings;

use stefans_libs::file_readers::bed_file;
use stefans_libs::file_readers::transfaq;

use FindBin;
my $plugin_path = "$FindBin::Bin";

my $VERSION = 'v1.0';


my ( $help, $debug, $database, $peaks, $separate, $motives, $outfile);

Getopt::Long::GetOptions(
	 "-peaks=s"    => \$peaks,
	 "-motives=s"    => \$motives,
	 "-outfile=s"    => \$outfile,
	 "-separate"     => \$separate,

	 "-help"             => \$help,
	 "-debug"            => \$debug
);

my $warn = '';
my $error = '';

unless ( defined $peaks) {
	$error .= "the cmd line switch -peaks is undefined!\n";
}
unless ( defined $motives) {
	$error .= "the cmd line switch -motives is undefined!\n";
}
unless ( defined $outfile) {
	$error .= "the cmd line switch -outfile is undefined!\n";
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

$task_description .= 'perl '.$plugin_path .'/map_motive_back_to_peak.pl';
$task_description .= " -peaks '$peaks'" if (defined $peaks);
$task_description .= " -motives '$motives'" if (defined $motives);
$task_description .= " -outfile '$outfile'" if (defined $outfile);
$task_description .= " -separate" if ( defined $separate);


open ( LOG , ">$outfile.log") or die $!;
print LOG $task_description."\n";
close ( LOG );


## Do whatever you want!
my $peaks_data =  stefans_libs::file_readers::bed_file->new( {'filename' => $peaks });
my $motives_data = stefans_libs::file_readers::transfaq ->new();
$motives_data->read_file( $motives );

my $data;
if ( $separate ) {
	my $fm = root->filemap( $outfile );
	for ( my $id = 0; $id < @{$motives_data->{'motives'}}; $id ++ ){
		$data = $motives_data->map_motives_to_peaks( $peaks_data, $id);
		print "I write file ". $fm->{'path'}."/".$fm->{'filename_core'}."_motiv_$id.bedGraph"."\n";
		$data->write_file( $fm->{'path'}."/".$fm->{'filename_core'}."_motiv_$id.bedGraph")
	}
}
else {
	$data = $motives_data->map_motives_to_peaks( $peaks_data );
	$data->write_file( $outfile );
}


