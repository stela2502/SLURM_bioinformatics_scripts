#! /usr/bin/perl -w

=head1 LICENCE

  Copyright (C) 2016-10-06 Stefan Lang

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

    download_fasta_from_Basespace_usingR.pl
       -access_token       :<please add some info!>
       -project_id       :<please add some info!>
       -opath       :<please add some info!>


       -help           :print this help
       -debug          :verbose output
   
=head1 DESCRIPTION

  Uses R to download basespace data that has been shared with you.

  To get further help use 'download_fasta_from_Basespace_usingR.pl -help' at the comman line.

=cut

use Getopt::Long;
use Pod::Usage;
use Cwd;

use strict;
use warnings;

use FindBin;
my $plugin_path = "$FindBin::Bin";

my $VERSION = 'v1.0';


my ( $help, $debug, $database, $access_token, $project_id, $opath);

Getopt::Long::GetOptions(
	 "-access_token=s"    => \$access_token,
	 "-project_id=s"    => \$project_id,
	 "-opath=s"    => \$opath,

	 "-help"             => \$help,
	 "-debug"            => \$debug
);

my $warn = '';
my $error = '';

unless ( defined $access_token) {
	warn "token for Stefan Lang set automaticly!\n";
	$access_token = '19b1dd5f50b344e3ac9f9c70e6ebb4fa';
}
unless ( defined $project_id) {
	$error .= "the cmd line switch -project_id is undefined!\n";
}
unless ( defined $opath) {
	$opath = getcwd;
	warn "opath set to $opath\n";
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

$task_description .= 'perl '.$plugin_path .'/download_fasta_from_Basespace_usingR.pl';
$task_description .= " -access_token '$access_token'" if (defined $access_token);
$task_description .= " -project_id '$project_id'" if (defined $project_id);
$task_description .= " -opath '$opath'" if (defined $opath);



open ( OUT , ">$opath/download.R" ) or die "I could not create the R script '$opath/download.R'\n$!\n";
print OUT "library(BaseSpaceR)
ACCESS_TOKEN<- '$access_token'
PROJECT_ID<- '$project_id'  ## Get proj ID from url of the project

aAuth<- AppAuth(access_token = ACCESS_TOKEN)
selProj <- Projects(aAuth, id = PROJECT_ID, simplify = TRUE)
sampl <- listSamples(selProj, limit= 1000)
inSample <- Samples(aAuth, id = Id(sampl), simplify = TRUE)
for(s in inSample){
    f <- listFiles(s, Extensions = '.gz')
    print(Name(f))
    getFiles(aAuth, id= Id(f), destDir = '$opath', verbose = TRUE)
} \n";

close( OUT );

system( "R CMD BATCH $opath/download.R &" );
print "The R script has been created and started - finshing here\n";

