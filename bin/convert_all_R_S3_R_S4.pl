#! /usr/bin/perl -w

=head1 LICENCE

  Copyright (C) 2018-09-17 Stefan Lang

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

    convert_all_R_S3_R_S4.pl
       -infiles     :<please add some info!> you can specify more entries to that
       -outpath       :<please add some info!>
       -className       :<please add some info!>


       -help           :print this help
       -debug          :verbose output
   
=head1 DESCRIPTION

  apply the script convert_R_S3_R_S4.pl on a list of R files

  To get further help use 'convert_all_R_S3_R_S4.pl -help' at the comman line.

=cut

use Getopt::Long;
use Pod::Usage;

use strict;
use warnings;

use stefans_libs::root;
use File::Spec;

use FindBin;
my $plugin_path = "$FindBin::Bin";

my $VERSION = 'v1.0';


my ( $help, $debug, $database, @infiles, $outpath, $className);

Getopt::Long::GetOptions(
       "-infiles=s{,}"    => \@infiles,
	 "-outpath=s"    => \$outpath,
	 "-className=s"    => \$className,

	 "-help"             => \$help,
	 "-debug"            => \$debug
);

my $warn = '';
my $error = '';

unless ( defined $infiles[0]) {
	$error .= "the cmd line switch -infiles is undefined!\n";
}
unless ( defined $outpath) {
	$error .= "the cmd line switch -outpath is undefined!\n";
}
unless ( defined $className) {
	$error .= "the cmd line switch -className is undefined!\n";
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

$task_description .= 'perl '.$plugin_path .'/convert_all_R_S3_R_S4.pl';
$task_description .= ' -infiles "'.join( '" "', @infiles ).'"' if ( defined $infiles[0]);
$task_description .= " -outpath '$outpath'" if (defined $outpath);
$task_description .= " -className '$className'" if (defined $className);



mkdir( $outpath ) unless ( -d $outpath );
open ( LOG , ">$outpath/".$$."_convert_all_R_S3_R_S4.pl.log") or die $!;
print LOG $task_description."\n";
close ( LOG );


## Do whatever you want!

my ( $outfile, $fm, $cmd) ;
foreach my $infile ( @infiles) {
	$fm = root->filemap($infile);
	$outfile = File::Spec->catfile($outpath, $fm->{'filename'} );
	$cmd = "convert_R_S3_R_S4.pl -R_source '$infile' -outfile '$outfile' -className '$className'";
	if ( $debug ){
		print "$cmd\n";
	}else {
		system( $cmd) ;
	}
}

print "Done\n";
