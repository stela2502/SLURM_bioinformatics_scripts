#! /usr/bin/perl -w

=head1 LICENCE

  Copyright (C) 2016-10-31 Stefan Lang

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

    transfaqMotive2SequenceLogo.pl
       -infile       :one transfaq formated zargos outfile
       -outpath      :the path where all figures should be plotted to (pdf)
       -options      :unused
                      format: key_1 value_1 key_2 value_2 ... key_n value_n


       -help           :print this help
       -debug          :verbose output
   
=head1 DESCRIPTION

  read a zargos output file and create sequence logos from that using the R seqLogo package.

  To get further help use 'transfaqMotive2SequenceLogo.pl -help' at the comman line.

=cut

use Getopt::Long;
use Pod::Usage;
use stefans_libs::file_readers::transfaq;

use strict;
use warnings;

use FindBin;
my $plugin_path = "$FindBin::Bin";

my $VERSION = 'v1.0';


my ( $help, $debug, $database, $infile, $outpath, $options, @options);

Getopt::Long::GetOptions(
	 "-infile=s"    => \$infile,
	 "-outpath=s"    => \$outpath,
       "-options=s{,}"    => \@options,

	 "-help"             => \$help,
	 "-debug"            => \$debug
);

my $warn = '';
my $error = '';

unless ( defined $infile) {
	$error .= "the cmd line switch -infile is undefined!\n";
}
unless ( defined $outpath) {
	$error .= "the cmd line switch -outpath is undefined!\n";
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

$task_description .= 'perl '.$plugin_path .'/transfaqMotive2SequenceLogo.pl';
$task_description .= " -infile '$infile'" if (defined $infile);
$task_description .= " -outpath '$outpath'" if (defined $outpath);
$task_description .= ' -options "'.join( '" "', @options ).'"' if ( defined $options[0]);


for ( my $i = 0 ; $i < @options ; $i += 2 ) {
	$options[ $i + 1 ] =~ s/\n/ /g;
	$options->{ $options[$i] } = $options[ $i + 1 ];
}
###### default options ########
#$options->{'something'} ||= 'default value';
##############################
mkdir( $outpath ) unless ( -d $outpath );
open ( LOG , ">$outpath/".$$."_transfaqMotive2SequenceLogo.pl.log") or die $!;
print LOG $task_description."\n";
close ( LOG );


## Do whatever you want!

my $OBJ = stefans_libs::file_readers::transfaq -> new({'debug' => $debug, 'filename' => $infile});
$OBJ -> create_logo( $outpath );

print "Done\n";

