#! /usr/bin/perl -w

=head1 LICENCE

  Copyright (C) 2019-03-27 Stefan Lang

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
   
   binCreate.pl from git@github.com:stela2502/Stefans_Lib_Esentials.git commit c35cfea822cac3435c5821897ec3976372a89673
   

=head1  SYNOPSIS

    backupFastq.pl
       -path       :the path to be backed up (default ./)
       -data_path  :the backup folder (default ~/DATA)
       -options    :options to be used later on (empty)
                    format: key_1 value_1 key_2 value_2 ... key_n value_n
            recursive 1: also recurse into downstream folders the downstream folders

       -help           :print this help
       -debug          :verbose output
   
=head1 DESCRIPTION

  copy all fastq.gz files to the DATA folder

  To get further help use 'backupFastq.pl -help' at the comman line.

=cut

use Getopt::Long;
use Pod::Usage;
use File::Spec;
use File::Path;

use strict;
use warnings;

use FindBin;
my $plugin_path = "$FindBin::Bin";

my $VERSION = 'v1.0';


my ( $help, $debug, $database, $path, $data_path, $options, @options, $outpath);

Getopt::Long::GetOptions(
	 "-path=s"    => \$path,
	 "-data_path=s"    => \$data_path,
       "-options=s{,}"    => \@options,
	 "-outpath=s"    => \$outpath,

	 "-help"             => \$help,
	 "-debug"            => \$debug
);

my $warn = '';
my $error = '';

unless ( defined $path) {
	open(PIPE, "pwd | " ) or die $!;
	$path = join("" , map{ chomp; $_} <PIPE>);
	close ( PIPE );
	$warn .= "The path was set to '$path'\n";
}
unless ( defined $data_path) {
	
	$data_path = File::Spec->catfile($path, 'NAS','');
#	open(PIPE, "ls -d ~/DATA | " ) or die $!;
#	$data_path = join("" , map{ chomp; $_} <PIPE>);
#	close ( PIPE );
	$warn .= "the data_path was set to '$data_path'\n";
}
unless ( defined $options[0]) {
	#$error .= "the cmd line switch -options is undefined!\n";
}
unless ( defined $outpath) {
	$outpath = File::Spec->catfile($path, 'DATA','');
#	open(PIPE, "whoami | " ) or die $!;
#	my $username = join("" , map{ chomp; $_} <PIPE>);
#	close ( PIPE );
#	my @tmp = split( "/", $path);
#	for (my $i = 0; $i < @tmp; $i ++){
#		if ( $tmp[$i] eq $username) {
#			splice(@tmp,0, $i+1);
#		}
#	}
#	$outpath = join("/", $data_path, @tmp);
	$warn .= "The outpath was set to '$outpath'\n";
#	unless ( -d $outpath ) {
#		system( "mkdir $outpath -p" );
#	}
}


if ( $help ){
	print helpString( ) ;
	exit;
}

if ( $error =~ m/\w/ ){
	helpString($error ) ;
	exit;
}

if ( $warn =~ m/\w/ ){
	warn $warn."\n";
	#exit;
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

$task_description .= 'perl '.$plugin_path .'/backupFastq.pl';
$task_description .= " -path '$path'" if (defined $path);
$task_description .= " -data_path '$data_path'" if (defined $data_path);
$task_description .= ' -options "'.join( '" "', @options ).'"' if ( defined $options[0]);
$task_description .= " -outpath '$outpath'" if (defined $outpath);


for ( my $i = 0 ; $i < @options ; $i += 2 ) {
	$options[ $i + 1 ] =~ s/\n/ /g;
	$options->{ $options[$i] } = $options[ $i + 1 ];
}
###### default options ########
$options->{'recursive'} ||= 0;
##############################
mkdir( $outpath ) unless ( -d $outpath );
open ( LOG , ">$outpath/".$$."_backupFastq.pl.log") or die $!;
print LOG $task_description."\n";
close ( LOG );


## Do whatever you want!
sub copyFastqsInPath {
	my ( $path, $outpath) = @_;
	my $i = 0;
	my ($DIR, $thisName, $targetName);
	opendir ($DIR, $path) or die $!;
	foreach my $fname (readdir( $DIR) ) {
		next if ( $fname =~ m/^\./);
		$thisName = File::Spec->catfile($path, $fname);
		$targetName = File::Spec->catfile($outpath, $fname);
		#print ( "processing $thisName\n");
		
		if ( -d  $thisName ){
			if ( $options->{'recursive'} ){
				$i += copyFastqsInPath($thisName, $targetName );
				#print ( "got a folder $thisName - recursion\n");
			}
			next;
		}
		next unless ( $fname =~ m/fastq\.gz$/);
		unless ( -f $targetName ){
			File::Path->make_path($outpath) unless ( -d $outpath);
			print   "cp \n\t$thisName \n\t$targetName\n";
			system( "cp $thisName $targetName" ) unless ( $debug );
		}
		if ( -f "$thisName.md5") {
			unless ( -f "$targetName.md5" ) {
				print "cp $fname.md5\n";
				system(  "cp $thisName.md5 $targetName.md5" ) unless ($debug);
			}
		}else {
			print ( "create the mds5sum for file $fname\n");
			system(  "md5sum $targetName > $targetName.md5") unless( $debug );
		}
		$i ++;
	}
	closedir($DIR);
	return $i;
}

my $i = 0;

$i = copyFastqsInPath( $data_path, $outpath);

print "$i files copied to the backup folder.\n";
