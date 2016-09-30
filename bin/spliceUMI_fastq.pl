#! /usr/bin/perl -w

=head1 LICENCE

  Copyright (C) 2016-07-05 Stefan Lang

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

    spliceUMI_fastq.pl
       -UMI_read         :the fastq file containing the UMI read
       -data_read        :the fastq file containing the data
       -sampleDefinition :the sample definitions to decode the UMI reads
       
       -options     : format: key_1 value_1 key_2 value_2 ... key_n value_n
       
       			UNUSED at the moment
       			
       -outpath          :the outpath to store the spliced fastq files in


       -help           :print this help
       -debug          :verbose output
   
=head1 DESCRIPTION

  uses two fastq files to splice one big fastq into several sample specific ones and on the same time adds the UMI definitions to the fastq name.

  To get further help use 'spliceUMI_fastq.pl -help' at the comman line.

=cut

use Getopt::Long;
use Pod::Usage;
use stefans_libs::flexible_data_structures::data_table;

use strict;
use warnings;

use FindBin;
my $plugin_path = "$FindBin::Bin";

my $VERSION = 'v1.0';


my ( $help, $debug, $database, $UMI_read, $data_read, $sampleDefinition, $options, @options, $outpath);

Getopt::Long::GetOptions(
	 "-UMI_read=s"    => \$UMI_read,
	 "-data_read=s"    => \$data_read,
	 "-sampleDefinition=s"    => \$sampleDefinition,
       "-options=s{,}"    => \@options,
	 "-outpath=s"    => \$outpath,

	 "-help"             => \$help,
	 "-debug"            => \$debug
);

my $warn = '';
my $error = '';

unless ( defined $UMI_read) {
	$error .= "the cmd line switch -UMI_read is undefined!\n";
}
unless ( defined $data_read) {
	$error .= "the cmd line switch -data_read is undefined!\n";
}
unless ( defined $sampleDefinition) {
	$error .= "the cmd line switch -sampleDefinition is undefined!\n";
}
unless ( defined $options[0]) {
	$warn .= "the cmd line switch -options is undefined!\n";
}
unless ( defined $outpath) {
	$error .= "the cmd line switch -outpath is undefined!\n";
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

$task_description .= 'perl '.$plugin_path .'/spliceUMI_fastq.pl';
$task_description .= " -UMI_read '$UMI_read'" if (defined $UMI_read);
$task_description .= " -data_read '$data_read'" if (defined $data_read);
$task_description .= " -sampleDefinition '$sampleDefinition'" if (defined $sampleDefinition);
$task_description .= ' -options "'.join( '" "', @options ).'"' if ( defined $options[0]);
$task_description .= " -outpath '$outpath'" if (defined $outpath);


for ( my $i = 0 ; $i < @options ; $i += 2 ) {
	if ( defined $options[$i+1] ){
		$options[ $i + 1 ] =~ s/\n/ /g;
		$options->{ $options[$i] } = $options[ $i + 1 ];
	}
}
###### default options ########
$options->{'BC_length'} ||= 6;
$options->{'UMI_length'} ||= 14;
$options->{'BC_column'} ||= 'BC sequence';
##############################
mkdir( $outpath ) unless ( -d $outpath );
open ( LOG , ">$outpath/".$$."_spliceUMI_fastq.pl.log") or die $!;
print LOG $task_description."\n";
close ( LOG );


## Do whatever you want!

my $samples = data_table->new({'filename' => $sampleDefinition , 'no_doubble_cross'=> 1});
unless ( defined $samples->Header_Position( $options->{'BC_column'}) ) {
	die &helpString ( "The sample table does not contain the required BC_column '$options->{'BC_column'}' (name taken from -options BC_column '$options->{'BC_column'}'\n".join(", ",@{$samples->{'header'}}));
}

my ($info,@str, $valid, @entry, $files, $tmp, $gzipped );
$valid = { map { $_ => 1} @{$samples->GetAsArray($options->{'BC_column'})} };
my ( $UMI, $D );
if ( $UMI_read =~ /gz$/ ) {
	$gzipped = 1;
	open ( $UMI , "zcat $UMI_read |") or die "I could not open the UMI read file $UMI_read\n$!\n";
}else {
	open ( $UMI , "<$UMI_read") or die "I could not open the UMI read file $UMI_read\n$!\n";
}
if ( $data_read =~ /gz$/ ) {
	$gzipped = 1;
	open ( $D , "zcat $data_read |") or die "I could not open the data read file $data_read\n$!\n";
}else {	
	open ( $D , "<$data_read") or die "I could not open the data read file '$data_read'\n$!\n";
}


## code taken from http://stackoverflow.com/questions/2498937/how-can-i-walk-through-two-files-simultaneously-in-perl
sub read_file_line {
  my $fh = shift;
  if ($fh and my $line = <$fh>) {
  	return $line
  }
  return undef;
}


## fastq 
## 1 -> useless ID
## 2 -> sequence
## 3 -> orientation (?)
## 4 -> quality
my $i = 0;
my @BC = (0..($options->{'BC_length'}-1));
my @UMI = ($options->{'BC_length'}..($options->{'BC_length'}+($options->{'UMI_length'}-1)) );

$data_read =~ s/.fa?s?t?q.?g?z?$/.fastq/;
my $fm = root->filemap( $data_read );
my $outbase = join("/", $outpath, $fm->{'filename_core'} );

my ($pair1,@tmp);

while ($pair1 = read_file_line($UMI) and $tmp = read_file_line($D)) {
	$i ++;
	push ( @entry, $tmp );
	if ( $i == 2 ) { ## sequence
		chomp($pair1);
		@str = split("", $pair1 );
		#print "Where is the uninitialized entry? \nA= '".join("', '", @str[@BC])."'\nB= '".join("', '", @str[@UMI])."\n";
		$info = [join("", @str[@BC]), join("", @str[@UMI]) ];
	}elsif ( $i == 4){
		$i = 0;
		unless ( $valid ->{@$info[0]} ){
			@$info[0] = "NoSample";
		}
		unless ( defined $files->{@$info[0]} ){
			if ( $gzipped ) {
				open ( $files->{@$info[0]} ,"| /bin/gzip -c > $outbase.@$info[0].fastq.gz" ) or die "I could not create the outfile '$outbase.@$info[0].fastq.gz'\n$!\n";
			}else {
				open ( $files->{@$info[0]} , ">$outbase.@$info[0].fastq" ) or die "I could not create the outfile '$outbase.@$info[0].fastq'\n$!\n";
			}
			
		}
		@tmp = split(" ", $entry[0]);
		$tmp[0] .= ":@$info[1]";
		$entry[0] = join(" ", @tmp )."\n";
		print ${$files->{@$info[0]}} join("", @entry );
		@entry = ();
		$i = 0;
	}
}

close ( $UMI );
close ( $D );
foreach ( values %$files ) { 
	close ( $_ );
}

print "Done!\n";


