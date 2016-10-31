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

    Zargos2pwm_plot.pl
       -infile       :<please add some info!>
       -outpath       :<please add some info!>
       -options     :<please add some info!> you can specify more entries to that
                         format: key_1 value_1 key_2 value_2 ... key_n value_n


       -help           :print this help
       -debug          :verbose output
   
=head1 DESCRIPTION

  read a zargos output file and create sequence logos from that using the R seqLogo package.

  To get further help use 'Zargos2pwm_plot.pl -help' at the comman line.

=cut

use Getopt::Long;
use Pod::Usage;
use stefans_libs::flexible_data_structures::data_table;

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
	$error .= "the cmd line switch -options is undefined!\n";
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

$task_description .= 'perl '.$plugin_path .'/Zargos2pwm_plot.pl';
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
open ( LOG , ">$outpath/".$$."_Zargos2pwm_plot.pl.log") or die $!;
print LOG $task_description."\n";
close ( LOG );


## Do whatever you want!

open ( IN , "<$infile" ) or die "I could not open the infile '$infile'\n";
my ($logo,@line, $logo_id);
$logo_id = 1;
while ( <IN> ) {
	if ( $_ =~ m/^[P\d]\d\s/ ) {
		print "Interesting line $_\n";
		@line = split( /\s+/, $_ );
		if ( $line[0] eq "P0" ){
			if ( defined $logo ) {
				&create_logo( $logo );
			}
			$logo = data_table->new();
			$logo -> Add_2_Header( [ qw(A C G T) ] );
		}
		else {
			$logo->AddDataset ( { 'A' => $line[1], 'C' => $line[2], 'G' => $line[3], 'T' => $line[4]} );
		}
	}
}
close ( IN );
if ( defined $logo ) {
	&create_logo( $logo );
}


print "Done\n";

sub create_logo {
	my ( $data_table ) = @_;
	for  ( my $i = 0; $i < $data_table->Lines(); $i ++ ){
		&div ( @{$data_table->{'data'}}[$i] );
	}
	my $sum = &sum(@{$data_table->{'data'}}[0] ); 
	$data_table = $data_table->Transpose($data_table);
	open ( OUTR, ">$outpath/$logo_id"."_$sum.R" ) or die $!;
	print OUTR "library(seqLogo)\n"."proportion <- function(x){\n"."   rs <- sum(x);\n"."   return(x / rs);\n"."}\n";
	foreach my $l ( @{$data_table->{'data'}} ) {
		print OUTR "@$l[0] <- c(".join(", ", @$l[1..(@$l-1)] ).")\n"; 
	}
	print  OUTR "df <- data.frame( ".join( ",", @{$data_table->GetAsArray(0)})." )\n"."pwm <- apply(df, 1, proportion)\n"."pwm <- makePWM(pwm)\n"
	."pdf( file = '$outpath/$logo_id"."_$sum.pdf', width=6, height=6)\n" . "seqLogo(pwm)\n". "dev.off()\n";
	close ( OUTR );
	
	system( "R CMD BATCH --no-save --no-restore --no-readline -- $outpath/$logo_id"."_$sum.R" );
	
	$logo_id++;
}

sub sum {
	my @val = @{$_[0]};
	my $ret = 0; 
	map { $ret += $_ } @val;
	return $ret;
}
sub div {
	my @val = @{$_[0]};
	my $sum = &sum($_[0]);
	map { $val[$_] = $val[$_]/$sum } 0..(@val-1);
}