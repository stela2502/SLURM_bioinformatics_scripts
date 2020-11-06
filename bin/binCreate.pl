#! /usr/bin/perl -w

#  Copyright (C) 2008 Stefan Lang

#  This program is free software; you can redistribute it
#  and/or modify it under the terms of the GNU General Public License
#  as published by the Free Software Foundation;
#  either version 3 of the License, or (at your option) any later version.

#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
#  See the GNU General Public License for more details.

#  You should have received a copy of the GNU General Public License
#  along with this program; if not, see <http://www.gnu.org/licenses/>.

=head1 binCreate.pl

This script is used to craete new scripts with a first pod description and a helpString function.

To get further help use 'binCreate.pl -help' at the comman line.

=cut

use Getopt::Long;
use stefans_libs::root;
use strict;
use warnings;

use stefans_libs::Version;

my $VERSION = "v1.1";

my ( $help, $debug, $name, $pod, $force, @commandLineSwitches,
	$additional_checks );

Getopt::Long::GetOptions(
	"-pod=s"                    => \$pod,
	"-name=s"                   => \$name,
	"-force"                    => \$force,
	"-commandLineSwitches=s{,}" => \@commandLineSwitches,
	"-help"                     => \$help,
	"-debug"                    => \$debug
);

if ($help) {
	print helpString();
	exit;
}
unless ( defined $name ) {
	print helpString("ERROR: no name defined");
	exit;
}
unless ( defined $pod ) {
	print helpString("ERROR: no description for the executable");
	exit;
}

my ( $exec_name, $path, @file );

print "We got the name $name and the pod $pod\n" if ($debug);

@file      = split( "/", $name );
$exec_name = pop(@file);
$path      = join( "/", @file );
unless ( -d $path ) {
	die "the path '$path' does not exist. Please create it!";
}
$exec_name = "$exec_name.pl" unless ( $exec_name =~ m/\.pl$/ );

my ( $add_2_variable_def, $add_2_variable_read, $add_2_help_string,
	$task_string, $options_string, $optionNames, );

$add_2_variable_def = $add_2_variable_read = $add_2_help_string =
  $options_string = '';
$task_string = "\$task_description .= 'perl '.\$plugin_path .'/$exec_name';\n";
my $error_check = '';
my $log_str     = '';
$additional_checks = '';
foreach my $variableStr (@commandLineSwitches) {
	if ( $variableStr eq "options#array" ) {
		$add_2_variable_def .= ", \$options";
		$options_string = join( "\n",
			"for ( my \$i = 0 ; \$i < \@options ; \$i += 2 ) {",
			"\t\$options[ \$i + 1 ] =~ s/\\n/ /g;",
			"\t\$options->{ \$options[\$i] } = \$options[ \$i + 1 ];",
			"}",
			"###### default options ########",
			'#$options->{\'something\'} ||= \'default value\';',
			"##############################" );
		$additional_checks .= "### initialize default options:\n\n"
		  . "#\$options->{'n'} ||= 10;\n\n" . "###\n";
	}	
	if ( lc($variableStr) eq "outfile" ) {
		$log_str =
		"use stefans_libs::Version;\n"
		. "my \$V = stefans_libs::Version->new();\n"
		. "my \$fm = root->filemap( \$outfile );\n"
		. "mkdir( \$fm->{'path'}) unless ( -d \$fm->{'path'} );\n\n"
		. "open ( LOG , \">\$outfile.log\") or die \$!;\n"
		. "print LOG '#library version '.\$V->version( AddRightLibNameHere ).\"\\n\";\n"
		. "print LOG \$task_description.\"\\n\";\n"
		. "close ( LOG );\n\n";
	}
	elsif ( lc($variableStr) eq "outpath" ) {
		$log_str =
		    "mkdir( \$outpath ) unless ( -d \$outpath );\n"
		  . "open ( LOG , \">\$outpath/\".\$\$.\"_$exec_name.log\") or die \$!;\n"
		  . "print LOG \$task_description.\"\\n\";\n"
		  . "close ( LOG );\n\n";
	}
	if ( $variableStr =~ s/#array$// ) {
		$add_2_variable_def .= ", \@$variableStr";
		$add_2_variable_read .=
		  "       \"-$variableStr=s{,}\"    => \\\@$variableStr,\n";
		$add_2_help_string .=
"       -$variableStr     :<please add some info!> you can specify more entries to that\n";
		if ( $variableStr eq "options" ) {
			$add_2_help_string .=
"                         format: key_1 value_1 key_2 value_2 ... key_n value_n\n";
		}
		$task_string .=
"\$task_description .= ' -$variableStr \"'.join( '\" \"', \@$variableStr ).'\"' if ( defined \$$variableStr"
		  . "[0]);\n";
		$error_check .=
		    "unless ( defined \$$variableStr" . "[0]" . ") {\n"
		  . "	\$error .= \"the cmd line switch -$variableStr is undefined!\\n\";\n"
		  . "}\n";
		$optionNames->{$variableStr} = 'ARRAY';
	}
	elsif ( $variableStr =~ s/#option$// ){
		$add_2_variable_def .= ", \$$variableStr";
		$add_2_variable_read .=
		  "       \"-$variableStr\"    => \\\$$variableStr,\n";
		$add_2_help_string .=
		  "       -$variableStr       :<please add some info!>\n";
		
		$task_string .=
"\$task_description .= \" -$variableStr \" if ( \$$variableStr);\n";
		$error_check .= "# $variableStr - no checks necessary\n";
		$optionNames->{$variableStr} = 'OPTION';
	}
	else {
		$add_2_variable_def .= ", \$$variableStr";
		$add_2_variable_read .=
		  "	 \"-$variableStr=s\"    => \\\$$variableStr,\n";
		$add_2_help_string .=
		  "       -$variableStr       :<please add some info!>\n";
		$task_string .=
"\$task_description .= \" -$variableStr '\$$variableStr'\" if (defined \$$variableStr);\n";
		$error_check .=
		    "unless ( defined \$$variableStr) {\n"
		  . "	\$error .= \"the cmd line switch -$variableStr is undefined!\\n\";\n"
		  . "}\n";
		$optionNames->{$variableStr} = 'STRING';
	}

}
my $V =stefans_libs::Version->new();
my $string = "#! /usr/bin/perl -w

=head1 LICENCE

  Copyright (C) " . root::Today() . " Stefan Lang

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
   
   binCreate.pl from ".$V->origin('Stefans_Lib_Esentials')." commit ".$V->version('Stefans_Lib_Esentials')."
   

=head1  SYNOPSIS

    EXECUTABLE
$add_2_help_string

       -help           :print this help
       -debug          :verbose output
   
=head1 DESCRIPTION

  INFO_STR

  To get further help use 'EXECUTABLE -help' at the comman line.

=cut

use Getopt::Long;
use Pod::Usage;

use strict;
use warnings;

use FindBin;
my \$plugin_path = \"\$FindBin::Bin\";

my \$VERSION = 'v1.0';


my ( \$help, \$debug, \$database$add_2_variable_def);

Getopt::Long::GetOptions(
$add_2_variable_read
	 \"-help\"             => \\\$help,
	 \"-debug\"            => \\\$debug
);

my \$warn = '';
my \$error = '';

$error_check

if ( \$help ){
	print helpString( ) ;
	exit;
}

if ( \$error =~ m/\\w/ ){
	helpString(\$error ) ;
	exit;
}

sub helpString {
	my \$errorMessage = shift;
	\$errorMessage = ' ' unless ( defined \$errorMessage); 
	print \"\$errorMessage.\\n\";
	pod2usage(q(-verbose) => 1);
}

$additional_checks

my ( \$task_description);

$task_string

$options_string
$log_str
## Do whatever you want!

";

$string =~ s/EXECUTABLE/$exec_name/g;
$string =~ s/INFO_STR/$pod/g;

eval {
&createTestFile( "$path/../t/", $exec_name, $optionNames );
};

unless ( -f "$path/$exec_name" && !$force ) {
	open( OUT, ">$path/$exec_name" )
	  or die "could not create file '$path/$exec_name'\n";
	print OUT $string;
	close OUT;

	print "\nnew executable written to '$path/$exec_name'\n\n";
}
else {
	warn
"\nthe file '$path/$exec_name' already exists!\nUse the -force switch to delete it\n\n";
}


sub helpString {
	my $errorMessage = shift;
	$errorMessage = ' ' unless ( defined $errorMessage );
	return "
 $errorMessage
 command line switches for binCreate.pl
   -pod            :a small text, that describes the function of the script (required)
   -name           :the name of the script as a full filename to the position of the script
   -force          :delete an existing script (no warning!)
   -commandLineSwitches
                   :a list of option switches that you want to have included into that script
   -database       :the databse to use for logging
   -help           :print this help
   -debug          :verbose output
 ";
}

sub createTestFile {
	my ( $testPath, $executable, $optionNames ) = @_;
	$executable =~ s/.pl$//;
	if ( -f "$testPath/xx_$executable.t" and !$force ) {
		print "test file is already present ($testPath/xx_$executable.t)\n";
		return "$testPath/xx_$executable.t";
	}
	open( Test, ">$testPath/xx_$executable.t" )
	  or die "could not open testFile $testPath/xx_$executable.t\n";

	## we have to create a stup for each sub in the INFILE!
	print Test "#! /usr/bin/perl
use strict;
use warnings;
use stefans_libs::root;
use Test::More tests => 2;
use stefans_libs::flexible_data_structures::data_table;

use FindBin;
my \$plugin_path = \"\$FindBin::Bin\";

my ( \$value, \@values, \$exp, ";
	foreach ( keys %$optionNames ) {
		if ( $optionNames->{$_} eq "ARRAY" ) {
			print Test "\@$_, ";
		}elsif ( $optionNames->{$_} eq "OPTION"){
			## do not need that!
		}
		else {
			print Test "\$$_, ";
		}
	}
	print Test ");\n
my \$exec = \$plugin_path . \"/../bin/$executable.pl\";
ok( -f \$exec, 'the script has been found' );
my \$outpath = \"\$plugin_path/data/output/$executable\";
if ( -d \$outpath ) {
	system(\"rm -Rf \$outpath\");
}


my \$cmd =
    \"perl -I \$plugin_path/../lib  \$exec \"\n";
	foreach ( keys %$optionNames ) {
		if ( $optionNames->{$_} eq "ARRAY" ) {
			print Test ". \" -$_ \" . join(' ', \@$_ )\n";
		}elsif ( $optionNames->{$_} eq "OPTION"){
			print Test ". \" -$_ \" # or not?\n";
		}
		else {
			print Test ". \" -$_ \" . \$$_ \n";
		}
	}
	print Test ". \" -debug\";\n". "my \$start = time;\n"
	  . "system( \$cmd );\n"
	  . "my \$duration = time - \$start;\n"
	  . "print \"Execution time: \$duration s\\n\";\n"
	  . "#print \"\\\$exp = \".root->print_perl_var_def(\$value ).\";\\n\";";

	close(OUT);
	return "$testPath/xx_$executable.t";
}
