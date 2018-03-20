#! /usr/bin/perl
use strict;
use warnings;
use stefans_libs::root;
use Test::More tests => 10;
BEGIN { use_ok 'stefans_libs::SLURM' }

use FindBin;
my $plugin_path = $FindBin::Bin;

my ( $value, @values, $exp );


my $outpath = "$plugin_path/data/outpath/stefans_libs_SLURM";
if ( -d $outpath ) {
	system("rm -Rf $outpath/*");
}else {
	mkdir ( $outpath );
}

chdir($outpath);

$value = { n => 2, N => 1, t => '01:00:00', A => 'snic2016-4-13' , w => 'nodelist=ls2-n3'};
my $SLURM = stefans_libs::SLURM->new( $value );

is_deeply($SLURM->{'options'}->options(), { n => 2, N => 1, t => '01:00:00', A => 'snic2016-4-13', w => 'nodelist=ls2-n3'  },"options stored correctly" );

ok( ref($SLURM) eq 'stefans_libs::SLURM',
	'simple test of function stefans_libs::SLURM -> new()' );

#is_deeply( $value, {}, "SLURM variables removed from hash");

$value = $SLURM->script( "some_command.pl", "my_name" );

#print "\$exp =  " . root->print_perl_var_def( [ split( "\n", $value ) ] ) . ";\n";
$exp = [
	'#! /bin/bash',
	'#SBATCH -w, --nodelist=ls2-n3',
	'#SBATCH -n 2',
	'#SBATCH -N 1',
	'#SBATCH -t 01:00:00',
	'#SBATCH -A snic2016-4-13',
	'#SBATCH -J my_name',
	'#SBATCH -o my_name%j.out',
	'#SBATCH -e my_name%j.err',
	'some_command.pl'
];

is_deeply([ split( "\n", $value ) ], $exp, "right script created" );
$SLURM->{'debug'} = 1;
$value = root->filemap( $plugin_path."/data/outpath/stefans_libs_SLURM.a.b.c.sh" );
unlink( $value->{'total'} ) if ( -f  $value->{'total'} );
$SLURM-> run ( "some_command.pl", $value );

ok ( -f $value->{'total'}, 'script file created'. " ($value->{'total'})" );
open ( IN, "<$value->{'total'}" );
$value = [];
@values = <IN>;
close ( IN );
$value = [ map { chomp; $_ } @values ];
$exp = [
	'#! /bin/bash',
	'#SBATCH -w, --nodelist=ls2-n3',
	'#SBATCH -n 2',
	'#SBATCH -N 1',
	'#SBATCH -t 01:00:00',
	'#SBATCH -A snic2016-4-13',
	'#SBATCH -J stefans_libs_SLURM.a.b.c',
	'#SBATCH -o stefans_libs_SLURM.a.b.c%j.out',
	'#SBATCH -e stefans_libs_SLURM.a.b.c%j.err',
	'some_command.pl'
];
is_deeply( $exp, $value, "the script file is as expected" );

#$SLURM->{'options'}->save();

$value = $SLURM->script( "echo 'this is fun'\nfailThis\ndf -h > $outpath/df_output.txt", "run_df" );

$exp = [
	'#! /bin/bash',
	'#SBATCH -w, --nodelist=ls2-n3',
	'#SBATCH -n 2',
	'#SBATCH -N 1',
	'#SBATCH -t 01:00:00',
	'#SBATCH -A snic2016-4-13',
	'#SBATCH -J run_df',
	'#SBATCH -o run_df%j.out',
	'#SBATCH -e run_df%j.err',
	"echo 'this is fun'",
	'failThis',
	"df -h > $outpath/df_output.txt"
];

is_deeply( [split("\n",$value)], $exp, "the script file is as expected" );

## OK now test the local run:
$SLURM->{'run_local'} = 1; ## I do not want to test the SLURM environment here ;-)
$SLURM->{'debug'} = 0;
$SLURM-> run( "df -h > $outpath/df_output.txt", "$outpath/run_df" );

ok ( -f "$outpath/run_df.sh" , "script created" );
ok ( -f "$outpath/df_output.txt" , "output created" );

opendir( DIR, $outpath ) or die $!;
my @slurm_outfiles = map { $_ =~ s/local\.\d+\./local.17111./; $_ }  grep { /run_df.*[eo][ru][rt]/ } readdir( DIR );
closedir( DIR );

is_deeply ( \@slurm_outfiles, ['run_df.local.17111.out', 'run_df.local.17111.err'], "No SLURM output created - as wanted" );


#print "\$exp = ".root->print_perl_var_def($value ).";\n";

