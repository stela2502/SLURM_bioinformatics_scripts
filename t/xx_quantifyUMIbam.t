#! /usr/bin/perl
use strict;
use warnings;
use stefans_libs::root;
use Test::More tests => 17;
use stefans_libs::flexible_data_structures::data_table;
use Digest::MD5;


use FindBin;
my $plugin_path = "$FindBin::Bin";

my ( $value, @values, $exp, $organism, $version, $genome_bed_file, @bam_files, $outpath, );

my $exec = $plugin_path . "/../bin/quantifyUMIbam.pl";
ok( -f $exec, 'the script has been found' );
$outpath = "$plugin_path/data/output/quantifyUMIbam";
if ( -d $outpath ) {
	system("rm -Rf $outpath/*.log $outpath/*.xls" );
}else {
	system( "mkdir $outpath" );
}
foreach ( 'AGATAG_top_1000000.sam.gz', 'CAGATA_900000.sam.gz' ) {
	ok ( -f $plugin_path ."/data/$_", "data file $_" );
	$value = $_;
	$value =~ s/.gz$//;
	unless ( -f "$outpath/$value" ) {
		system( "cp $plugin_path/data/$_  $outpath " );
		system( "gunzip $outpath/$_");
	}
	push ( @bam_files, "$outpath/$value" );
}

$version = 'BUILD.38.1';
$organism = 'M_musculus';
$genome_bed_file = '';
my $cmd =
    "perl -I $plugin_path/../lib -I $plugin_path/../../stefans_libs-GenomeDB/lib $exec "
. " -organism " . $organism 
. " -version " . $version 
#. " -genome_bed_file " . $genome_bed_file 
. " -bam_files " . join(' ', @bam_files )
. " -outpath " . $outpath 
. " -debug";
print ( $cmd."\n" );

my $start_run = time();

system ( $cmd );

my $end_run = time();
my $run_time = $end_run - $start_run;

print "Job took $run_time seconds\n";


my @digests = (
'a419625558b322643cf41fce16ff2a70',  'AGATAG_top_1000000.sam',
'6a424b16bc86d4904c2466a9d1116820',  'AGATAG_top_1000000.simple',
'2026edecd352dfebe9be14426da36ddb',  'AGATAG_top_1000000.txt',
'9e55f03b4340bc6536c564b0f1180dbc',  'CAGATA_900000.sam',
'd8ac181758deb5176f941302dc08329a',  'CAGATA_900000.simple',
'054727d186dea03f311b93be737e9365',  'CAGATA_900000.txt',
'c4a8822da307839b29cb2c9b31134a2b',  'genome.bed'
);

foreach ( qw(  AGATAG_top_1000000.simple  CAGATA_900000.sam     CAGATA_900000.txt
AGATAG_top_1000000.sam      AGATAG_top_1000000.txt     CAGATA_900000.simple  genome.bed ) ) {
	ok ( -f $outpath."/".$_, "outfile $_" );
	if ( $_ =~ m/txt$/ ) {
		system( "sort $outpath/$_ > $outpath/$_.sorted" );
		system( "mv $outpath/$_.sorted $outpath/$_");
	}
}

my ($ctx, $filehandle, $tmp);
for ( my $i = 0; $i < @digests; $i +=2 ) {
	$ctx = Digest::MD5->new;
	open ( $filehandle, "<$outpath/$digests[$i+1]" );
	$ctx->addfile($filehandle);
	close ( $filehandle );
	ok( $digests[$i] eq  ($tmp = $ctx->hexdigest), "outfile $digests[$i+1] ($tmp)");
}


#print "\$exp = ".root->print_perl_var_def($value ).";\n