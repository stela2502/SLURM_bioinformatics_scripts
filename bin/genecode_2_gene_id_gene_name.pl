#! /usr7bin/perl -w

use strict;
use warnings;

#problems:
#gencode.vM9.annotation.gff3:chrM	ENSEMBL	exon	15356	15422	.	-	.	ID=exon:ENSMUST00000082423.1:1;Parent=ENSMUST00000082423.1;gene_id=ENSMUSG00000064372.1;transcript_id=ENSMUST00000082423.1;gene_type=Mt_tRNA;gene_status=KNOWN;gene_name=mt-Tp;transcript_type=Mt_tRNA;transcript_status=KNOWN;transcript_name=mt-Tp-201;exon_number=1;exon_id=ENSMUSE00000521550.1;level=3;transcript_support_level=NA;tag=basic

my ($geneid, $seen);
print "gene_id\tgene_name\tgene_type\n";
while(my $line = <> ) { 
	if (  $geneid = &match_and_ret($line,'gene_id') ){ #.*gene_name=([\w\.\_\-]+);/) {
		next if ( $seen -> { $geneid } );
		$seen -> { $geneid } =1;
		print "$geneid";
		foreach ( 'gene_name', 'gene_type' ) {
			print "\t".&match_and_ret($line,$_);
		}
		print "\n";
	} 
}


sub match_and_ret{
	my ($str, $tag) = @_;
	# gene_id "ENSMUSG00000102693.1";
	if ( $str =~ m/$tag=([\w\.\-]+);/){ 
		return $1;
	}
	return '';
}
