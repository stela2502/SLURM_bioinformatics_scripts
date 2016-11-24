#! /usr/bin/perl -w

=head1 LICENCE

  Copyright (C) 2016-11-23 Stefan Lang

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

    reformat_bed_gtf_intersect.pl
       -infile       :<please add some info!>
       -names       :<please add some info!>
       -options     :<please add some info!> you can specify more entries to that
                         format: key_1 value_1 key_2 value_2 ... key_n value_n
       -outfile       :<please add some info!>


       -help           :print this help
       -debug          :verbose output
   
=head1 DESCRIPTION

  Reformate a file created by the annotate_bed_with_gtf_genome.pl perl script.

  To get further help use 'reformat_bed_gtf_intersect.pl -help' at the comman line.

=cut

use Getopt::Long;
use Pod::Usage;

use stefans_libs::flexible_data_structures::data_table;

use strict;
use warnings;

use FindBin;
my $plugin_path = "$FindBin::Bin";

my $VERSION = 'v1.0';


my ( $help, $debug, $database, $infile, $names, $options, @options, $outfile);

Getopt::Long::GetOptions(
	 "-infile=s"    => \$infile,
       "-names"    => \$names,
       "-options=s{,}"    => \@options,
	 "-outfile=s"    => \$outfile,

	 "-help"             => \$help,
	 "-debug"            => \$debug
);

my $warn = '';
my $error = '';

unless ( defined $infile) {
	$error .= "the cmd line switch -infile is undefined!\n";
}
# names - no checks necessary
unless ( defined $options[0]) {
	$error .= "the cmd line switch -options is undefined!\n";
}
unless ( defined $outfile) {
	$error .= "the cmd line switch -outfile is undefined!\n";
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

$task_description .= 'perl '.$plugin_path .'/reformat_bed_gtf_intersect.pl';
$task_description .= " -infile '$infile'" if (defined $infile);
$task_description .= " -names " if ( $names);
$task_description .= ' -options "'.join( '" "', @options ).'"' if ( defined $options[0]);
$task_description .= " -outfile '$outfile'" if (defined $outfile);


for ( my $i = 0 ; $i < @options ; $i += 2 ) {
	$options[ $i + 1 ] =~ s/\n/ /g;
	$options->{ $options[$i] } = $options[ $i + 1 ];
}
###### default options ########
#$options->{'something'} ||= 'default value';
##############################

my $fm = root->filemap( $outfile );
unless ( -d $fm->{'path'} ){
	system( "mkdir -p $fm->{'path'}");
}
open ( LOG , ">$outfile.log") or die $!;
print LOG $task_description."\n";
close ( LOG );


## Do whatever you want!

#chromosome    start    end    name    orig.bed.file.col.4    orig.bed.file.col.5    orig.bed.file.col.6    orig.bed.file.col.7    

#gtf_seqname    gtf_source    gtf_feature    gtf_start    gtf_end    gtf_score    gtf_strand    gtf_frame    gtf_attribute    gtf_gene_id    gtf_transcript_id    
#gtf_gene_type    gtf_gene_status    gtf_gene_name    gtf_transcript_type    gtf_transcript_status    gtf_transcript_name    gtf_exon_number    gtf_exon_id    gtf_level    gtf_protein_id    gtf_tag    gtf_transcript_support_level    gtf_ccdsid    gtf_havana_gene    gtf_havana_transcript    gtf_ont



#chr1    93847572    93847608    X    301.667    +    0.000296039    0.231056    

#gtf_seqname   #chr1 // chr1 // chr1 // chr1    
#gtf_source    #HAVANA // HAVANA // HAVANA // MANUAL    
#gtf_feature   #exon // gene // transcript // gene    
#gtf_start     #93847174 // 93847174 // 93847174 // 93847572    
#gtf_end       #93848939 // 93848939 // 93848939 // 93847657    
#gtf_score     #. // . // . // .    
#gtf_strand    #+ // + // + // +    
#gtf_frame     #. // . // . // .    
#gtf_attribute #ENSG00000260464.1 // ENSG00000260464.1 // ENSG00000260464.1 // hg38:chr1-93847572:93847657    
#gtf_gene_id   #ENST00000565336.1 //  // ENST00000565336.1 //     
#gtf_transcript_id    #lincRNA // lincRNA // lincRNA // tRNA    
#gtf_gene_type #KNOWN // KNOWN // KNOWN // KNOWN    
#gtf_gene_status      #RP4-561L24.3 // RP4-561L24.3 // RP4-561L24.3 // tRNA-Arg-AGA    
#gtf_gene_name #lincRNA //  // lincRNA //     
#gtf_transcript_type  #KNOWN //  // KNOWN //     
#gtf_transcript_status #RP4-561L24.3-001 //  // RP4-561L24.3-001 //     
#gtf_transcript_name#1 //  //  //     
##ENSE00002588542.1 //  //  //     
##2 // 2 // 2 //      
##//  //  //     
##basic //  // basic //     
##NA //  // NA //      
##//  //  //     
##OTTHUMG00000175883.1 // OTTHUMG00000175883.1 // OTTHUMG00000175883.1 //     
##OTTHUMT00000431234.1 //  // OTTHUMT00000431234.1 //      
##//  //  //      
##//  //  //


#chr1    145459620    145459656    X    102    +    1.1625e-12    0.346579    

#chr1    MANUAL    gene    145459657    145459729    .+    .    hg38:chr1-145459657:145459729        tRNA    KNOWN    tRNA-Gln-CAG
    

my ( @line );
open ( IN , "<$infile" ) or die "could not open infile '$infile'\n$!\n";

my $data_table = data_table->new();
my ( @orig_cols_names, @orig_col_ids, @gtf_col_name, @gtf_ids,  $hash, $tmp_table );
while ( <IN> ) {
	chomp($_);
	@line = split ( "\t", $_ );
    if ( $_ =~ m/^#/ ) { ## is header line
    	$line[0] =~ s/^#//;
    	my $colname;
    	for (my $i = 0; $i < @line; $i ++ ) {
    		$colname = $line[$i];
    		if ( $colname =~ m/^gtf_/ ){
    			push ( @gtf_ids , $i );
    			push ( @orig_cols_names, $colname);
    		}else {
    			$data_table->Add_2_Header( $colname );
    			push ( @orig_col_ids , $i );
    			push ( @orig_cols_names, $colname);
    		}
    	}
    	next;
    }
    ## not a header line
    $hash = {};
    map { $hash->{$orig_cols_names[$_]} = $line[$_] } @orig_col_ids;
    $tmp_table = data_table->new();
    
    foreach ( @gtf_ids ) {
    	unless ( defined  $line[$_]){
    		$line[$_] = '' 
    	}else {
    	if ( $line[$_] =~ m/^ \/\// ) {
    		$line[$_] = "---".$line[$_];
    	}
    	if ( $line[$_] =~ m/ \/\/ $/ ){
    		 $line[$_] .= "---";
    	}
    	}
    	my @values = map { if ( defined $_ ) {$_} else{ '---' } } split( " // ",$line[$_]);
    	if ( @values == 0 ){
    		@values = ( '' );
    	}
    #	print "$orig_cols_names[$_]\t".join("\t", @values)."\n";
    	$tmp_table->add_column( $orig_cols_names[$_] , @values );
    }
    
    $tmp_table->write_table( $outfile );

    #print "\$exp = ".root->print_perl_var_def({ 'hash' => $hash, 'orig_col_ids' => \@orig_col_ids, 'gtf_ids' => \@gtf_ids, 'gtf_col_name' => \@gtf_col_name }).";\n";
	#die "does that make sense??\n".$tmp_table->AsString();
	## we want to process that table into a 
	#gene	trna	rna exon
	#1	0	1	1
	#like structure
}





