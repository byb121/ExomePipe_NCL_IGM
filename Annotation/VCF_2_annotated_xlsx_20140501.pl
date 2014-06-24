#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use Math::CDF;
use Excel::Writer::XLSX;
use Spreadsheet::ParseExcel;

print "\n";
print "VCF_2_annotated_excel_20140219.pl says\n";
print "############################################################################################\n";
print "# YX: in house annotation part was amended from HG's                                       #\n";
print "# YX: the script will invoke Annovar and produce annotated output in Excel XLSX format     #\n";
print "# YX: input shoud be a single vcf file of a family or single sample                        #\n";                                 
print "# HG edits: to take (standalone) annovar output (with vcf file lines appended after)       #\n";
print "# vcf output (i.e. chr) starts at column 37 (fileline array element 36)                    #\n";
print "# column number of 'vcf' chr will depend on how many annotations are included from annovar #\n";
print "############################################################################################\n";
print "\n";
print "##########################################################################################################\n";
print "# Changes in this version:                                                                               #\n";
print "# 1. Minor change of the last version: will output \'Everything\' table into a separate file if request, #\n";
print "#                                                                       check help!                      #\n";
print "# 2. Minor change of the this version: extra sample column to flag low qual genotype calls.              #\n";
print "# 3. Minor change of the this version: gene_type_call_qual raised to 20                                  #\n";
print "##########################################################################################################\n";

my $vcf_in;
my $InterestedGenefile=""; #Ensembl gene ids on each line
my $exomiserXLSfiles="";
my $CNV_file="";
my $output_excel=""; #if empty then inHouse annotated file will not converted to excel file
my $output_excel_everything=""; # file to output table 'everything'
my $add_genotypeCall_flags="No";
my $annovarDIR;
my $help;


usage() if ( @ARGV < 1 || ! GetOptions('help|?' => \$help, "vcf=s"=>\$vcf_in, "InterestedGenes=s"=>\$InterestedGenefile, 'CNV=s' => \$CNV_file, 
'exomiserXLS=s' => \$exomiserXLSfiles, 'out=s' => \$output_excel, 'outAll=s' => \$output_excel_everything, 'AnnovarDIR=s' => \$annovarDIR,
		'add_genotypeCall_flags=s' => \$add_genotypeCall_flags) || defined $help );

unless (defined $vcf_in && -e $vcf_in) {
	die "You have not supplied input vcf file using --vcf or the file does not exist\n\n";
}

unless (defined $annovarDIR && -e $vcf_in) {
        die "You have not supplied the valid directory where Annovar is installed and annovar_YX.sh is located, use --AnnovarDIR to set it up.\n\n";
}

print "Start to run Annovar on $vcf_in\n";
`$annovarDIR/annovar_YX.sh $vcf_in`;
print "Annovar is done!\n\n";

my $annovar_output_filename = $vcf_in.".hg19_multianno.txt";
my $vcf = $annovar_output_filename;
my $output_file=$annovar_output_filename."_inHouse_annotated.txt"; #in house annotation script output file

my $input_file;
if ($output_excel ne "") {$input_file=$output_file;} #input file for the script part to convert it to excel file

#### input part for converting excel file part
my $bash_output=`grep '#CHR' $vcf_in`;
chomp $bash_output;
my @sample_names;
my @sample_columns;
my @bash_output_split = split("\t", $bash_output);
my $first_sample_column_index = 47; # change if add or reduced Annovar annotation 0_based
# it is the sample column not the genotype column eg: 0/0:..... not the R/R column

print "Looking for sample names:\n";
for(my $i=9;$i<scalar @bash_output_split;$i++){
	print "Found: ".$bash_output_split[$i]."\n";
	push @sample_names, $bash_output_split[$i];
}

for(my $i=1;$i<=scalar @sample_names;$i++){
	push @sample_columns, $first_sample_column_index;
	$first_sample_column_index += 2;
}


#### Requiring input in house db files to add as annotation
my $genefile="/users/data/Files_HG/vcf_annotation_november2013/inHouse_db/Ensembl_Genes_GRCh37_75.txt";
my $sharedfile="/users/data/Files_HG/vcf_annotation_november2013/inHouse_db/InHouse_OnTarget_Variants_MAFs.txt_418exomes";
my $sharedfile_GATK="/users/data/Files_HG/vcf_annotation_november2013/inHouse_db/InHouse_OnTarget_GATK_MAFs.txt_179exomes";
my $OMIMfile="/users/data/Files_HG/vcf_annotation_november2013/inHouse_db/Ensembl_GRCh37_75_OMIM_AllGeneNames.txt";
my $ens_gene_OMIM_Uniprot_Acc = "/users/data/Files_HG/vcf_annotation_november2013/inHouse_db/ens_gene_symbol_omim_id_uniprot_id.txt";

print "\n";
print "Require db file: Ensembl_Genes_GRCh37_75.txt\n";
if (-e $genefile) {
	print "Found $genefile\n";
	print "Consider update when necessary\n";
} else {
	print "$genefile not exist\n exit\n";
	exit;
}

print "Require db file: InHouse_OnTarget_Variants_MAFs.txt\n";
if (-e $sharedfile) {
	print "Found $sharedfile\n";
	print "Consider update when necessary\n";
} else {
	print "$sharedfile not exist\n exit\n";
	exit;
}

print "Require db file: InHouse_OnTarget_GATK_MAFs.txt\n";
if (-e $sharedfile_GATK) {
	print "Found $sharedfile_GATK\n";
	print "Consider update when necessary\n";
} else {
	print "$sharedfile_GATK not exist\n exit\n";
	exit;
}

print "Require db file: Ensembl_GRCh37_75_OMIM_AllGeneNames.txt\n";
if (-e $OMIMfile) {
	print "Found $OMIMfile\n";
	print "Consider update when necessary\n";
} else {
	print "$OMIMfile not exist\n exit\n";
	exit;
}

print "Require db file: ens_gene_symbol_omim_id_uniprot_id.txt\n";
if (-e $ens_gene_OMIM_Uniprot_Acc) {
	print "Found $ens_gene_OMIM_Uniprot_Acc\n";
	print "Consider update when necessary\n";
} else {
	print "$ens_gene_OMIM_Uniprot_Acc not exist\n exit\n";
	exit;
}

# An comlumn to tell if the variant is a SNP or INDEL
my $SNPorINDEL; # this will add a column to indicate the variant is a SNP/indel; easy for filtering

### Read in all annotations first ###
my %en_genes;
my %ens_symbol_map;
my %inHouse_MAF;
my %inHouse_MAF_GATK;
my %OMIM;
my %isInterestedGenes;
my %OmimAcc;
my %UniprotAcc;

my $en_genes_ref = GetGeneCoords($genefile);
%en_genes = %$en_genes_ref;
my $ens_symbol_map_ref = GetEnsSymbolMap($genefile);
%ens_symbol_map = %$ens_symbol_map_ref;
my $inHouse_MAF_ref = GetInHouseMaf($sharedfile);
%inHouse_MAF = %$inHouse_MAF_ref;
my $inHouse_MAF_GATK_ref = GetInHouseMafGATK($sharedfile_GATK);
%inHouse_MAF_GATK = %$inHouse_MAF_GATK_ref;
my $OMIM_ref = GetOMIManno($OMIMfile);
%OMIM = %$OMIM_ref;
if ($InterestedGenefile ne "") {
	my $isInterestedGenes_ref = GetIsInterestedGenes($InterestedGenefile);
	%isInterestedGenes = %$isInterestedGenes_ref;
}

my $OmimAcc_ref = GetOmimAcc($ens_gene_OMIM_Uniprot_Acc);
%OmimAcc = %$OmimAcc_ref;
my $UniprotAcc_ref = GetUniprotAcc($ens_gene_OMIM_Uniprot_Acc);
%UniprotAcc = %$UniprotAcc_ref;

my @output;

open VCF, "$vcf" or die "Can not open the file $vcf";
while (my $line = <VCF> ) {
	chomp $line;
	if($line =~ m/^Chr/) {
		push @output, $line."\n";
		next;
	} else {
		my @elements  = split ("\t", $line);

		my $annovar_chr=$elements[0];
		my $annovar_pos=$elements[1];
		my $annovar_R=$elements[3];
		my $annovar_A=$elements[4];
		
		my $vcf_chr=$elements[38];
		my $vcf_pos=$elements[39];
		my $vcf_R=$elements[41];
		my @temp_temp_temp=split(/\,/, $elements[42]);
		my $vcf_A=$temp_temp_temp[0];
		 		
		my $FORMAT=$elements[46];
		my $Sample_Call="";
		for(my $i=47;$i<scalar(@elements);$i++) {
			$Sample_Call=$Sample_Call."\t".$elements[$i];
		}
		$Sample_Call =~ s/^\t//;
		
		
		if ($annovar_R =~ m/^[atgcATGC]$/ && $annovar_A =~ m/^[atgcATGC]$/) { # a SNP with just one alternative allele
			$SNPorINDEL = "SNP"; 
		} else {
			$SNPorINDEL="INDEL";
		}
		
		my $genecard_link = "";
		my $gene_name;
		my $ens_id;
		
		if ($elements[5] eq "intergenic") {
			$ens_id = "NA";
		} elsif ($elements[5] eq "splicing") {
			if ($elements[6] =~ m/^(.*?)\(.*\)$/) {
				$ens_id = $1;
			} elsif ($elements[6] !~ m/(\(|\,)/) {
				$ens_id = $elements[6];
			} elsif ($elements[6] =~ m/(EN\w+\d+?)\,EN.*/) {
				$ens_id = $1;
			} else {
				print "fatal match error on ensembl annotation splicing match;\n";
				print "$line.\n\n";
				exit 2;
			}
		} elsif ($elements[5] eq 'exonic;splicing') {
			if ($elements[6] =~ m/^(.*?)\;/) {
				$ens_id = $1;
			} else {
				$ens_id = $elements[6];
			}
		} else {
			if ($elements[6] ne "NA" &&  $elements[6] !~ m/(\(|\,)/) {
				$ens_id = $elements[6];	
			} else {
				$ens_id = "NA";
			}
		}
		
		if ($ens_id eq "NA") {
			$genecard_link = "NA";
		} else {
			$genecard_link = "=HYPERLINK(\'http://www.genecards.org/cgi-bin/carddisp.pl?id=$ens_id&id_type=ensembl\', \'GeneCard Link\')";
		}
		
		if (exists $ens_symbol_map{$ens_id}) {
			$gene_name = $ens_symbol_map{$ens_id};
		} else {
			$gene_name = "NA";
		}
		
			
		# add OMIM links and Uniprot links
		############# all links added need to have single quotes replaced bu double quotes after wAnnovar ############
		my $omim_link = "";
		if (exists $OmimAcc{$ens_id}) {
			my $omim_accesion =  $OmimAcc{$ens_id};
			$omim_link = "=HYPERLINK(\'http://omim.org/entry/$omim_accesion\', \'OMIM Link\')";
		} else {
			$omim_link = "NA";
		}
		
		my $uniprot_link = "";
		if (exists $UniprotAcc{$ens_id}) {
			my $uniprot_accession = $UniprotAcc{$ens_id};
			$uniprot_link = "=HYPERLINK(\'http://www.uniprot.org/uniprot/$uniprot_accession\', \'Uniprot Link\')";
		} else {
			$uniprot_link = "NA";
		}
		
		# find if the gene is in the interested gene list
		my $isInterested = "NO";
		if (exists $isInterestedGenes{$ens_id}) {
			$isInterested = "YES";
		}
		
		#find the corresponding omim annotation
		my $omim_anno = "NA\tNA\tNA\tNA";
		if (exists $OMIM{$ens_id}) {
			$omim_anno = $OMIM{$ens_id};
		}		
		
		############################### find in house MAF for variants
		my $maf = "";
		my $pos1;
		my $ref1;
		my $alt1;
		if ($vcf_R =~ m/^[atgcATGC]$/ && $vcf_A =~ m/^[atgcATGC]$/) { # a SNP with just one alternative allele
			$pos1 = $vcf_pos;
			$ref1 = $vcf_R;
			$alt1 = $vcf_A;
		} elsif ($vcf_R =~ m/^[atgcATGC]$/) { # insertions
			$pos1 = $vcf_pos;
			$ref1 = $vcf_R;
			$alt1 = $vcf_A;
			$alt1 =~ s/^./\+/;
		} else { ##deletion(s) !!and in some cases of multiple (,) vars insertions!!
			if 	(length $vcf_A < length $vcf_R){
				my $diff=(length $vcf_A) - 1; 
				$alt1 = substr($vcf_R,$diff);
				$ref1 = substr($vcf_R,$diff, 1);
				$alt1 =~s/^./\-/; 
				$pos1 = $vcf_pos+$diff;
			} elsif(length $vcf_A > length $vcf_R){
				my $diff=(length $vcf_R)-1; 
				$alt1 = substr($vcf_A,$diff);
				$ref1 = substr($vcf_A,$diff,1); 
				$alt1 =~s/^./\+/;
				$pos1 = $vcf_pos+$diff;			
			} else {
				$alt1 =  substr($vcf_A, 0, 1);
				$ref1 = substr($vcf_R, 0, 1);
				$pos1 = $vcf_pos;
			}
		}
		
		if (exists $inHouse_MAF{$vcf_chr}{$pos1}{$ref1}{$alt1} ) {
			$maf = $inHouse_MAF{$vcf_chr}{$pos1}{$ref1}{$alt1};
		} else {
			$maf = "NA";
		}

		#find in house 'GATK' MAF for variants
		############################################## need to map v with Annovar format
		my $maf_GATK = "";
		if (exists $inHouse_MAF_GATK{$annovar_chr}{$annovar_pos}{$annovar_R}{$annovar_A}) {
			$maf_GATK = $inHouse_MAF_GATK{$annovar_chr}{$annovar_pos}{$annovar_R}{$annovar_A};
		} else {
			$maf_GATK = "NA";
		}
		############################################## vcf format changed done #######
		
		my $Sample_Call_processed = AddGenoTypeToSampleCalls($FORMAT, $Sample_Call);
		for (my $i=0;$i<=46;$i++){ push @output, $elements[$i]."\t";}
		push @output, $Sample_Call_processed."\t".$SNPorINDEL ."\t".$gene_name."\t".$ens_id."\t".$isInterested."\t"
		.$maf."\t".$maf_GATK."\t".$omim_anno."\t".$genecard_link."\t".$omim_link."\t".$uniprot_link."\n";
	}
}
close (VCF);

open OUTPUT, ">$output_file" or die "Cannot open file $output_file to output. \n";
print OUTPUT @output;
close OUTPUT;



################################################
#
#  output to excel part
#
#


if ($output_excel eq "") {print "#####  No excel file is required. ######\nDone!\n"; exit;}

if(scalar @sample_columns != scalar @sample_names) {
	die "Number of sample names and columns are inconsistent.\n";
}
my @output_everything;
my @output_all;
my @output_filtered;
my @output_possibleHits;
my @very_rare_interested_vars;  #to filter very rare < 0.01 variants on interested genes
my @output_Xlinked;
my @output_ARmodelHits;
my @output_CNV;
my @output_exp_compound;
my %CNV_regions;
#my %gene_name_count; #it's for selecting Xlinked rare variants
my $ArrayOfExomiserXLS_ref;
my @ArrayOfExomiserXLS;
my %HashOfExomiserXLS;
my @ArrayOfExomiserXLS_Output;

if ($exomiserXLSfiles ne "") {
	my @exomiserFileList = split(",", $exomiserXLSfiles);
	foreach my $exomiserFile (@exomiserFileList){
		$ArrayOfExomiserXLS_ref = ExomiserXLS2Array ($exomiserFile);
		push  @ArrayOfExomiserXLS, @{$ArrayOfExomiserXLS_ref};
	}
	
	my %HashofExomiserRecordCount;
	foreach my $exomiserRecord (@ArrayOfExomiserXLS) {
		my @temp = split("\t", $exomiserRecord);
		if (!exists $HashofExomiserRecordCount{$temp[0]}{$temp[1]}{$temp[2]}{$temp[3]}) {
			$HashofExomiserRecordCount{$temp[0]}{$temp[1]}{$temp[2]}{$temp[3]} = 1;
		} else {
			$HashofExomiserRecordCount{$temp[0]}{$temp[1]}{$temp[2]}{$temp[3]} = $HashofExomiserRecordCount{$temp[0]}{$temp[1]}{$temp[2]}{$temp[3]} + 1;
		}
	}
	
	foreach my $exomiserRecord (@ArrayOfExomiserXLS) {
		my @temp = split("\t", $exomiserRecord);
		if ($HashofExomiserRecordCount{$temp[0]}{$temp[1]}{$temp[2]}{$temp[3]} == scalar @exomiserFileList) {
			$HashOfExomiserXLS{$temp[0]}{$temp[1]}{$temp[2]}{$temp[3]} = $exomiserRecord;
		}
	}
}


# To find extended compound hetero: multiple variant hits and CNVs on the same gene 
my %all_gene_id; # to record gene_id => gene_name
my %CNV_on_samples; # to record all CNVs in a new manner [Gene][sample]=cnv_id
my %Variants_on_samples; # to record all rare variants in a new manner [Gene][sample]=chr_v-position_genoType_alleleFreq(ESP)_alleleFreq(1000G)_alleleFreq(in house);

if ($CNV_file ne "") {
	print "CNV file provided:$CNV_file\n";
	open CNV, $CNV_file or die "Cannot open the CNV file $CNV_file\n";
	my  %CNV_ids;
	while (my $line=<CNV>) {
		chomp $line;
	 	if ($line eq "") {
			next;
		}
		# remove \t within double quotes
		my @elements=split("@", $line);
		
		if($elements[0] =~ m/SampleID/) {
			next;
		}
		#print scalar @elements."\n";
		my $id = $elements[0];
		if (ElementInArray(\@sample_names, $id) eq "YES") { # if the ID is in sample names
			my $temp_temp_string="";
			foreach my $cell (@elements) {
				#print $cell."\n";
				$temp_temp_string = $temp_temp_string."\t".$cell;
			}
			$temp_temp_string =~ s/^\t//;
			push @output_CNV, $temp_temp_string;
			my $start=$elements[5];
			my $end=$elements[6];
			$elements[7] =~ s/"//g;
			my $chr="chr".$elements[7];
			$CNV_regions{$id}{$chr}{$start}{$end} = $elements[8] ;
			unless( exists $CNV_ids{$id} ) {
				$CNV_ids{$id} =0;
				print "The file contains info for $id.\n";
			}
		}
	}
}
close CNV;

open INPUT, $input_file or die "Cannot open the file $input_file.\n";
while (my $line=<INPUT>) {
	if ($line !~ m/^Chr\t/) {
		chomp $line;
		my @elements = split("\t", $line);
		my $length = scalar @elements -1;
		$elements[$length-2] =~ s/\'/\"/g;
		$elements[$length-1] =~ s/\'/\"/g;
		$elements[$length] =~ s/\'/\"/g;
		
		###### for add pnorm test in a extra column
		my $format = $elements[46];
		my @format_fields = split(":", $format);
		my $GQ_index;
		my $AD_index;
		my $GT_index;
		for (my $i=0;$i<scalar(@format_fields);$i++){
			if ($format_fields[$i] =~ m/GQ/) {
				$GQ_index = $i;
				#print $line."$GQ_index\n";
			}
			if ($format_fields[$i] =~ m/GT/) {
				$GT_index = $i;
				#print $line."$GT_index\n";
			}
			if ($format_fields[$i] =~ m/AD/) {
				$AD_index = $i;
				#print $line."$AD_index\n";
			}
		}
		####
		
		#if ($elements[26] ne "unknown") { # insert "unknown" to make every row in the same length
		#	splice @elements, 26, 0, 'unknown'; 
		#}
		my $sample_string="";
		my $rr_na_sum=0;
		for (my $i=0;$i<scalar @sample_names;$i++) {
			my ($on_CNV, $geno_flag);
			
			if ($CNV_file ne ""){ #when CNV file is available
				$on_CNV="NO";
				start_loop: foreach my $start ( sort {$a<=>$b} keys %{$CNV_regions{$sample_names[$i]}{$elements[38]}} ) {
					foreach my $end (keys %{$CNV_regions{$sample_names[$i]}{$elements[38]}{$start}}) {
						if ( $elements[39] >= $start && $elements[39] <= $end) {
							$on_CNV=$CNV_regions{$sample_names[$i]}{$elements[38]}{$start}{$end};
							last start_loop;
						}
					}
				}
			} 
			
			###### for add pnorm test in a extra column
			if ($add_genotypeCall_flags =~ m/^y/i) { #Anything less 10 times coverge, we consider it is a bad genotype call
				my ($p_bino, $GQ, $GT, $AD);
				$geno_flag = "NA";
				my $sample = $elements[$sample_columns[$i]];
				if (defined $GQ_index && defined $AD_index && defined $GT_index) {
					if ($sample !~ m/\.\/\./) {
						my @fields = split(":", $sample);
						$GT = $fields[$GT_index];
						$GQ = $fields[$GQ_index];
						$AD = $fields[$AD_index];
						#print $line."\t$GT\t$GQ\t$AD\n";
						if ($GT =~ m/0\/1|1\/1/ && $AD ne ".") {
							my @cov = split(",", $AD);
							if ($cov[0]==0 && $cov[1]==0) {
								$geno_flag = "E";
								#print "#WARNING: Call made on 0 depth is found, possible bug of GATK VQSR. the line is :\n$line\n";
							} else {
								my $test_v = $cov[1];
								$test_v = $cov[0] if $cov[0] < $cov[1];
								$p_bino = 2*&Math::CDF::pbinom($test_v, $cov[0]+$cov[1], 0.5);
								if ($GT =~ m/0\/1/) {
									if ($p_bino < 0.05) {
										$geno_flag = "B"; #bad call for hetero
									} else {
										if ($cov[0]+$cov[1] >= 10) {
											$geno_flag = "G";
										} else {
											$geno_flag = "B";
										}
									}
								} else {
									if ($p_bino <= 0.05) {
										if ($cov[0]+$cov[1] >= 10) {
											$geno_flag = "G"; #bad call for hetero
										}else {
											$geno_flag = "B";
										}
									} else {
										$geno_flag = "B";
									}
								}
							}
						}
					}
				}
			}
			####
			if ($CNV_file ne "" && $add_genotypeCall_flags =~ m/^y/i) {
				$sample_string=$sample_string.$elements[$sample_columns[$i]]."\t".$elements[$sample_columns[$i]+1]."\t".$on_CNV."\t".$geno_flag."\t";	
			} elsif ($CNV_file ne ""){
				$sample_string=$sample_string.$elements[$sample_columns[$i]]."\t".$elements[$sample_columns[$i]+1]."\t".$on_CNV."\t";
			} elsif ($add_genotypeCall_flags =~ m/^y/i) {
				$sample_string=$sample_string.$elements[$sample_columns[$i]]."\t".$elements[$sample_columns[$i]+1]."\t".$geno_flag."\t";
			} else {
				$sample_string=$sample_string.$elements[$sample_columns[$i]]."\t".$elements[$sample_columns[$i]+1]."\t";
			}			
			
			if($elements[$sample_columns[$i]+1] eq 'R/R' || $elements[$sample_columns[$i]+1] eq 'NA') {
				$rr_na_sum+=1;
			}
			
		}
		
		#Change 1000genome, ESP6500 and cg69 var freq empty value to "0"
		if($elements[20] eq "NA" ) {#1000G
			$elements[20] = 0;
		}
		if($elements[19] eq "NA" ) { #esp6500
			$elements[19] = 0;
		}
		if($elements[21] eq "NA" ) { #cg69
			$elements[21] = 0;
		}
		
		#Change in-house MAF and GATK MAF 'NA' to '0'
		if ($elements[$length-8] eq "NA") { #inHouse maf
			$elements[$length-8] = 0;
		}
		if ($elements[$length-7] eq "NA") { #inhouse gatk maf
			$elements[$length-7] = 0;
		}
		
		
		#for Exomiser output spreadsheet
		if ($exomiserXLSfiles ne "" && exists $HashOfExomiserXLS{$elements[38]}{$elements[1]}{$elements[3]}{$elements[4]}) {
			my @ExomiserSplit = split("\t", $HashOfExomiserXLS{$elements[38]}{$elements[1]}{$elements[3]}{$elements[4]});
			my $ExomiserString = "";
			for (my $i=4;$i < scalar @ExomiserSplit;$i++){
				$ExomiserString=$ExomiserString."\t".$ExomiserSplit[$i];
			}
			if ($elements[19] <= 0.05 && $elements[20] <= 0.05) { #Remove high freq vars in ESP6500 and 1000G
				push @ArrayOfExomiserXLS_Output, $elements[38]."\t".$elements[39]."\t".$elements[41]."\t".$elements[42].$ExomiserString."\t".
				$elements[43]."\t".$elements[5]."\t".$elements[7]."\t".$elements[8]."\t".$elements[9]
				."\t".$elements[11]."\t".$elements[12]."\t".$elements[13]."\t".$elements[15]."\t".$elements[16]."\t".$elements[$length-11]."\t"
				.$elements[$length-10]."\t".$elements[$length-9]."\t".$elements[18]."\t".$elements[19]."\t".$elements[20]."\t".$elements[21]."\t".
				$elements[$length-8]."\t".$elements[$length-7]."\t".$elements[22]."\t".$elements[25].":".$elements[27].":".$elements[29].":".$elements[31]."\t".$elements[46]."\t".$sample_string.$elements[$length-2]."\t".
				$elements[$length-1]."\t".$elements[$length]."\t".$elements[17]."\t".$elements[23]."\t".$elements[24]."\t".$elements[25]."\t".$elements[26]."\t".$elements[27]."\t"
				.$elements[28]."\t".$elements[29]."\t".$elements[30]."\t".$elements[31]."\t".$elements[32]."\t".$elements[33]."\t".$elements[34]."\t".$elements[35]."\t"
				.$elements[36]."\t".$elements[37];
			}
			
		}
		
#my @header = ('Chr','Start','Ref','Obs','Func.ensGene','ExonicFunc.ensGene','AAChange.ensGene','GeneName','Predictions')

#my @header_2 = ('InHouseMAF_GATK','OMIM Link','GeneCard Link','Uniprot Link','OMIM_Gene_Description','WikiGene_Description','GoTerm', 'isInterested',
#'Format',);

#my $superDupsRow=19;

#my @header_3 = ('Func.knownGene','ExonicFunc.knownGene','AAChange.knownGene','Func.refGene','ExonicFunc.refGene','AAChange.refGene',
#'ESP6500_ALL','1000g2012feb_ALL','dbSNP137','VariantCall_quality', 'GeneID', 'genomicSuperDups', 'SNPorINDEL', 'cg69', 'InHouseMAF', 
#'phastConsElements46way','LJB2_SIFT','LJB2_PolyPhen2_HDIV','LJB2_PP2_HDIV_Pred','LJB2_PolyPhen2_HVAR','LJB2_PolyPhen2_HVAR_Pred',
#'LJB2_LRT','LJB2_LRT_Pred','LJB2_MutationTaster','LJB2_MutationTaster_Pred','LJB_MutationAssessor','LJB_MutationAssessor_Pred',
#'LJB2_FATHMM','LJB2_GERP++','LJB2_PhyloP','LJB2_SiPhy');

		#re-arrange columns
		my $string_1 = $elements[38]."\t".$elements[39]."\t".$elements[41]."\t".$elements[42]."\t".$elements[5]."\t".$elements[7]."\t".$elements[8]."\t"
						.$elements[$length-11]."\t".$elements[25].":".$elements[27].":".$elements[29].":".$elements[31]."\t";
		
		my $string_2 =	$elements[$length-7]."\t".$elements[$length-1]."\t".$elements[$length-2]."\t".$elements[$length]."\t".$elements[$length-3]."\t"
						.$elements[$length-5]."\t".$elements[$length-6]."\t".$elements[$length-9]."\t".$elements[46]."\t";	
		
		my $string_3 = $elements[9]."\t".$elements[11]."\t".$elements[12]."\t".$elements[13]."\t".$elements[15]."\t".$elements[16]."\t".$elements[19]."\t"
						.$elements[20]."\t".$elements[22]."\t".$elements[43]."\t".$elements[$length-10]."\t".$elements[18]."\t".$elements[$length-12]."\t"
						.$elements[21]."\t".$elements[$length-8]."\t".$elements[17]."\t".$elements[23]."\t".$elements[24]."\t".$elements[25]."\t"
						.$elements[26]."\t".$elements[27]."\t".$elements[28]."\t".$elements[29]."\t".$elements[30]."\t".$elements[31]."\t".$elements[32]."\t"
						.$elements[33]."\t".$elements[34]."\t".$elements[35]."\t".$elements[36]."\t".$elements[37]; #re-arrange columns
		
		my $sample_string_final;
		my @string_temp = split("\t", $sample_string);
		my $sample_string_temp="";
		for (my $i=0;$i<scalar @sample_names;$i++) {
			if ($CNV_file ne "" && $add_genotypeCall_flags =~ m/^y/i) {
				my $j = $i*4;
				$sample_string_temp=$sample_string_temp.$string_temp[$j+1]."\t".$string_temp[$j+3]."\t";
			} elsif ($CNV_file ne ""){
				my $j = $i*3;
				$sample_string_temp=$sample_string_temp.$string_temp[$j+1]."\t";
			} elsif ($add_genotypeCall_flags =~ m/^y/i) {
				my $j = $i*3;
				$sample_string_temp=$sample_string_temp.$string_temp[$j+1]."\t".$string_temp[$j+2]."\t";
			} else {
				my $j = $i*2;
				$sample_string_temp=$sample_string_temp.$string_temp[$j+1]."\t";
			}
		}
		$sample_string_final=$string_1.$sample_string_temp.$string_2;
		$sample_string_temp="";
		for (my $i=0;$i<scalar @sample_names;$i++) {
			if ($CNV_file ne "" && $add_genotypeCall_flags =~ m/^y/i) {
				my $j = $i*4;
				$sample_string_temp=$sample_string_temp.$string_temp[$j]."\t".$string_temp[$j+2]."\t";
			} elsif ($CNV_file ne ""){
				my $j = $i*3;
				$sample_string_temp=$sample_string_temp.$string_temp[$j]."\t".$string_temp[$j+2]."\t";
			} elsif ($add_genotypeCall_flags =~ m/^y/i) {
				my $j = $i*3;
				$sample_string_temp=$sample_string_temp.$string_temp[$j]."\t";
			} else {
				my $j = $i*2;
				$sample_string_temp=$sample_string_temp.$string_temp[$j]."\t";
			}
		}
		$sample_string_final=$sample_string_final.$sample_string_temp.$string_3;
		
			
		if($rr_na_sum < scalar @sample_names) { # RR and NA are removed
			push @output_everything, $sample_string_final;
		}
		
		if($elements[7] eq "synonymous SNV" && $elements[11] eq "synonymous SNV" && $elements[15] eq "synonymous SNV"  ) { # eliminate synonymous SNV
			next;
		}
		if($elements[5] =~ m/(intronic|intergenic|upstream|downstream)/i &&  $elements[9] =~ m/(intronic|intergenic|upstream|downstream)/i 
						&& $elements[13] =~ m/(intronic|intergenic|upstream|downstream)/i  ) { # eliminate intronic or intergenic variants
			next;
		}
		
		if($rr_na_sum < scalar @sample_names) { # RR and NA are removed
			push @output_all, $sample_string_final;
			
			# for filtered spreadsheet

			if( $elements[19] <= 0.05 && $elements[20] <= 0.05 && $elements[21] <= 0.05 && 
			$elements[$length-8] <= 0.05 && $elements[$length-7] <= 0.05) {
				push @output_filtered, $sample_string_final;
				
				#unless(exists $gene_name_count{$elements[$length-9]}) { # could be deleted
				#	$gene_name_count{$elements[$length-9]} = 0;
				#}
				
				# for deleterious hits spreadsheet
				my $number_deleterious_predictions=0;
				if($elements[25] eq "D" || $elements[25] eq "P") {
					$number_deleterious_predictions+=1;
				}
				if($elements[27] eq "D" || $elements[27] eq "P") {
					$number_deleterious_predictions+=1;
				}
				if($elements[29] eq "D") {
					$number_deleterious_predictions+=1;
				}
				if($elements[31] eq "A" || $elements[31] eq "D") {
					$number_deleterious_predictions+=1;
				}
				
				my $number_CompoundHet=0;
				foreach my $sample_col (@sample_columns) {
					if($elements[$sample_col+1] eq "V1/V2") {
						$number_CompoundHet +=1;
					}
				}
				
				for( my $i=0; $i< scalar @sample_columns;$i++) {
					my $sample_col = $sample_columns[$i];
					#if($elements[$sample_col+1] eq "V1/V2") {
					#	$number_CompoundHet +=1;
					#}
					
					# extend the hash to record variants for samples and genes for %Variants_on_samples: $variants_on_samples{sample_name}{gene_id}=chr_pos_genoType
					unless (exists $all_gene_id{$elements[$length-10]}) { 
						$all_gene_id{$elements[$length-10]} = $elements[$length-11]."\t".$elements[$length-9]."\t".
						$elements[$length-6]."\t".$elements[$length-5]."\t".$elements[$length-4]."\t".$elements[$length-3]."\t".
						$elements[$length-2]."\t".$elements[$length-1]."\t".$elements[$length]; 
					}
					if ( exists $Variants_on_samples{$sample_names[$i]}{$elements[$length-10]}  && $elements[$sample_col+1] ne "R/R" && $elements[$sample_col+1] ne "NA" ) {
						$Variants_on_samples{$sample_names[$i]}{$elements[$length-10]} = 
						$Variants_on_samples{$sample_names[$i]}{$elements[$length-10]}.";".
						$elements[38]."_".$elements[39]."_".$elements[$sample_col+1]."_".sprintf( "%.4f",$elements[19])."_".sprintf( "%.4f",$elements[20])."_".sprintf( "%.4f",$elements[$length-8]);
					} elsif( $elements[$sample_col+1] ne "R/R" && $elements[$sample_col+1] ne "NA" ) {
						$Variants_on_samples{$sample_names[$i]}{$elements[$length-10]} = 
						$elements[38]."_".$elements[39]."_".$elements[$sample_col+1]."_".sprintf( "%.4f",$elements[19])."_".sprintf( "%.4f",$elements[20])."_".sprintf( "%.4f",$elements[$length-8]);				
					}
				}
				
				#for rare variants on a same gene
				#if($number_CompoundHet < scalar @sample_names) {
				#	$gene_name_count{$elements[$length-11]} =$gene_name_count{$elements[$length-11]}+1;
				#} #retired
				
				if ($elements[$length-9] eq "YES" && $number_deleterious_predictions >=1) {  #those that were predicted as “Deleterious” by at least two predictors
					push @output_possibleHits, $sample_string_final;
				} elsif ($elements[$length-9] eq "NO" && $number_deleterious_predictions >=3) { #those that were predicted as “Deleterious” by all 4 predictors
					push @output_possibleHits, $sample_string_final;
				} elsif($elements[$length-12] eq "SNP" && $number_CompoundHet >=1) { #Compound heterozygous “V1/V2”* of single base polymorphisms
					push @output_possibleHits, $sample_string_final;
				} elsif($elements[$length-9] eq "YES" && $elements[19] <= 0.01 && $elements[20] <= 0.01) { # those that have AAF less than 0.01 in both 1000Genomes project and ESP6500 
					push @output_possibleHits, $sample_string_final;
				}
				
				if ($elements[$length-9] eq "YES" && $elements[19] <= 0.01 && $elements[20] <= 0.01) {
					push @very_rare_interested_vars, $sample_string_final;
				}
				
				#for Xlinked rare variants
				if ($elements[38] =~ m/chrx/i) {
					push @output_Xlinked, $sample_string_final;
				}
				
				# AR model Hits
				my $number_Homo=0;
				foreach my $sample_col (@sample_columns) {
					if($elements[$sample_col+1] eq "V/V") {
						$number_Homo +=1;
					}
				}
				if ($number_Homo == scalar @sample_names) {
					if ( $elements[38] !~ m/chr[XYMxym]/ ) {
						push @output_ARmodelHits, $sample_string_final;
					}
				}
			}
		}
	}
}
close(INPUT);


# record CNVs into the hash %CNV_on_samples
if ($CNV_file ne "") {
	####
	# need to read in the whole list of insterested genes in ENSEMBL id
	####
	foreach my $sample (keys %CNV_regions) {
		foreach my $chr (keys %{$CNV_regions{$sample}}) {
			foreach my $start_CNV (keys %{$CNV_regions{$sample}{$chr}}) {
				foreach my $end_CNV (keys %{$CNV_regions{$sample}{$chr}{$start_CNV}}) {
					start_loop: foreach my $start_gene ( sort {$a<=>$b} keys %{$en_genes{$chr}} ) {
						foreach my $end_gene (keys %{$en_genes{$chr}{$start_gene}}) {
							if ( $start_gene <= $end_CNV && $end_gene >= $start_CNV ) {
								my @temp = split("\t", $en_genes{$chr}{$start_gene}{$end_gene});
								my $gene_name = $temp[0];
								my $gene_id = $temp[1];
								unless (exists $all_gene_id{$gene_id}) {
									### to add more annotation to the genes as requested
									my $temp_string;
									if ( ! exists $isInterestedGenes{$gene_id}) {
										$temp_string = $gene_name."\t"."NO"; 
									} else {
										$temp_string = $gene_name."\t"."YES"; 
									}
									my $ens_id = $gene_id;
									my $omim_anno = "NA\tNA\tNA\tNA";
									if (exists $OMIM{$ens_id}) {
										$omim_anno = $OMIM{$ens_id};
									}
									my $genecard_link = "=HYPERLINK(\'http://www.genecards.org/cgi-bin/carddisp.pl?id=$ens_id&id_type=ensembl\', \'GeneCard Link\')";
									my $omim_link = "";
									if (exists $OmimAcc{$ens_id}) {
										my $omim_accesion =  $OmimAcc{$ens_id};
										$omim_link = "=HYPERLINK(\'http://omim.org/entry/$omim_accesion\', \'OMIM Link\')";
									} else {
										$omim_link = "NA";
									}
									
									my $uniprot_link = "";
									if (exists $UniprotAcc{$ens_id}) {
										my $uniprot_accession = $UniprotAcc{$ens_id};
										$uniprot_link = "=HYPERLINK(\'http://www.uniprot.org/uniprot/$uniprot_accession\', \'Uniprot Link\')";
									} else {
										$uniprot_link = "NA";
									}
									$genecard_link =~ s/\'/\"/g;
									$omim_link =~ s/\'/\"/g;
									$uniprot_link =~ s/\'/\"/g;
									$all_gene_id{$gene_id} = $temp_string."\t".$omim_anno."\t".$genecard_link."\t".$omim_link."\t".$uniprot_link;
					
								}
								if(exists $CNV_on_samples{$sample}{$gene_id}) {
									$CNV_on_samples{$sample}{$gene_id}=$CNV_on_samples{$sample}{$gene_id}.";".$CNV_regions{$sample}{$chr}{$start_CNV}{$end_CNV};
								} else {
									$CNV_on_samples{$sample}{$gene_id} = $CNV_regions{$sample}{$chr}{$start_CNV}{$end_CNV};
								}
							} elsif ( $end_gene < $start_CNV) {
								next start_loop;
							} elsif ($start_gene > $end_CNV) {
								last start_loop;
							}
						}
					}
				}
			}
		}
	}
}


###
###########
####################
################################################ output in excel format

print "Creating excel output...\n";

my @header = ('Chr','Start','Ref','Obs','Func.ensGene','ExonicFunc.ensGene','AAChange.ensGene','GeneName','Predictions');

my @header_2 = ('InHouseMAF_GATK','OMIM Link','GeneCard Link','Uniprot Link','OMIM_Gene_Description','WikiGene_Description','GoTerm', 'isInterested',
'Format');

my @header_3 = ('Func.knownGene','ExonicFunc.knownGene','AAChange.knownGene','Func.refGene','ExonicFunc.refGene','AAChange.refGene',
'ESP6500_ALL','1000g2012feb_ALL','dbSNP137','VariantCall_quality', 'GeneID', 'genomicSuperDups', 'SNPorINDEL', 'cg69', 'InHouseMAF', 
'phastConsElements46way','LJB2_SIFT','LJB2_PolyPhen2_HDIV','LJB2_PP2_HDIV_Pred','LJB2_PolyPhen2_HVAR','LJB2_PolyPhen2_HVAR_Pred',
'LJB2_LRT','LJB2_LRT_Pred','LJB2_MutationTaster','LJB2_MutationTaster_Pred','LJB_MutationAssessor','LJB_MutationAssessor_Pred',
'LJB2_FATHMM','LJB2_GERP++','LJB2_PhyloP','LJB2_SiPhy');

my $superDupsRow=8;

foreach my $sample (@sample_names) {	
	push @header, $sample.".anno";
	$superDupsRow += 1;
	if ($add_genotypeCall_flags =~ m/^y/i) {
		push @header, $sample.".flag";
		$superDupsRow += 1;
	}
}
push @header, @header_2;
foreach my $sample (@sample_names) {
	push @header, $sample;
	$superDupsRow += 1;	
	if ($CNV_file ne ""){
		push @header, $sample.".cnv";
		$superDupsRow += 1;
	}
	
}
push @header, @header_3;
$superDupsRow = $superDupsRow + 21;


my $workbook = Excel::Writer::XLSX->new($output_excel);

#format settings
my $format_header = $workbook->add_format();
$format_header->set_bold();
my $format_segmental_dup = $workbook->add_format(color => 10);

#my $format_possibleHits_Y = $workbook->add_format();
#$format_possibleHits_Y->set_color('lime');
#$format_possibleHits_Y->set_bg_color('green');
#my $format_possibleHits_N = $workbook->add_format();
#$format_possibleHits_N->set_color('red');
#my $format_possibleHits_C = $workbook->add_format();
#$format_possibleHits_C->set_bg_color('yellow');

if ($output_excel_everything ne "") {
	my $workbook_everything = Excel::Writer::XLSX->new($output_excel_everything);
	my $format_header_everything = $workbook_everything->add_format();
	$format_header_everything->set_bold();
	my $format_segmental_dup_everything = $workbook_everything->add_format(color => 10);
	
	my $worksheet_everything = $workbook_everything->add_worksheet('ALL');
	$worksheet_everything->write_row(0,0,\@header, $format_header_everything);
	for (my $i=0;$i< scalar @output_everything;$i++) {
		my @row = split("\t", $output_everything[$i]);
		#print "$row[$superDupsRow]\n";
		if ($row[$superDupsRow] ne "NA") {
			$worksheet_everything->write_row($i+1,0,\@row, $format_segmental_dup_everything);
		} else {
			$worksheet_everything->write_row($i+1,0,\@row);
		}
	}
	$workbook_everything->close();
} else {
	my $worksheet_everything = $workbook->add_worksheet('ALL');
	$worksheet_everything->write_row(0,0,\@header, $format_header);
	for (my $i=0;$i< scalar @output_everything;$i++) {
		my @row = split("\t", $output_everything[$i]);
		if ($row[$superDupsRow] ne "NA") {
			$worksheet_everything->write_row($i+1,0,\@row, $format_segmental_dup);
		} else {
			$worksheet_everything->write_row($i+1,0,\@row);
		}
	}
}

my $worksheet_all = $workbook->add_worksheet('ExonicVariants');
$worksheet_all->write_row(0,0,\@header, $format_header);
for (my $i=0;$i< scalar @output_all;$i++) {
	my @row = split("\t", $output_all[$i]);
	if ($row[$superDupsRow] ne "NA") {
		$worksheet_all->write_row($i+1,0,\@row, $format_segmental_dup);
	} else {
		$worksheet_all->write_row($i+1,0,\@row);
	}
}

my $worksheet_filtered = $workbook->add_worksheet('RareExonicVariants');
$worksheet_filtered->write_row(0,0,\@header, $format_header);
for (my $i=0;$i< scalar @output_filtered;$i++) {
	my @row = split("\t", $output_filtered[$i]);
	if ($row[$superDupsRow] ne "NA") {
		$worksheet_filtered->write_row($i+1,0,\@row, $format_segmental_dup);
	} else {
		$worksheet_filtered->write_row($i+1,0,\@row);
	}
}

my $worksheet_xlinked = $workbook->add_worksheet('XLinked');
$worksheet_xlinked->write_row(0,0,\@header, $format_header);
for (my $i=0;$i< scalar @output_Xlinked;$i++) {
	my @row = split("\t", $output_Xlinked[$i]);
	if ($row[$superDupsRow] ne "NA") {
		$worksheet_xlinked->write_row($i+1,0,\@row, $format_segmental_dup);
	} else {
		$worksheet_xlinked->write_row($i+1,0,\@row);
	}
}

my $worksheet_possibleHits = $workbook->add_worksheet('Deleterious_Hits');
$worksheet_possibleHits ->write_row(0,0,\@header, $format_header);
for (my $i=0;$i< scalar @output_possibleHits;$i++) {
	my @row = split("\t", $output_possibleHits[$i]);
	if ($row[$superDupsRow] ne "NA") {
		$worksheet_possibleHits->write_row($i+1,0,\@row, $format_segmental_dup);
	} else {
		$worksheet_possibleHits->write_row($i+1,0,\@row);
	}
	
}

my $worksheet_rareInterestedHits = $workbook->add_worksheet('Rare_Interested_Hits');
$worksheet_rareInterestedHits ->write_row(0,0,\@header, $format_header);
for (my $i=0;$i< scalar @very_rare_interested_vars;$i++) {
	my @row = split("\t", $very_rare_interested_vars[$i]);
	if ($row[$superDupsRow] ne "NA") {
		$worksheet_rareInterestedHits->write_row($i+1,0,\@row, $format_segmental_dup);
	} else {
		$worksheet_rareInterestedHits->write_row($i+1,0,\@row);
	}
}

my $worksheet_ARmodelHits = $workbook->add_worksheet('Homozygous_Hits');
$worksheet_ARmodelHits ->write_row(0,0,\@header, $format_header);
for (my $i=0;$i< scalar @output_ARmodelHits;$i++) {
	my @row = split("\t", $output_ARmodelHits[$i]);
	if ($row[$superDupsRow] ne "NA") {
		$worksheet_ARmodelHits->write_row($i+1,0,\@row, $format_segmental_dup);
	} else {
		$worksheet_ARmodelHits->write_row($i+1,0,\@row);
	}
}


if ($CNV_file ne "") {
	my @CNV_header=('SampleID', 'start.p', 'end.p', 'type', 'nexons', 'start', 'end', 'chromosome', 'id', 'BF', 'reads.expected', 'reads.observed', 	
	'reads.ratio', 'Conrad.hg19', 'exons.hg19', 'Averagecontrols_DupDel_acrossCNV', 	
	'min_controls', 'max_controls', 'FullGenesNames', 'GO-terms', 'OMIM');
	my $worksheet_CNV = $workbook->add_worksheet('CNV');
	$worksheet_CNV -> write_row(0,0,\@CNV_header, $format_header);
	for (my $i=0;$i< scalar @output_CNV;$i++) {
		my @row = split("\t", $output_CNV[$i]);
		$worksheet_CNV -> write_row($i+1,0,\@row);
	}
}

if ($exomiserXLSfiles ne "") {
	my @exomiser_header=('Chr','Start','Ref','Obs','Exomiser_gene','Score_no_prior','Score_mouse_exist','nn_change','Pathpogenicity','Path_score','Mouse_pheno','OMIM_anno'	
	,'VariantCall_quality','Func.ensGene','ExonicFunc.ensGene','AAChange.ensGene','Func.knownGene','ExonicFunc.knownGene',
	'AAChange.knownGene','Func.refGene','ExonicFunc.refGene','AAChange.refGene','GeneName','GeneID','isInterested','genomicSuperDups',
	'ESP6500_ALL','1000g2012feb_ALL','cg69','InHouseMAF','InHouseMAF_GATK','dbSNP137','Predictions', 'Format');
	foreach my $sample (@sample_names) {
		push @exomiser_header, $sample;
		push @exomiser_header, $sample.".anno";
		if ($CNV_file ne "") {
			push @exomiser_header, $sample.".cnv";
		}
	}
	my @hula_header = ('GeneCard Link', 'OMIM Link', 'Uniprot Link','phastConsElements46way','LJB2_SIFT','LJB2_PolyPhen2_HDIV','LJB2_PP2_HDIV_Pred','LJB2_PolyPhen2_HVAR',
	'LJB2_PolyPhen2_HVAR_Pred','LJB2_LRT','LJB2_LRT_Pred','LJB2_MutationTaster','LJB2_MutationTaster_Pred','LJB_MutationAssessor','LJB_MutationAssessor_Pred',
	'LJB2_FATHMM','LJB2_GERP++','LJB2_PhyloP','LJB2_SiPhy');
	push @exomiser_header, @hula_header;
	my $worksheet_exomiser = $workbook->add_worksheet('Exomiser');
	$worksheet_exomiser -> write_row(0,0,\@exomiser_header, $format_header);
	my @Sorted_ArrayOfExomiserXLS_Output = map  { $_->[0] }
										sort { $b->[1]<=>$a->[1] }
										map {my @temp=split("\t",$_); [$_, $temp[5]]} @ArrayOfExomiserXLS_Output; #Schwartzian Transform
	for (my $i=0;$i< scalar @Sorted_ArrayOfExomiserXLS_Output;$i++) {
		my @row = split("\t", $Sorted_ArrayOfExomiserXLS_Output[$i]);
		if ($row[25] ne "NA") {
			$worksheet_exomiser -> write_row($i+1,0,\@row, $format_segmental_dup);
		} else {
			$worksheet_exomiser -> write_row($i+1,0,\@row);
		}			
	}
}


############################# output compound heterozygous sheet ###################################################
my @compound_header = ('Gene_ID', 'Gene_Name', 'IsInterested', 'GoTerm', 'WikiGene_Description', 'MIM_Gene_Description', 'OMIM_Gene_Description',
									 'GeneCard Link', 'OMIM Link', 'Uniprot Link');
foreach my $sample (@sample_names) {
	push @compound_header, $sample."_vars";
	push @compound_header, $sample."_CNV";
}
my $worksheet_compound = $workbook->add_worksheet('Compound_Heterozygous');
$worksheet_compound -> write_row(0,0,\@compound_header, $format_header);
my $number_of_wide_cols = (scalar @sample_names)*2+9;
$worksheet_compound -> set_column( 0, 0 , 19 );
$worksheet_compound -> set_column( 1, 1 , 13 );
$worksheet_compound -> set_column( 2, 2 , 13 );
$worksheet_compound -> set_column( 10, $number_of_wide_cols , 40 );


my %gene_sample_CNV_count; #number of how many CNVs on the gene
my %gene_sample_VAR_count; #number of how many snp/indel on the gene

foreach  my $sample (@sample_names) {
	foreach my $gene_id (keys %{$CNV_on_samples{$sample}}) {
		my @temp = split(";",$CNV_on_samples{$sample}{$gene_id});
		if(exists $gene_sample_CNV_count{$sample}{$gene_id}) {
			$gene_sample_CNV_count{$sample}{$gene_id} = $gene_sample_CNV_count{$sample}{$gene_id} + scalar @temp;
		} else {
			$gene_sample_CNV_count{$sample}{$gene_id} = scalar @temp;
		}
	}
}

foreach  my $sample (@sample_names) {
	foreach my $gene_id (keys %{$Variants_on_samples{$sample}}) {
		$gene_sample_VAR_count{$sample}{$gene_id} = 0;
		my @temp = split(";",$Variants_on_samples{$sample}{$gene_id});
		# Remove 'R/R' 'NA' and same position variants
		my %positions;
		foreach my $var_info (@temp) {
			if ($var_info !~ m/\_R\/R\_/ && $var_info !~ m/NA\_/) {
				#print "var is $var_info \n";
				my @temp_temp = split("_", $var_info);
				#print @temp_temp."\n";
				######### testing here  ##################
				#foreach my $wolegequ (@temp_temp){
				#	print $wolegequ."  ";
				#}
				#print "\n";
				##########################################
				if (!exists $positions{$temp_temp[1]}) {
					$positions{$temp_temp[1]} = 0;
					if(exists $gene_sample_VAR_count{$sample}{$gene_id}) {
						#print "+1\n";
						$gene_sample_VAR_count{$sample}{$gene_id} = $gene_sample_VAR_count{$sample}{$gene_id} + 1;
					} 
				}
			} 
		}
		#testing here
		#print $sample."_".$gene_id."\t".$Variants_on_samples{$sample}{$gene_id}." Count result ".$gene_sample_VAR_count{$sample}{$gene_id}."\n";
	}
}

# ready to output to excel sheet
my $excel_line_index=1;
my $compond_format= $workbook->add_format(text_wrap => 1);
foreach my $gene_id (keys %all_gene_id) {
	my @output;
	my $output_light = "red";
	push  @output, $gene_id;
	my @temp_yuanyuan = split("\t", $all_gene_id{$gene_id});
	foreach my $temp_yuanyuan_string (@temp_yuanyuan) {
	  push  @output, $temp_yuanyuan_string;
	}
	my $row_height =1;
	foreach my $sample (@sample_names) {
		if( exists $CNV_on_samples{$sample}{$gene_id} && exists $Variants_on_samples{$sample}{$gene_id}) {
			my @gene_sample_var = split(";", $Variants_on_samples{$sample}{$gene_id});
			my @gene_sample_cnv = split(";", $CNV_on_samples{$sample}{$gene_id});
			my $temp = "";
			for(my $i=0;$i<scalar @gene_sample_var;$i++) {
				if ($i == scalar @gene_sample_var -1 ) {
					$temp = $temp.$gene_sample_var[$i];
				} else {
					$temp = $temp.$gene_sample_var[$i]."\n";
				}
			}
			################################# need to filter the $temp before push in to array.
			push @output, $temp;
			$temp = "";
			for(my $i=0;$i<scalar @gene_sample_cnv;$i++) {
				if ($i == scalar @gene_sample_cnv -1 ) {
					$temp = $temp.$gene_sample_cnv[$i];
				} else {
					$temp = $temp.$gene_sample_cnv[$i]."\n";
				}
			}
			push @output, $temp;
			if ($gene_sample_CNV_count{$sample}{$gene_id} + $gene_sample_VAR_count{$sample}{$gene_id} >= 2) {
				$output_light = "green";
			}
			
			# decided the largest row_height for excel cell
			if (scalar @gene_sample_cnv >= scalar @gene_sample_var) {
				if (scalar @gene_sample_cnv > $row_height) {
					$row_height = scalar @gene_sample_cnv;
				}
			} else {
				if (scalar @gene_sample_var > $row_height) {
					$row_height = scalar @gene_sample_var;
				}
			}
			
		} elsif (exists $CNV_on_samples{$sample}{$gene_id} && !exists $Variants_on_samples{$sample}{$gene_id}) {
			my @gene_sample_cnv = split(";", $CNV_on_samples{$sample}{$gene_id});
			push @output, 'NA';
			my $temp = "";
			for(my $i=0;$i<scalar @gene_sample_cnv;$i++) {
				if ($i == scalar @gene_sample_cnv -1 ) {
					$temp = $temp.$gene_sample_cnv[$i];
				} else {
					$temp = $temp.$gene_sample_cnv[$i]."\n";
				}
			}
			push @output, $temp;
			if ($gene_sample_CNV_count{$sample}{$gene_id} >= 2) {
				$output_light = "green";
			}
			if (scalar @gene_sample_cnv > $row_height) {
				$row_height = scalar @gene_sample_cnv;
			}
		} elsif(!exists $CNV_on_samples{$sample}{$gene_id} && exists $Variants_on_samples{$sample}{$gene_id}) {
			my @gene_sample_var = split(";", $Variants_on_samples{$sample}{$gene_id});
			my $temp = "";
			for(my $i=0;$i<scalar @gene_sample_var;$i++) {
				if ($i == scalar @gene_sample_var -1 ) {
					$temp = $temp.$gene_sample_var[$i];
				} else {
					$temp = $temp.$gene_sample_var[$i]."\n";
				}
			}
			push @output, $temp;
			push @output, 'NA';
			if ($gene_sample_VAR_count{$sample}{$gene_id} >= 2) {
				$output_light = "green";
			}
			
			if (scalar @gene_sample_var > $row_height) {
				$row_height = scalar @gene_sample_var;
			}
		} elsif(!exists $CNV_on_samples{$sample}{$gene_id} && !exists $Variants_on_samples{$sample}{$gene_id}) {
			push @output, 'NA';
			push @output, 'NA';
		}
	}
	
	if ($output_light eq "green") {
		######## adjust row height according to the number of lines in the cell
		my $final_row_height = $row_height*17;
		$worksheet_compound -> set_row ($excel_line_index, $final_row_height);
		$worksheet_compound -> write_row($excel_line_index,0,\@output, $compond_format);
		$excel_line_index += 1;
	}
}
############################# output compound heterozygous sheet  end here ###################################################

############################# output compound heterozygous sheet for Autosomal recessive Model ###################################################
if (scalar @sample_names >= 2) {

	my $worksheet_compound_AR = $workbook->add_worksheet('Compound_Heterozygou_AR_Model');
	$worksheet_compound_AR -> write_row(0,0,\@compound_header, $format_header);
	$worksheet_compound_AR -> set_column( 0, 0 , 19 );
	$worksheet_compound_AR -> set_column( 1, 1 , 13 );
	$worksheet_compound_AR -> set_column( 2, 2 , 13 );
	$worksheet_compound_AR -> set_column( 10, $number_of_wide_cols , 40 );
	
	# ready to output to excel sheet
	$excel_line_index=1;
	foreach my $gene_id (keys %all_gene_id) {
		my @output;
		my $output_light = "red";
		my %gene_sample_count; # to record variants for each sample
		push  @output, $gene_id;
		my @temp_yuanyuan = split("\t", $all_gene_id{$gene_id});
		foreach my $temp_yuanyuan_string (@temp_yuanyuan) {
		  push  @output, $temp_yuanyuan_string;
		}
		my $row_height =1;
		foreach my $sample (@sample_names) {
			if( exists $CNV_on_samples{$sample}{$gene_id} && exists $Variants_on_samples{$sample}{$gene_id}) {
				my @gene_sample_var = split(";", $Variants_on_samples{$sample}{$gene_id});
				my @gene_sample_cnv = split(";", $CNV_on_samples{$sample}{$gene_id});
				my $temp = "";
				for(my $i=0;$i<scalar @gene_sample_var;$i++) {
					if ($i == scalar @gene_sample_var -1 ) {
						$temp = $temp.$gene_sample_var[$i];
					} else {
						$temp = $temp.$gene_sample_var[$i]."\n";
					}
				}
				push @output, $temp;
				$temp = "";
				for(my $i=0;$i<scalar @gene_sample_cnv;$i++) {
					if ($i == scalar @gene_sample_cnv -1 ) {
						$temp = $temp.$gene_sample_cnv[$i];
					} else {
						$temp = $temp.$gene_sample_cnv[$i]."\n";
					}
				}
				push @output, $temp;
				
				$gene_sample_count{$sample} = $gene_sample_CNV_count{$sample}{$gene_id} + $gene_sample_VAR_count{$sample}{$gene_id};
				
				# decided the largest row_height for excel cell
				if (scalar @gene_sample_cnv >= scalar @gene_sample_var) {
					if (scalar @gene_sample_cnv > $row_height) {
						$row_height = scalar @gene_sample_cnv;
					}
				} else {
					if (scalar @gene_sample_var > $row_height) {
						$row_height = scalar @gene_sample_var;
					}
				}
				
			} elsif (exists $CNV_on_samples{$sample}{$gene_id} && !exists $Variants_on_samples{$sample}{$gene_id}) {
				my @gene_sample_cnv = split(";", $CNV_on_samples{$sample}{$gene_id});
				push @output, 'NA';
				my $temp = "";
				for(my $i=0;$i<scalar @gene_sample_cnv;$i++) {
					if ($i == scalar @gene_sample_cnv -1 ) {
						$temp = $temp.$gene_sample_cnv[$i];
					} else {
						$temp = $temp.$gene_sample_cnv[$i]."\n";
					}
				}
				push @output, $temp;
				$gene_sample_count{$sample} = $gene_sample_CNV_count{$sample}{$gene_id};
				
				if (scalar @gene_sample_cnv > $row_height) {
					$row_height = scalar @gene_sample_cnv;
				}
			} elsif(!exists $CNV_on_samples{$sample}{$gene_id} && exists $Variants_on_samples{$sample}{$gene_id}) {
				my @gene_sample_var = split(";", $Variants_on_samples{$sample}{$gene_id});
				my $temp = "";
				for(my $i=0;$i<scalar @gene_sample_var;$i++) {
					if ($i == scalar @gene_sample_var -1 ) {
						$temp = $temp.$gene_sample_var[$i];
					} else {
						$temp = $temp.$gene_sample_var[$i]."\n";
					}
				}
				push @output, $temp;
				push @output, 'NA';
				$gene_sample_count{$sample} = $gene_sample_VAR_count{$sample}{$gene_id};
				
				if (scalar @gene_sample_var > $row_height) {
					$row_height = scalar @gene_sample_var;
				}
			} elsif(!exists $CNV_on_samples{$sample}{$gene_id} && !exists $Variants_on_samples{$sample}{$gene_id}) {
				push @output, 'NA';
				push @output, 'NA';
				$gene_sample_count{$sample} = 0;
			}
		}
		
		# if each sample has more than 2 hits then output it
		$output_light = "green";
		#print $gene_id."\t";
		foreach my $sample_sample (keys %gene_sample_count) {
			if ($gene_sample_count{$sample_sample} < 2) {
				$output_light = "red";
			}
			#print $sample_sample."_".$gene_sample_count{$sample_sample}.$output_light."\t";
		}
		#print "\n";
		
		if ($output_light eq "green") {
			######## adjust row height according to the number of lines in the cell
			my $final_row_height = $row_height*17;
			$worksheet_compound_AR -> set_row ($excel_line_index, $final_row_height);
			$worksheet_compound_AR -> write_row($excel_line_index,0,\@output, $compond_format);
			$excel_line_index += 1;
		}
	}
}
############################# output compound heterozygous sheet for Autosomal recessive Model Ends Here ###################################################
$workbook->close();

exit;

sub ElementInArray {
	my ($array_ref, $item) = @_;
	my $result="NO";
	my @array_test = @{$array_ref};
	foreach my $element (@array_test) {
		if($element eq $item) {
			$result = "YES";
		}
	}
	return $result;
}

sub usage {
    print "Unknown option: @_\n" if ( @_ );
    print "\nusage: VCF_2_annotated_excel_20131120.pl \n";
    print "--vcf input vcf file (of a single sample or a family;\n";
    print "--InterestedGenes file of interested gene names list (optional); Format: Ensembl gene IDs on the 1st column.\n"; 
    print "--out output excel file name (optional, only when required);\n";
    print "--outAll output \'Everthing\' table to a sparate file to reduce memory loading excel file (optional, only when required).\n";
    print "--CNV CNV result file, output of HG's Annotate_CNVs_combine_multiple_files.pl, sample names must be consistent with the vcf.\n";
    print "--exomiserXLS Exomiser result files, separated with ',' it must be in XLS format.\n";
    print "--add_genotypeCall_flags=(Yes/No), anyt word starting with letter \"y\" will be considered as yes.\n";
    print "--AnnovarDIR=(where Annovar is installed and annovar_YX.sh is located.) \n";
	return(1);
}

sub GetGeneCoords {
	my ($gene_file) = @_;
	my %gene_coords;
	open INPUT, $genefile or die "Cannot open $genefile\n";
	while (my $Line = <INPUT>){
		chomp $Line;
		my @linesplit1 = split(/\t/,$Line);
		if($linesplit1[0] eq 'MT'){
			$linesplit1[0]='M'
		}
		my $chr="chr".$linesplit1[0];
		my $st=$linesplit1[1];
		my $end=$linesplit1[2];
		my $gen=$linesplit1[3]."\t".$linesplit1[4]; # gene name \t gene id
		$gene_coords{$chr}{$st}{$end} = $gen;
	}
	close INPUT;
	return \%gene_coords;
}

sub GetEnsSymbolMap{
	my ($gene_file) = @_;
	my %map;
	open INPUT, $genefile or die "Cannot open $genefile\n";
	while (my $Line = <INPUT>){
		chomp $Line;
		my @linesplit1 = split(/\t/,$Line);
		$map{$linesplit1[4]}=$linesplit1[3]; #ens_id = symbol
	}
	close INPUT;
	return \%map;
}

sub GetOmimAcc {
	my ($sharedfile) = @_;
	my %OmimAcc;
	open INPUT, $sharedfile or die "Cannot open $sharedfile\n";
	while (my $line = <INPUT>) {
		if ($line =~ m/^ENSG/) {
			chomp $line;
			my @elements = split(/\t/, $line ,-1);
			if ($elements[2] ne "") {
				$OmimAcc{$elements[0]} = $elements[2];
			}
		} else {
			next;
		}
	}
	close INPUT;
	return \%OmimAcc;
}

sub GetUniprotAcc {
	my ($sharedfile) = @_;
	my %UniprotAcc;
	open INPUT, $sharedfile or die "Cannot open $sharedfile\n";
	while (my $line = <INPUT>) {
		if ($line =~ m/^ENSG/) {
			chomp $line;
			my @elements = split(/\t/, $line, -1);
			if ($elements[3] ne "") {
				$UniprotAcc{$elements[0]} = $elements[3];
			}
		} else {
			next;
		}
	}
	close INPUT;
	return \%UniprotAcc;
}

sub GetInHouseMaf {
	my ($sharedfile) = @_;
	my %IHcontrols;
	open INPUT3, $sharedfile or die "Cannot open $sharedfile\n"; 
	shloop: while (my $Line = <INPUT3>){
		chomp $Line;
		my @linesplit2 = split(/\t/,$Line);
		my $c=$linesplit2[0];
		my $p=$linesplit2[1];
		my $r=$linesplit2[2];
		my $v=$linesplit2[3];
		my $maf=$linesplit2[7];
		$IHcontrols{$c}{$p}{$r}{$v}=$maf;
	}
	close INPUT3;
	return \%IHcontrols;
}

sub GetInHouseMafGATK {
	my ($sharedfile) = @_;
	my %IHcontrols;
	open INPUT3, $sharedfile or die "Cannot open $sharedfile\n"; 
	shloop: while (my $Line = <INPUT3>){
		chomp $Line;
		my @linesplit2 = split(/\t/,$Line);
		my $c=$linesplit2[0];
		my $p=$linesplit2[1];
		my $r=$linesplit2[2];
		my $v=$linesplit2[3];
		my $maf=$linesplit2[7];
		$IHcontrols{$c}{$p}{$r}{$v}=$maf;
	}
	close INPUT3;
	return \%IHcontrols;
}

sub GetOMIManno {
	my ($OMIMfile) = @_;
	my %OMIM;
	my %add_OMIM;
	my %seen_already;
	open INPUT, $OMIMfile or die "Cannot open $OMIMfile\n"; 
	while (my $Line = <INPUT>){
		chomp $Line;
		my @linesplit = split(/\t/,$Line);
		my $ens_id = $linesplit[0];
		my $go;
		my $wiki;
		my $MIM;
		my $oMIM;
		if ( exists $linesplit[2] ) {
			if ($linesplit[2] eq "") {
				$go = "NA";
			} else {
				$go = $linesplit[2];
			}
		} else {
			$go = "NA";
		}
			
		
		if ( exists $linesplit[3] ){
			if ($linesplit[3] eq "") {
				$wiki = "NA";
			} else {
				$wiki = $linesplit[3];
			}
		} else {
			$wiki = "NA";
		}
		
		if ( exists $linesplit[4] ){
			if ($linesplit[4] eq "") {
				$MIM = "NA";
			} else {
				$MIM = $linesplit[4];
			}
		} else {
			$MIM = "NA";
		}
		
		if ( exists $linesplit[5] ){
			if ($linesplit[5] eq "") {
				$oMIM = "NA";
			} else {
				$oMIM = $linesplit[5];
			}
		} else {
			$oMIM = "NA";
		}
		
		#group all individual terms per gene
		
		if(!exists $add_OMIM{$ens_id}){
			$add_OMIM{$ens_id}{'go'}="";
			$add_OMIM{$ens_id}{'wiki'}="";
			$add_OMIM{$ens_id}{'MIM'}="";
			$add_OMIM{$ens_id}{'oMIM'}="";
			}
		
		if(!exists $seen_already{$ens_id}{$go}){
		$add_OMIM{$ens_id}{'go'} = $add_OMIM{$ens_id}{'go'}.", ".$go;
		$seen_already{$ens_id}{$go}=0;
		}
		if(!exists $seen_already{$ens_id}{$wiki}){
		$add_OMIM{$ens_id}{'wiki'} = $add_OMIM{$ens_id}{'wiki'}.", ".$wiki;
		$seen_already{$ens_id}{$wiki}=0;
		}
		if(!exists $seen_already{$ens_id}{$MIM}){
		$add_OMIM{$ens_id}{'MIM'} = $add_OMIM{$ens_id}{'MIM'}.", ".$MIM;
		$seen_already{$ens_id}{$MIM}=0;
		}
		if(!exists $seen_already{$ens_id}{$oMIM}){
		$add_OMIM{$ens_id}{'oMIM'} = $add_OMIM{$ens_id}{'oMIM'}.", ".$oMIM;
		$seen_already{$ens_id}{$oMIM}=0;
		}
	}
	close INPUT;
		
	#loop to add all individual terms per ENSID together
	foreach my $e (keys %add_OMIM){
		$add_OMIM{$e}{'go'}=~s/^\,\s//;
		$add_OMIM{$e}{'wiki'}=~s/^\,\s//;
		$add_OMIM{$e}{'MIM'}=~s/^\,\s//;
		$add_OMIM{$e}{'oMIM'}=~s/^\,\s//;
		
		$OMIM{$e} = $add_OMIM{$e}{'go'}."\t".$add_OMIM{$e}{'wiki'}."\t".$add_OMIM{$e}{'MIM'}."\t".$add_OMIM{$e}{'oMIM'};
	}

	return \%OMIM;
}

sub GetIsInterestedGenes { # one ensembl gene_id on each line
	my ($InterestedGenefile) = @_;
	my %Mito=();
	if (-e $InterestedGenefile) {
		open MF, $InterestedGenefile or die "cannot open $InterestedGenefile";
		while (my $line = <MF>) {
			chomp $line;
			my @genam=split(/\t/,$line);
			if(!exists $Mito{$genam[4]}){
				$Mito{$genam[4]}=0;
			}	
		}
		close MF;
	} else {
		print "No insterested gene list is provided.\n";
	}
	return \%Mito;
}

sub AddGenoTypeToSampleCalls {
	my ($format, $sample_call) = @_;
	my $gene_type_call_qual = 20; ##### genotype call quality cut off
	my @format_fields = split(":", $format);
	my $GT_index;
	my $GQ_index;
	for (my $i=0;$i<scalar(@format_fields);$i++){
		if ($format_fields[$i] =~ m/GT/) {
			$GT_index = $i;
		}
		if ($format_fields[$i] =~ m/GQ/) {
			$GQ_index = $i;
		}
	}
	my $sample_call_processed=""; # vcf_sample1 \t geno_sample1 \t vcf_sample2 \t geno_sample2 ....
	my @sample_call_split = split("\t", $sample_call);
	for(my $i=0;$i<scalar(@sample_call_split);$i++){
		my $sample = $sample_call_split[$i];
		if ($sample !~ m/\.\/\./) {
			my @fields = split(":", $sample);
			my $GT = $fields[$GT_index];
			my $GQ = $fields[$GQ_index];
			if ($GQ < $gene_type_call_qual) {
				#$GT = "R/R";
				$GT = "NA"; #if low quality make a null call "NA"
			} elsif ($GT =~ m/0\/0/) {
				$GT = "R/R";
			} elsif ($GT =~ m/0\/[123456789]/) {
				$GT = "R/V";
			} elsif ($GT =~ m/([123456789])\/([123456789])/) {
				my $left = $1;
				my $right = $2;
				if($left!=$right) {
					$GT = "V1/V2"; 
				} else {
					$GT = "V/V";
				}
			}
			$sample_call_processed = $sample_call_processed."\t".$sample."\t".$GT;
		} else {
			$sample_call_processed = $sample_call_processed."\t".$sample."\t"."NA";
		}
	}
	$sample_call_processed =~ s/^\t//;
	return $sample_call_processed;
}

sub AddGenoTypeToSampleCalls_CompondHet {
	my ($first_V, $variants, $format, $sample_call) = @_;
	my $gene_type_call_qual = 13; ##### genotype call quality cut off
	my @format_fields = split(":", $format);
	my $GT_index;
	my $GQ_index;
	my $AD_index;
	for (my $i=0;$i<scalar @format_fields ;$i++){
		if ($format_fields[$i] =~ m/GT/) {
			$GT_index = $i;
		}
		if ($format_fields[$i] =~ m/GQ/) {
			$GQ_index = $i;
		}
		if ($format_fields[$i] =~ m/AD/) {
			$AD_index = $i;
		}		
	}
	my $sample_call_processed=""; # vcf_sample1 \t geno_sample1 \t vcf_sample2 \t geno_sample2 ....
	my @sample_call_split = split("\t", $sample_call);
	for(my $i=0;$i<scalar(@sample_call_split);$i++){
		my $sample = $sample_call_split[$i];
		if ($sample !~ m/\.\/\./) {
			my @fields = split(":", $sample);
			my $GT = $fields[$GT_index];
			my $GQ = $fields[$GQ_index];
			my $AD;
			if(defined $AD_index) {
				$AD = $fields[$AD_index];
			}
			
			if( $first_V eq $variants) {
				if ($GQ < $gene_type_call_qual) {
					#$GT = "R/R";
					$GT = "NA"; #low qual genotype ... make a null call "NA"
				} elsif ($GT =~ m/0\/0/) {
					$GT = "R/R";
				} elsif ($GT =~ m/0\/[123456789]/) {
					$GT = "R/V";
				} elsif ($GT =~ m/([123456789])\/([123456789])/) {
					my $left = $1;
					my $right = $2;
					if($left!=$right) {
					$GT = "V1/V2"; 
					} else {
					$GT = "V/V";
					}
				}
				$sample_call_processed = $sample_call_processed."\t".$sample."\t".$GT;
			} else {
				#determine which alternative allele is moved to the front
				my @V = split("," , $variants);
				my @AD_array;
				if(defined $AD_index) {
					@AD_array = split(",", $AD);
				}
				my $V_index;
				my $the_first_V_AD_index; 
				for(my $h=0;$h<scalar @V;$h++) {
					if($first_V eq $V[$h]) {
						$V_index = $h+1;
						$the_first_V_AD_index = $V_index+1;
					}
				}
				
				my %temp_hash;
				$temp_hash{$V_index} = 1;
				my $new_AD;
				if(defined $AD_index) {
					#print "sample_".$sample."\n";
					#print "AD_ARRAY_".$AD_array[0]."\n";
					#print "the first index _".$the_first_V_AD_index."\n";
					#print "ad array [index]_".$AD_array[$the_first_V_AD_index-1]."\n";
					if(scalar @AD_array == 1 && $AD_array[0] =~ m/\./ ) {
						$new_AD = '.';
					} else {
						$new_AD = $AD_array[0].",".$AD_array[$the_first_V_AD_index-1];
						for(my $h=1;$h < scalar @AD_array;$h++) {
							if($h != $the_first_V_AD_index-1) {
								$new_AD = $new_AD.",".$AD_array[$h];
							}
						}
					}
					
				}
				
				my $j = 2;
				for(my $h=1;$h<=scalar @V;$h++) {
					if ($h != $V_index) {
						$temp_hash{$h} = $j;
						$j+=1;
					}
				}
				
				#adding genotye
				if ($GQ < $gene_type_call_qual) {
					#$GT = "R/R";
					$GT = "NA"; #low qual null call!
				} elsif ($GT =~ m/0\/0/) {
					$GT = "R/R";
				} elsif ($GT =~ m/0\/[123456789]/) {
					$GT = "R/V";
				} elsif ($GT =~ m/([123456789])\/([123456789])/) {
					my $left = $1;
					my $right = $2;
					if($left!=$right) {
						$GT = "V1/V2"; 
					} else {
						$GT = "V/V";
					}
				}
				
				my $left_number;
				my $right_number;
				if ($fields[$GT_index] =~ m/(\d)\/(\d)/ ) {
					my $A_left = $1;
					my $A_right = $2;
					#print "left $A_left\n";
					#print "right $A_right\n";
					if ($A_left == $A_right && $A_left ==0 ) {
						$left_number = "0";
						$right_number = "0";
					} elsif ( $A_left ==0 ) {
						$left_number = 0;
						$right_number = $temp_hash{$A_right};
					} else {
						$left_number = $temp_hash{$A_left};
						$right_number = $temp_hash{$A_right};
					}
				}
				###replace sample GT numbers
				$sample = $left_number."/".$right_number;
				for(my $h=1;$h < scalar @fields;$h++) {
					if(defined $AD_index) {
						if($h == $AD_index) {
							$sample=$sample.":".$new_AD;
						} else {
							$sample=$sample.":".$fields[$h];
						}
					} else {
						$sample=$sample.":".$fields[$h];
					}
					my $teno = scalar @fields;
					#print "h is $h   total is $teno   fields[h] is $fields[$h]\n";
					#print "$sample\n";
				}
				
				$sample_call_processed = $sample_call_processed."\t".$sample."\t".$GT;
			}
		} else {
			$sample_call_processed = $sample_call_processed."\t".$sample."\t"."NA";
		}
	}
	$sample_call_processed =~ s/^\t//;
	return $sample_call_processed;
}

sub ExomiserXLS2Array {
	my ($xlsFile) = @_;
	my @output;
	#my $xlsFile = "/home/yaobo/Downloads/test.xls";
	my $parser = Spreadsheet::ParseExcel->new();
	my $workbook = $parser->parse($xlsFile);
	
	my @readinArray;
	
	die $parser->error(), ".\n" if ( !defined $workbook );
	
	# Following block is used to Iterate through all worksheets
	# in the workbook and print the worksheet content 
	
	my $i=1; #to read only the first worksheet
	for my $worksheet ( $workbook->worksheets() ) {
		# Find out the worksheet ranges
		my ( $row_min, $row_max ) = $worksheet->row_range();
		my ( $col_min, $col_max ) = $worksheet->col_range();
		my ($gene, $score_no_prior, $score_mouse_exists, $chr, $pos, $ref, $var, $nn_change, $pathpogenicity, $path_score, $mouse_pheno, $omim_anno);
		# $pos is as the wAnnovar format start position
		for my $row ( $row_min+1 .. $row_max ) {
			#print "Row $row:\n";
			for my $col ( $col_min .. $col_max ) {
				# Return the cell object at $row and $col
				my $cell = $worksheet->get_cell( $row, $col ); 
				if ($col == 0 ) {
					$gene = $cell->value();
				} elsif($col == 1) {
					$score_no_prior = $cell->value();
				} elsif($col == 2) {
					$score_mouse_exists = $cell->value();
				}elsif ($col ==3) {
					my @lines = split("\n", $cell->value());
					$nn_change = "";
					for (my $j=0;$j<scalar @lines;$j++) {
						if ($j == 0) {					
							my @temp = split(":", $lines[$j]);
							$chr = $temp[0];
							my @temp2 = split(/\./, $temp[1]);
							if ($temp2[1] =~ m/(\d+)(\w+)\>(\w+)\s\[(.*)\]/) {
								$pos = $1;
								$ref = "$2";
								$var = "$3";
							} elsif ($temp2[1] =~ m/(\d+)(\w+)\>(\-)\s\[(.*)\]/) { #deletion
								$pos = $1;
								$ref = "$2";
								$var = "$3";
							} elsif ($temp2[1] =~ m/(\d+)(\-)\>(\w+)\s\[(.*)\]/) { #insertion
								$pos = $1+1;  ###########################might be a bug of exomiser to record positions
								$ref = "$2";
								$var = "$3";
							} else {
								print "shit, a Ghost!!\n";
							}
						} elsif ($j>2) {
							chomp $lines[$j];
							if ($lines[$j] ne "" && $lines[$j] !~ m/^PHRED/) {
								$nn_change = $nn_change."\n".$lines[$j];
							}
						}
					}
					$nn_change =~ s/^\n//;
				} elsif($col==4) {
					my @lines = split("\n", $cell->value());
					$pathpogenicity = $lines[1];
					if ($lines[2] =~ m/(\d\.\d+)$/) {
						$path_score = $1;
					}
				} elsif($col==5) {
					my @lines = split("\n", $cell->value());
					#print scalar @lines."\n";
					$mouse_pheno = $lines[0];
					$mouse_pheno =~ s/No\sOMIM\sdisease\sentry//;
					#chomp $mouse_pheno;
					if (scalar @lines >= 2) {
						$omim_anno="";
						for (my $j=1;$j<scalar @lines;$j++) {
							#chomp $lines[$j];
							if ($lines[$j] =~ m/^OMIM\:/){
								$lines[$j] =~ s/^OMIM\://;
								$omim_anno = $omim_anno.";".$lines[$j];
								#chomp $omim_anno;
							}
							$lines[$j] =~ s/^OMIM\://;
							#chomp $omim_anno;
						}
						$omim_anno =~ s/^\;//;
					} else {
						$omim_anno="NA";
					}
				}
			}
			
			#print "$gene, $score_no_prior, $score_mouse_exists, $chr, $pos, $ref, $var, $nn_change, $pathpogenicity, $path_score, $mouse_pheno, $omim_anno\n";
			push @output, "$chr\t$pos\t$ref\t$var\t$gene\t$score_no_prior\t$score_mouse_exists\t$nn_change\t$pathpogenicity\t$path_score\t$mouse_pheno\t$omim_anno"
		}
		$i+=1;
		if ($i >=2 ){last;}
	}
	return \@output;
}

