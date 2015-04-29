#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use Math::CDF;
use Excel::Writer::XLSX;
use Spreadsheet::ParseExcel;

print "\n";
print "VCF_2_annotated_excel_20141010.pl says\n";
print "############################################################################################\n";
print "# YX: in house annotation part was amended from HG's                                       #\n";
print "# YX: the script will invoke Annovar and produce annotated output in Excel XLSX format     #\n";
print "# YX: input shoud be a single vcf file of a family or single sample                        #\n";                                 
print "# HG edits: to take (standalone) annovar output (with vcf file lines appended after)       #\n";
print "############################################################################################\n";
print "\n";
print "############################################################################################\n";
print "# Changes in this version:                                                                 #\n";
print "# *Major change*:                                                                          #\n";
print "# 1. CNVs can now be output to vcf files, but such vcf can contain only ONE sample!        #\n";
print "############################################################################################\n";

my $vcf_in;
my $InterestedGenefile=""; #Ensembl gene ids on each line
#my $exomiserXLSfiles="";
my $CNV_file="";
my $output_excel=""; #if empty then inHouse annotated file will not converted to excel file
my $output_excel_everything=""; # file to output table 'everything'
my $output_new_vcf=""; #specify the output of the new VCF file output
my $batch_MAF_file="";
my $reference_fasta="";
my $add_genotypeCall_flags="No";
my $help;
my $compount_var_supporting_reads_threshold = 5; #variants will have at least this number of supporting reads, otherwise will not count as a het/homo for a gene
my $annovar_dir ="/sharedlustre/IGM/annovar_2014jul14";

usage() if ( @ARGV < 1 || ! GetOptions('help|?' => \$help, "vcf=s"=>\$vcf_in, "InterestedGenes=s"=>\$InterestedGenefile, 'CNV=s' => \$CNV_file, 
			'out=s' => \$output_excel, 'outAll=s' => \$output_excel_everything, 'outVCF=s' => \$output_new_vcf,
			'add_genotypeCall_flags=s' => \$add_genotypeCall_flags, 
			'AnnovarDIR=s' => \$annovar_dir,  'batchMAF=s' => \$batch_MAF_file ) || defined $help );


unless (defined $vcf_in && -e $vcf_in) {
	die "You have not supplied input vcf file using --vcf or the file does not exist\n\n";
}


print "Start to run Annovar on $vcf_in\n";
`$annovar_dir/annovar_YX.sh $vcf_in`;
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
my $first_sample_column_index = 78; # change if add or reduced Annovar annotation 0_based
# it is the sample column not the genotype column eg: 0/0:..... not the R/R column
my @output_vcf_lines;
########### working here############
my $new_vcf_header=VCFHeader();

print "Looking for sample names:\n";
for(my $i=9;$i<scalar @bash_output_split;$i++){
	print "Found: ".$bash_output_split[$i]."\n";
	push @sample_names, $bash_output_split[$i];
	if($output_new_vcf ne "") {
		$new_vcf_header=$new_vcf_header."\t".$bash_output_split[$i];
	}
}

if($output_new_vcf ne "") {
	push @output_vcf_lines, $new_vcf_header."\n";
}

for(my $i=0;$i<scalar @sample_names;$i++){
	push @sample_columns, $first_sample_column_index+2*$i;
}


#### Requiring input in house db files to add as annotation
my $genefile="$annovar_dir/inHouse_db/Ensembl_Genes_GRCh37_75.txt";
my $sharedfile="$annovar_dir/inHouse_db/InHouse_OnTarget_Variants_MAFs.txt_418exomes";
my $sharedfile_GATK="$annovar_dir/inHouse_db/InHouse_OnTarget_GATK_MAFs.txt_179exomes";
my $OMIMfile="$annovar_dir/inHouse_db/Ensembl_GRCh37_75_OMIM_AllGeneNames.txt";
my $ens_gene_OMIM_Uniprot_Acc = "$annovar_dir/inHouse_db/ens_gene_symbol_omim_id_uniprot_id.txt";

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

if ($batch_MAF_file ne "") {
	print "Require db file: $batch_MAF_file\n";
	if (-e $batch_MAF_file) {
		print "Found $batch_MAF_file\n";
	} else {
		print "$batch_MAF_file not exist\n exit\n";
		exit;
	}
}


if ($reference_fasta ne "") {
	print "Require reference file: $reference_fasta\n";
	if (-e $reference_fasta) {
		print "Found $reference_fasta\n";
	} else {
		print "$reference_fasta not exist\n exit\n";
		exit;
	}
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
my %batchMAF;

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
my $batchMAF_ref;
if ($batch_MAF_file ne "") {
	my $batchMAF_ref = GetInHouseMafGATK($batch_MAF_file);
	%batchMAF = %$batchMAF_ref;
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
		
		my $vcf_chr=$elements[69];
		my $vcf_pos=$elements[70];
		my $vcf_R=$elements[72];
		my @temp_temp_temp=split(/\,/, $elements[73]);
		my $vcf_A=$temp_temp_temp[0];
		 		
		my $FORMAT=$elements[77];
		my $Sample_Call="";
		for(my $i=78;$i<scalar(@elements);$i++) {
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
		my $query_chr;
		if ($vcf_chr !~ m/^chr/) {
			$query_chr = "chr".$vcf_chr;
		} else {
			$query_chr = $vcf_chr;
		}
		
		if ($query_chr =~ m/chrMT/i) {
			$query_chr = "chrM";
		}
		
		if (exists $inHouse_MAF{$query_chr}{$pos1}{$ref1}{$alt1} ) {
			$maf = $inHouse_MAF{$query_chr}{$pos1}{$ref1}{$alt1};
		} else {
			$maf = "0";
		}

		#find in house 'GATK' MAF and batch MAF for variants
		############################################## need to map v with Annovar format
		my $maf_GATK = "";
		if (exists $inHouse_MAF_GATK{$annovar_chr}{$annovar_pos}{$annovar_R}{$annovar_A}) {
			$maf_GATK = $inHouse_MAF_GATK{$annovar_chr}{$annovar_pos}{$annovar_R}{$annovar_A};
		} else {
			$maf_GATK = "0";
		}
		
		my $maf_batch = "";
		if ($batch_MAF_file ne "") {
			if (exists $batchMAF{$annovar_chr}{$annovar_pos}{$annovar_R}{$annovar_A}) {
				$maf_batch = $batchMAF{$annovar_chr}{$annovar_pos}{$annovar_R}{$annovar_A};
			} else {
				$maf_batch = "0";
			}
		}
		############################################## vcf format changed done #######
		my $output_vcf_string;
		if ($output_new_vcf ne "") {
			$output_vcf_string="$vcf_chr\t$vcf_pos\t.\t$vcf_R\t$elements[73]\t$elements[74]\t$elements[75]\t$elements[76]";
			$output_vcf_string=$output_vcf_string.";ANNO_INT=$isInterested";
			if($elements[56] eq "NA" ) { #PopFreqMax
				$output_vcf_string=$output_vcf_string.";ANNO_POPMAX=0";
			} else {
				$output_vcf_string=$output_vcf_string.";ANNO_POPMAX=".sprintf( "%.4f",$elements[56]);
			}
			if($elements[54] eq "NA" ) { #CADD gt 10
				$output_vcf_string=$output_vcf_string.";ANNO_CADD=0";
			} else {
				$output_vcf_string=$output_vcf_string.";ANNO_CADD=".sprintf( "%.4f",$elements[54]);
			}
				
		}
		my $Sample_Call_processed = AddGenoTypeToSampleCalls($FORMAT, $Sample_Call);
		for (my $i=0;$i<=77;$i++){ push @output, $elements[$i]."\t";}
		if ($batch_MAF_file ne "") {
			push @output, $Sample_Call_processed."\t".$SNPorINDEL ."\t".$gene_name."\t".$ens_id."\t".$isInterested."\t"
			.$maf_batch."\t".$maf_GATK."\t".$omim_anno."\t".$genecard_link."\t".$omim_link."\t".$uniprot_link."\n";
			if ($output_new_vcf ne "") {
				$output_vcf_string=$output_vcf_string.";ANNO_B_MAF=".sprintf( "%.4f",$maf_batch);
			}
		} else {
			push @output, $Sample_Call_processed."\t".$SNPorINDEL ."\t".$gene_name."\t".$ens_id."\t".$isInterested."\t"
			.$maf."\t".$maf_GATK."\t".$omim_anno."\t".$genecard_link."\t".$omim_link."\t".$uniprot_link."\n";
			if ($output_new_vcf ne "") {
				$output_vcf_string=$output_vcf_string.";ANNO_InHo_MAF=".sprintf( "%.4f",$maf);
			}
		}
		if ($output_new_vcf ne "") {
			$output_vcf_string=$output_vcf_string.";ANNO_InHoG_MAF=".sprintf( "%.4f",$maf_GATK);
			$output_vcf_string=$output_vcf_string."\t$elements[77]";
			for( my $i=0; $i < scalar @sample_columns;$i++) {
				$output_vcf_string=$output_vcf_string."\t".$elements[$first_sample_column_index+$i];
			}
			$output_vcf_string=$output_vcf_string."\n";
			push @output_vcf_lines, $output_vcf_string;
		}
	}
}
close (VCF);

open OUTPUT, ">$output_file" or die "Cannot open file $output_file to output. \n";
print OUTPUT @output;
close OUTPUT;

if ($output_new_vcf ne "") {
	print "Output customized VCF files to $output_new_vcf....";
	open OUTPUTVCF, ">$output_new_vcf" or die "Cannot open file $output_new_vcf to output.\n";
	print OUTPUTVCF @output_vcf_lines;
	close OUTPUTVCF;
	print "Done!\n";
}

############### when CNV is provided, export an extra vcf with CNV added for EACH sample!
if ($CNV_file ne "" && $output_new_vcf ne "") { 
	print "CNV file provided:$CNV_file, start to output a vcf with CNV inserted.\n";
	### process the CNV hash and output vcf array for each sample
	# read in SNP and indels into the CNV hash
	my %vcf_with_CNV_out;
	for(my $i=0;$i<scalar @sample_names;$i++){
		# reorder the vcf of mixture of snp indel and CNV
		foreach my $vcf_line (@output_vcf_lines) {
			chomp $vcf_line;
			if ($vcf_line !~ m/^\#/) {
				my @temp_temp_temp = split("\t", $vcf_line);
				my $sample_column = 9+$i;
				my $new_vcf_line ="$temp_temp_temp[0]\t$temp_temp_temp[1]\t$temp_temp_temp[2]\t$temp_temp_temp[3]\t$temp_temp_temp[4]\t$temp_temp_temp[5]\t$temp_temp_temp[6]\t$temp_temp_temp[7]\t$temp_temp_temp[8]\t$temp_temp_temp[$sample_column]" ;
				if (! exists $vcf_with_CNV_out{$sample_names[$i]}{$temp_temp_temp[0]}{$temp_temp_temp[1]}) {
					$vcf_with_CNV_out{$sample_names[$i]}{$temp_temp_temp[0]}{$temp_temp_temp[1]} = $new_vcf_line;
				} else {
					$vcf_with_CNV_out{$sample_names[$i]}{$temp_temp_temp[0]}{$temp_temp_temp[1]} = 
					$vcf_with_CNV_out{$sample_names[$i]}{$temp_temp_temp[0]}{$temp_temp_temp[1]}."\n".$new_vcf_line;
				}
				
			}
		}
	}
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
			my $start=$elements[5];
			my $end=$elements[6];
			$elements[7] =~ s/"//g;
			my $chr="chr".$elements[7];
			# in order to avoid identical start point with a snp/indel
			# if CNV start is smaller than 50, then the start will add 1, minus 1 othewise
			# CNV exome sequencing data does not have accurate bondaries anyway.
			#if ($start <= 50) {
			#	while (exists $vcf_with_CNV_out{$id}{$chr}{$start}) {
			#		$start = $start + 1;
			#	}
			#} else {
			#	while (exists $vcf_with_CNV_out{$id}{$chr}{$start}) {
			#		$start = $start - 1;
			#	}
			#}
			
			my $ref = `sh $annovar_dir/py_Get_the_allele.sh $chr $start`;
			chomp $ref;
			if ($ref !~ m/A|C|G|T|N/i) {
				print "Reference returned from py_Get_the_allele.sh is $ref, please make sure the right reference fasta is used.\n";
				exit;
			}
			my ($alt, $svtype, $svlen, $CN);
			if ( $elements[3] eq "deletion"  ) { #select which CNVs are to used to determine Compount Het
				$alt='<CNV>';
				$svtype="CNV";
				$svlen=$start-$end;
				$CN=1;
			} else {
				$alt='<CNV>';
				$svtype="CNV";
				$svlen=$end-$start;
				$CN=3;
			}
			
			### Note: Genotype of CNV line is absent as ./.
			### Genotyping quality score is set to 31
			# ANNO_CNV_BF=$elements[9];
			# ANNO_CNV_EXP_READS=$elements[10];
			# ANNO_CNV_OBS_READS=$elements[11];
			# ANNO_CNV_RATIO=$elements[12];
			# ANNO_TRUE_START=$elements[5]
			my $temp_string="$chr\t$start\t.\t$ref\t$alt\t31\tPASS\tIMPRECISE;SVTYPE=$svtype;END=$end;SVLEN=$svlen;ANNO_CNV_BF=$elements[9];ANNO_CNV_EXP_READS=$elements[10];ANNO_CNV_OBS_READS=$elements[11];ANNO_CNV_RATIO=$elements[12];ANNO_TRUE_START=$elements[5]\tGT:GQ:CN\t./.:31:$CN";
			if (! exists $vcf_with_CNV_out{$id}{$chr}{$start}) {
				$vcf_with_CNV_out{$id}{$chr}{$start} = $temp_string;
			} else {
				$vcf_with_CNV_out{$id}{$chr}{$start} = 
				$vcf_with_CNV_out{$id}{$chr}{$start}."\n".$temp_string;
			}
			
			unless( exists $CNV_ids{$id} ) {
				$CNV_ids{$id} =0;
				print "The file contains info for $id.\n";
			}
		}
	}
	###sort and output vcf file for each sample
	
	### This is script works with human reference only
	my @chromosomes = ("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11",
	"chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", 
	"chrX", "chrY", "chrM" );
	for(my $i=0;$i<scalar @sample_names;$i++){
		my $cnv_vcf_file = $output_new_vcf."."."$sample_names[$i]".".withCNV.vcf";
		open OUTPUT_CNVVCF, ">$cnv_vcf_file" or die "Cannot open file $cnv_vcf_file to output.\n";
		print OUTPUT_CNVVCF CNV_VCFHeader();
		print OUTPUT_CNVVCF "\t".$sample_names[$i]."\n";
		foreach my $chr_here (@chromosomes) {
			foreach my $start ( sort {$a<=>$b} keys %{$vcf_with_CNV_out{$sample_names[$i]}{$chr_here}} ) {
				print OUTPUT_CNVVCF $vcf_with_CNV_out{$sample_names[$i]}{$chr_here}{$start}."\n";
			}
		}
		close OUTPUT_CNVVCF;
		print "Done for sample $sample_names[$i]!\n";
	}
}
close CNV;

################################################
#
#  output to excel part
#
#


if ($output_excel eq "") {print "#####  No excel file is required. ######\nDone!\n"; exit;}

if(scalar @sample_columns != scalar @sample_names) {
	die "Number of sample names and columns are inconsistent.\n";
}

for(my $i=0;$i<scalar @sample_columns;$i++){
	print "Found sample: $sample_names[$i] on column $sample_columns[$i]\n";
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
my %CNV_regions_compound; # it is used to recored a selection of (Deletions only or all) CNV thus used to determine compound het;
#my %gene_name_count; #it's for selecting Xlinked rare variants
my $ArrayOfExomiserXLS_ref;
my @ArrayOfExomiserXLS;
my %HashOfExomiserXLS;
my @ArrayOfExomiserXLS_Output;

# To find extended compound hetero: multiple variant hits and CNVs on the same gene 
my %all_gene_id; # to record gene_id => gene_name
my %CNV_on_samples; # to record all CNVs in a new manner [Gene][sample]=cnv_id
my %Variants_on_samples; # to record all rare variants;
my %Variants_on_samples_count; # to record all rare variants in a new manner [Gene][sample]=chr_v-position_genoType_alleleFreq(ESP)_alleleFreq(1000G)_alleleFreq(in house);

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
			$CNV_regions{$id}{$chr}{$start}{$end} = $elements[8];
			if ( $elements[3] eq "deletion"  ) { #select which CNVs are to used to determine Compount Het
				$CNV_regions_compound{$id}{$chr}{$start}{$end} = $elements[3]."\t".$elements[7]."\t".$elements[8];
			}
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
		my $format = $elements[77];
		my @format_fields = split(":", $format);
		my $GQ_index;
		my $AD_index;
		my $GT_index;
		my $DP_index;
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
				my $query_chr;
				if ($elements[69] !~ m/^chr/) {
					$query_chr = "chr".$elements[69];
				} else {
					$query_chr = $elements[69];
				}
				
				if ($query_chr =~ m/chrMT/i) {
					$query_chr = "chrM";
				}
				if (exists $CNV_regions{$sample_names[$i]}) {
					start_loop: foreach my $start ( sort {$a<=>$b} keys %{$CNV_regions{$sample_names[$i]}{$query_chr}} ) {
						foreach my $end (keys %{$CNV_regions{$sample_names[$i]}{$query_chr}{$start}}) {
							if ( $elements[70] >= $start && $elements[70] <= $end) {
								$on_CNV=$CNV_regions{$sample_names[$i]}{$query_chr}{$start}{$end};
								last start_loop;
							}
						}
					}
				} else {
					$on_CNV="NA";
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
						if ($GT =~ m/0\/[123456789]|1\/1|2\/2|3\/3|4\/4|5\/5|6\/6|7\/7/ && $AD ne ".") {
							my @bases=split("/", $GT);
							my @allele_cov=split(",", $AD);
						
							if ($allele_cov[0] == 0 && $allele_cov[$bases[1]] == 0) {
								$geno_flag = "E";
								#print "#WARNING: Call made on 0 depth is found, possible bug of GATK VQSR. the line is :\n$line\n";
							} else {
								my $test_v = $allele_cov[$bases[1]];
								$test_v = $allele_cov[0] if $allele_cov[0] < $allele_cov[$bases[1]];
								$p_bino = 2*&Math::CDF::pbinom($test_v, $allele_cov[0]+$allele_cov[$bases[1]], 0.5);
								if ($GT =~ m/0\/[123456789]/) {
									if ($p_bino < 0.05) {
										$geno_flag = "B"; #bad call for hetero
									} else {
										if ($allele_cov[0]+$allele_cov[$bases[1]] >= 10) {
											$geno_flag = "G";
										} else {
											$geno_flag = "B";
										}
									}
								} else {
									if ($p_bino <= 0.05) {
										if ($allele_cov[0]+$allele_cov[$bases[1]] >= 10) {
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
		if($elements[57] eq "NA" || $elements[57] eq "." ) {#1000G
			$elements[57] = 0;
		}
		if($elements[62] eq "NA" || $elements[62] eq ".") { #esp6500
			$elements[62] = 0;
		}
		if($elements[24] eq "NA" || $elements[24] eq "." ) { #cg69
			$elements[24] = 0;
		}
		if($elements[56] eq "NA" ) { #PopFreqMax
			$elements[56] = 0;
		}
		if($elements[54] eq "NA" ) { #CADD gt 10
			$elements[54] = 0;
		}
		
		#Change in-house MAF and GATK MAF 'NA' to '0'
		if ($elements[$length-8] eq "NA") { #inHouse maf
			$elements[$length-8] = 0;
		}
		if ($elements[$length-7] eq "NA") { #inhouse gatk maf
			$elements[$length-7] = 0;
		}

#my @header = ('Chr','Start','Ref','Obs','Func.ensGene','ExonicFunc.ensGene','AAChange.ensGene','GeneName', 'caddgt10', 'Predictions');

#my @header_2 = ('InHouseMAF_GATK','OMIM Link','GeneCard Link','Uniprot Link','OMIM_Gene_Description','WikiGene_Description','GoTerm', 'isInterested',
#'Format',);

#my @header_3 = ('Func.knownGene','ExonicFunc.knownGene','AAChange.knownGene','Func.refGene','ExonicFunc.refGene','AAChange.refGene', 'PopFreqMax',
#'InHouseMAF', 'SNPorINDEL', 'dbSNP138', 'VariantCall_quality', 'GeneID', 'genomicSuperDups', 'ESP6500_ALL','1000g2012feb_ALL', 'cg69',
#'clinvar_20140303', 'gwasCatalog', 'tfbsConsSites', 'cosmic68',
#'phastConsElements46way', 'LJB23_SIFT_score_converted', 'LJB23_Polyphen2_HDIV_score', 'LJB23_Polyphen2_HDIV_pred',
#'LJB23_Polyphen2_HVAR_score', 'LJB23_Polyphen2_HVAR_pred', 'LJB23_LRT_score_converted', 'LJB23_LRT_pred', 
#'LJB23_MutationTaster_score_converted', 'LJB23_MutationTaster_pred', 'LJB23_MutationAssessor_score_converted', 
#'LJB23_MutationAssessor_pred', 'LJB23_FATHMM_score_converted', 'LJB23_FATHMM_pred', 'LJB23_RadialSVM_score_converted', 
#'LJB23_RadialSVM_pred', 'LJB23_LR_score', 'LJB23_LR_pred', 'LJB23_GERP++', 'LJB23_PhyloP', 'LJB23_SiPhy');

#my $superDupsRow=9;
#$superDupsRow= $superDupsRow + 22;

		#re-arrange columns
		my $string_1 = $elements[69]."\t".$elements[70]."\t".$elements[72]."\t".$elements[73]."\t".$elements[5]."\t".$elements[8]."\t".$elements[9]."\t"
						.$elements[$length-11]."\t".$elements[54]."\t".$elements[30].":".$elements[32].":".$elements[35].":".$elements[38].":"
						.$elements[41].":".$elements[44].":".$elements[47].":".$elements[49]."\t";
		
		my $string_2 =	$elements[$length-7]."\t".$elements[$length-1]."\t".$elements[$length-2]."\t".$elements[$length]."\t".$elements[$length-3]."\t"
						.$elements[$length-5]."\t".$elements[$length-6]."\t".$elements[$length-9]."\t".$elements[77]."\t";	
		
		my $string_3 = $elements[10]."\t".$elements[13]."\t".$elements[14]."\t".$elements[15]."\t".$elements[18]."\t".$elements[19]."\t".$elements[56]
						."\t".$elements[$length-8]."\t".$elements[$length-12]."\t".$elements[25]."\t".$elements[74]."\t".$elements[$length-10]."\t".
						$elements[21]."\t".$elements[62]."\t".$elements[57]."\t".$elements[24]."\t".$elements[53]."\t".$elements[66]."\t".$elements[68]
						."\t".$elements[55]."\t".$elements[20]."\t".$elements[27]."\t".$elements[29]."\t".$elements[30]."\t".$elements[31]."\t".$elements[32]
						."\t".$elements[34]."\t".$elements[35]."\t".$elements[37]."\t".$elements[38]."\t".$elements[40]."\t".$elements[41]."\t".$elements[43]
						."\t".$elements[44]."\t".$elements[46]."\t".$elements[47]."\t".$elements[48]."\t".$elements[49]."\t".$elements[50]."\t".$elements[51]
						."\t".$elements[52]; #re-arrange columns
		
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
		
		if($elements[8] eq "synonymous SNV" && $elements[13] eq "synonymous SNV" && $elements[18] eq "synonymous SNV"  ) { # eliminate synonymous SNV
			next;
		}
		if($elements[5] =~ m/(intronic|intergenic|upstream|downstream)/i &&  $elements[10] =~ m/(intronic|intergenic|upstream|downstream)/i 
						&& $elements[15] =~ m/(intronic|intergenic|upstream|downstream)/i  ) { # eliminate intronic or intergenic variants
			next;
		}
		
		if($rr_na_sum < scalar @sample_names) { # RR and NA are removed
			push @output_all, $sample_string_final;
			
			# for filtered spreadsheet

			if( $elements[56] <= 0.05 && $elements[$length-7] <= 0.05 && $elements[$length-8] <= 0.25) { #batch MAF filtering is very relax
				push @output_filtered, $sample_string_final;
				
				# for deleterious hits spreadsheet
				my $number_deleterious_predictions=0;
				if($elements[54] >= 0) { # use CADD most great 10 percent  scores
					$number_deleterious_predictions+=1;
				}
				if($elements[30] eq "D" || $elements[30] eq "P") { # use Polyphen2 HDIV prediction
					$number_deleterious_predictions+=1;
				}
				if($elements[32] eq "D" || $elements[32] eq "P") { # use Polyphen2 HVAR prediction
					$number_deleterious_predictions+=1;
				}
				if($elements[35] eq "D") { # use LRT prediction
					$number_deleterious_predictions+=1;
				}
				if($elements[38] eq "A" || $elements[38] eq "D") { # use MutationTaster prediction
					$number_deleterious_predictions+=1;
				}
				if($elements[41] =~ m/H/ || $elements[41] eq "M") { # use MutationAccessor prediction
					$number_deleterious_predictions+=1;
				}
				if($elements[44] eq "D") { # use FATHMM prediction
					$number_deleterious_predictions+=1;
				}
				if($elements[47] eq "D") { # use MetaSVM prediction
					$number_deleterious_predictions+=1;
				}
				if($elements[49] eq "D") { # use MetaLR prediction
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
					my $sample = $elements[$sample_columns[$i]]; 
					my $var_supporting_reads = 0;
					my ($GT,$AD);
					#print $sample."\n";
					if ($sample !~ m/^\./) {
						my @fields = split(":", $sample);
						$GT = $fields[$GT_index];
						$AD = $fields[$AD_index];
						my @bases=split("/", $GT);
						my @allele_cov=split(",", $AD);
						$var_supporting_reads = $allele_cov[$bases[1]];
						#print "GT:$GT AD:$AD bases1:$bases[1] var_supporting_reads:$var_supporting_reads"."\n";
					} #else {
						#print "var_supporting_reads:$var_supporting_reads"."\n";
					#}
					
					# extend the hash to record variants for samples and genes for %Variants_on_samples: $variants_on_samples{sample_name}{gene_id}=chr_pos_genoType
					unless ($elements[$length-10] eq "NA" || exists $all_gene_id{$elements[$length-10]}) { 
						$all_gene_id{$elements[$length-10]} = $elements[$length-11]."\t".$elements[$length-9]."\t".
						$elements[$length-6]."\t".$elements[$length-5]."\t".$elements[$length-4]."\t".$elements[$length-3]."\t".
						$elements[$length-2]."\t".$elements[$length-1]."\t".$elements[$length]; 
					}
					if ( exists $Variants_on_samples{$sample_names[$i]}{$elements[$length-10]}  && $elements[$sample_col+1] ne "R/R" 
						&& $elements[$sample_col+1] ne "NA" && $var_supporting_reads >= $compount_var_supporting_reads_threshold) {
						$Variants_on_samples{$sample_names[$i]}{$elements[$length-10]} = 
						$Variants_on_samples{$sample_names[$i]}{$elements[$length-10]}."\n".
						$sample_string_final;
					} elsif( $elements[$length-10] ne "NA" && $elements[$sample_col+1] ne "R/R" && $elements[$sample_col+1] ne "NA" 
						&& $var_supporting_reads >= $compount_var_supporting_reads_threshold ) {
						$Variants_on_samples{$sample_names[$i]}{$elements[$length-10]} = $sample_string_final;				
					}
					
					if ( exists $Variants_on_samples_count{$sample_names[$i]}{$elements[$length-10]}  && $elements[$sample_col+1] ne "R/R" 
						&& $elements[$sample_col+1] ne "NA" && $var_supporting_reads >= $compount_var_supporting_reads_threshold) {
						$Variants_on_samples_count{$sample_names[$i]}{$elements[$length-10]} = 
						$Variants_on_samples_count{$sample_names[$i]}{$elements[$length-10]}.";".
						$elements[69]."_".$elements[70]."_".$elements[$sample_col+1]."_".sprintf( "%.4f",$elements[56])."_".sprintf( "%.4f",$elements[$length-8]);
					} elsif( $elements[$length-10] ne "NA" && $elements[$sample_col+1] ne "R/R" && $elements[$sample_col+1] ne "NA" 
						&& $var_supporting_reads >= $compount_var_supporting_reads_threshold ) {
						$Variants_on_samples_count{$sample_names[$i]}{$elements[$length-10]} = 
						$elements[69]."_".$elements[70]."_".$elements[$sample_col+1]."_".sprintf( "%.4f",$elements[56])."_".sprintf( "%.4f",$elements[$length-8]);				
					}
				}
				
				if ($elements[$length-9] eq "YES" && $number_deleterious_predictions >=1) {  #those that were predicted as “Deleterious” by at least two predictors
					push @output_possibleHits, $sample_string_final;
				} elsif ($elements[$length-9] eq "NO" && $number_deleterious_predictions >=4) { #those that were predicted as “Deleterious” by half of the predictors
					push @output_possibleHits, $sample_string_final;
				} elsif($elements[$length-12] eq "SNP" && $number_CompoundHet >=1) { #Compound heterozygous “V1/V2”* of single base polymorphisms
					push @output_possibleHits, $sample_string_final;
				} elsif($elements[$length-9] eq "YES" && $elements[56] <= 0.01) { # those that have AAF less than 0.01 in popFreqMax 
					push @output_possibleHits, $sample_string_final;
				}
				
				if ($elements[$length-9] eq "YES" && $elements[56] <= 0.01) {
					push @very_rare_interested_vars, $sample_string_final;
				}
				
				#for Xlinked rare variants
				if ($elements[69] =~ m/^(chrx|X)/i) {
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
					if ( $elements[69] !~ m/[XYMxym]/ ) {
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
	foreach my $sample (keys %CNV_regions_compound) {
		foreach my $chr (keys %{$CNV_regions_compound{$sample}}) {
			foreach my $start_CNV (keys %{$CNV_regions_compound{$sample}{$chr}}) {
				foreach my $end_CNV (keys %{$CNV_regions_compound{$sample}{$chr}{$start_CNV}}) {
					start_loop: foreach my $start_gene ( sort {$a<=>$b} keys %{$en_genes{$chr}} ) {
						foreach my $end_gene (keys %{$en_genes{$chr}{$start_gene}}) {
							if ( $start_gene <= $end_CNV && $end_gene >= $start_CNV ) {
								my @temp = split("\t", $en_genes{$chr}{$start_gene}{$end_gene});
								my $gene_name = $temp[0];
								my $gene_id = $temp[1];
								my $isOnInterestedGenes = "NO";
								unless ( $gene_id eq "NA" || exists $all_gene_id{$gene_id}) {
									### to add more annotation to the genes as requested
									my $temp_string;
									if ( ! exists $isInterestedGenes{$gene_id}) {
										$temp_string = $gene_name."\t"."NO";
									} else {
										$temp_string = $gene_name."\t"."YES";
										$isOnInterestedGenes = "YES";
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
									$CNV_on_samples{$sample}{$gene_id}=$CNV_on_samples{$sample}{$gene_id}.";".$CNV_regions_compound{$sample}{$chr}{$start_CNV}{$end_CNV}."\t".$isOnInterestedGenes;
								} else {
									$CNV_on_samples{$sample}{$gene_id} = $CNV_regions_compound{$sample}{$chr}{$start_CNV}{$end_CNV}."\t".$isOnInterestedGenes;
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

my @header = ('Chr','Start','Ref','Obs','Func.ensGene','ExonicFunc.ensGene','AAChange.ensGene','GeneName', 'caddgt10', 'Predictions');

my @header_2 = ('InHouseMAF_GATK','OMIM Link','GeneCard Link','Uniprot Link','OMIM_Gene_Description','WikiGene_Description','GoTerm', 'isInterested',
'Format',);

my @header_3 = ('Func.knownGene','ExonicFunc.knownGene','AAChange.knownGene','Func.refGene','ExonicFunc.refGene','AAChange.refGene', 'PopFreqMax',
'InHouseMAF', 'SNPorINDEL', 'dbSNP138', 'VariantCall_quality', 'GeneID', 'genomicSuperDups', 'ESP6500_ALL','1000g2012feb_ALL', 'cg69',
'clinvar_20140303', 'gwasCatalog', 'tfbsConsSites', 'cosmic68',
'phastConsElements46way', 'LJB23_SIFT_score_converted', 'LJB23_Polyphen2_HDIV_score', 'LJB23_Polyphen2_HDIV_pred',
'LJB23_Polyphen2_HVAR_score', 'LJB23_Polyphen2_HVAR_pred', 'LJB23_LRT_score_converted', 'LJB23_LRT_pred', 
'LJB23_MutationTaster_score_converted', 'LJB23_MutationTaster_pred', 'LJB23_MutationAssessor_score_converted', 
'LJB23_MutationAssessor_pred', 'LJB23_FATHMM_score_converted', 'LJB23_FATHMM_pred', 'LJB23_RadialSVM_score_converted', 
'LJB23_RadialSVM_pred', 'LJB23_LR_score', 'LJB23_LR_pred', 'LJB23_GERP++', 'LJB23_PhyloP', 'LJB23_SiPhy');

if ($batch_MAF_file ne "") {
	$header_3[7] = "Batch_MAF";
}

my $superDupsRow=9;

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
$superDupsRow= $superDupsRow + 22;

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


############################# output compound heterozygous sheet ###################################################
my @compound_header = ('Sample', 'Gene_ID', 'Count', 'Helper');
my $start_output_index=3; # the col number for CNVs info to start to add NA
foreach my $sample (@sample_names) {
	push @compound_header, "in $sample";
	$start_output_index += 1;
}
push @compound_header, @header;
my $worksheet_compound = $workbook->add_worksheet('Compound_Heterozygous');
$worksheet_compound -> write_row(0,0,\@compound_header, $format_header);


my %gene_sample_CNV_count; #number of how many CNVs on the gene
my %gene_sample_VAR_count; #number of how many snp/indel on the gene

foreach  my $sample (@sample_names) {
	if (exists $CNV_on_samples{$sample}) {
		foreach my $gene_id (keys %{$CNV_on_samples{$sample}}) {
			my @temp = split(";",$CNV_on_samples{$sample}{$gene_id});
			if(exists $gene_sample_CNV_count{$sample}{$gene_id}) {
				$gene_sample_CNV_count{$sample}{$gene_id} = $gene_sample_CNV_count{$sample}{$gene_id} + scalar @temp;
			} else {
				$gene_sample_CNV_count{$sample}{$gene_id} = scalar @temp;
			}
		}
	}
}

foreach  my $sample (@sample_names) {
	foreach my $gene_id (keys %{$Variants_on_samples_count{$sample}}) {
		$gene_sample_VAR_count{$sample}{$gene_id} = 0;
		my @temp = split(";",$Variants_on_samples_count{$sample}{$gene_id});
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

# just to count how many rows need to output
my $range=1;
foreach my $gene_id (keys %all_gene_id) {
	my $output_light = "red";
	foreach my $sample (@sample_names) {
		if( exists $CNV_on_samples{$sample}{$gene_id} && exists $Variants_on_samples{$sample}{$gene_id}) {
			my @gene_sample_var = split("\n", $Variants_on_samples{$sample}{$gene_id});
			my @gene_sample_cnv = split(";", $CNV_on_samples{$sample}{$gene_id});
			if ($gene_sample_CNV_count{$sample}{$gene_id} + $gene_sample_VAR_count{$sample}{$gene_id} >= 2) {
				$output_light = "green";
			}
		} elsif (exists $CNV_on_samples{$sample}{$gene_id} && !exists $Variants_on_samples{$sample}{$gene_id}) {
			my @gene_sample_cnv = split(";", $CNV_on_samples{$sample}{$gene_id});
			if ($gene_sample_CNV_count{$sample}{$gene_id} >= 2) {
				$output_light = "green";
			}
		} elsif(!exists $CNV_on_samples{$sample}{$gene_id} && exists $Variants_on_samples{$sample}{$gene_id}) {
			my @gene_sample_var = split("\n", $Variants_on_samples{$sample}{$gene_id});
			if ($gene_sample_VAR_count{$sample}{$gene_id} >= 2) {
				$output_light = "green";
			}
		}
		
		if ($output_light eq "green") {
			if( exists $Variants_on_samples{$sample}{$gene_id}) {
				my @gene_sample_var = split("\n", $Variants_on_samples{$sample}{$gene_id});
				foreach my $var_string (@gene_sample_var) {
					$range += 1;
				}
			}
			if( exists $CNV_on_samples{$sample}{$gene_id}) {
				my @gene_sample_cnv = split(";", $CNV_on_samples{$sample}{$gene_id});
				foreach my $cnv_string (@gene_sample_cnv) {
					$range += 1;
				}
			}
		}
	}
}


my $excel_line_index=1;
foreach my $gene_id (keys %all_gene_id) {
	my $output_light = "red";
	foreach my $sample (@sample_names) {
		if( exists $CNV_on_samples{$sample}{$gene_id} && exists $Variants_on_samples{$sample}{$gene_id}) {
			my @gene_sample_var = split("\n", $Variants_on_samples{$sample}{$gene_id});
			my @gene_sample_cnv = split(";", $CNV_on_samples{$sample}{$gene_id});
			if ($gene_sample_CNV_count{$sample}{$gene_id} + $gene_sample_VAR_count{$sample}{$gene_id} >= 2) {
				$output_light = "green";
			}
		} elsif (exists $CNV_on_samples{$sample}{$gene_id} && !exists $Variants_on_samples{$sample}{$gene_id}) {
			my @gene_sample_cnv = split(";", $CNV_on_samples{$sample}{$gene_id});
			if ($gene_sample_CNV_count{$sample}{$gene_id} >= 2) {
				$output_light = "green";
			}
		} elsif(!exists $CNV_on_samples{$sample}{$gene_id} && exists $Variants_on_samples{$sample}{$gene_id}) {
			my @gene_sample_var = split("\n", $Variants_on_samples{$sample}{$gene_id});
			if ($gene_sample_VAR_count{$sample}{$gene_id} >= 2) {
				$output_light = "green";
			}
		}
		
		#output
		if ($output_light eq "green") {
			if( exists $Variants_on_samples{$sample}{$gene_id}) {
				my @gene_sample_var = split("\n", $Variants_on_samples{$sample}{$gene_id});
				foreach my $var_string (@gene_sample_var) {
					#print "here is the string: $var_string\n";
					my @output;
					push @output, $sample;
					push @output, $gene_id;
					my $excel_row=$excel_line_index+1;
					push @output, "=COUNTIFS(A2:A$range, A$excel_row, B2:B$range, B$excel_row, D2:D$range, 1)";
					push @output, "=SUBTOTAL(3, B$excel_row)";
					foreach my $sample (@sample_names) {
						push @output, "=IF(COUNTIFS(A2:A$range, \"=$sample\", D2:D$range, 1, B2:B$range, B$excel_row, C2:C$range, \">=2\")>=1, \"T\", \"F\")";
					}
					#push formula into
					#push @output;
					my @row = split("\t", $var_string);
					foreach my $var_ele (@row){
						#print "here is the ele: $var_ele\n";
						#print "superDupsRow: $superDupsRow\n";
						push @output, $var_ele;
					}
					if ($row[$superDupsRow] ne "NA") {
						$worksheet_compound->write_row($excel_line_index,0,\@output, $format_segmental_dup);
						$excel_line_index += 1;
					} else {
						$worksheet_compound->write_row($excel_line_index,0,\@output);
						$excel_line_index += 1;
					}
				}
			}
			if( exists $CNV_on_samples{$sample}{$gene_id}) {
				my @gene_sample_cnv = split(";", $CNV_on_samples{$sample}{$gene_id});
				foreach my $cnv_string (@gene_sample_cnv) {
					my @output;
					push @output, $sample;
					push @output, $gene_id;
					my $excel_row=$excel_line_index+1;
					push @output, "=COUNTIFS(A2:A$range, A$excel_row, B2:B$range, B$excel_row, D2:D$range, 1)";
					push @output, "=SUBTOTAL(3, B$excel_row)";
					foreach my $sample (@sample_names) {
						push @output, "=IF(COUNTIFS(A2:A$range, \"=$sample\", D2:D$range, 1, B2:B$range, B$excel_row, C2:C$range, \">=2\")>=1, \"T\", \"F\")";
					}
					#push formula into
					#push @output;
					my @row = split("\t", $cnv_string);
					my @extra_gene_info=split("\t", $all_gene_id{$gene_id});
					my $is_interested_index;
					for(my $j=$start_output_index+1; $j<scalar @compound_header; $j++) {
						if ($compound_header[$j] eq "Chr") {
							push @output,"chr".$row[1];
							next;
						}
						if ($compound_header[$j] eq "Func.ensGene") {
							push @output,$row[0];
							next;
						}
						if ($compound_header[$j] eq "ExonicFunc.ensGene") {
							push @output,$row[2];
							next;
						}
						if ($compound_header[$j] eq "AAChange.ensGene") {
							push @output,$row[2];
							next;
						}
						if ($compound_header[$j] eq "isInterested") {
							push @output,$row[3];
							next;
						}
						if ($compound_header[$j] eq "GeneName") {
							push @output,$extra_gene_info[0];
							next;
						}
						#GoTerm	WikiGene_Description	MIM_Gene_Description	OMIM_Gene_Description	GeneCard Link	OMIM Link	Uniprot Link
						if ($compound_header[$j] eq "GoTerm") {
							push @output,$extra_gene_info[2];
							next;
						}
						if ($compound_header[$j] eq "WikiGene_Description") {
							push @output,$extra_gene_info[3];
							next;
						}
						if ($compound_header[$j] eq "MIM_Gene_Description") {
							push @output,$extra_gene_info[4];
							next;
						}
						if ($compound_header[$j] eq "OMIM_Gene_Description") {
							push @output,$extra_gene_info[5];
							next;
						}
						if ($compound_header[$j] eq "GeneCard Link") {
							push @output,$extra_gene_info[6];
							next;
						}
						if ($compound_header[$j] eq "OMIM Link") {
							push @output,$extra_gene_info[7];
							next;
						}
						if ($compound_header[$j] eq "Uniprot Link") {
							push @output,$extra_gene_info[8];
							next;
						}
						
						
						push @output,"0";
					}
					$worksheet_compound->write_row($excel_line_index,0,\@output);
					$excel_line_index += 1;
				}
			}
		}
	}
}
my @dummy = ('dummy_row');
$worksheet_compound->write_row($excel_line_index,0,\@dummy);
############################# output compound heterozygous sheet  end here ###################################################

$workbook->close();

exit;



#Subroutines


sub usage {
    print "Unknown option: @_\n" if ( @_ );
    print "\nusage: VCF_2_annotated_excel_20131120.pl \n";
    print "--vcf input vcf file (of a single sample or a family;\n";
    print "--InterestedGenes file of interested gene names list (optional); Format: Ensembl gene IDs on the 1st column.\n"; 
    print "--out output excel file name (optional, only when required);\n";
    print "--outAll output \'Everthing\' table to a sparate file to reduce memory loading excel file (optional, only when required).\n";
    print "--CNV CNV result file, output of HG's Annotate_CNVs_combine_multiple_files.pl, sample names must be consistent with the vcf.\n";
    print "--add_genotypeCall_flags=(Yes/No), anyt word starting with letter \"y\" will be considered as yes.\n";
    print "--outVCF filename of the output VCF file that has extra infor in the INFO column.\n";
    print "--AnnovarDIR installation directory of annovar. \n\n";
	return(1);
}

sub VCFHeader {
	my $header_string=
'##fileformat=VCFv4.1
##FILTER=<ID=LowQual,Description="Low quality">
##FILTER=<ID=VQSRTrancheINDEL99.00to99.90,Description="Truth sensitivity tranche level for INDEL model at VQS Lod: -4.6047 <= x < 0.4197">
##FILTER=<ID=VQSRTrancheINDEL99.90to100.00+,Description="Truth sensitivity tranche level for INDEL model at VQS Lod < -28779.2906">
##FILTER=<ID=VQSRTrancheINDEL99.90to100.00,Description="Truth sensitivity tranche level for INDEL model at VQS Lod: -28779.2906 <= x < -4.6047">
##FILTER=<ID=VQSRTrancheSNP99.90to100.00+,Description="Truth sensitivity tranche level for SNP model at VQS Lod < -2742.7724">
##FILTER=<ID=VQSRTrancheSNP99.90to100.00,Description="Truth sensitivity tranche level for SNP model at VQS Lod: -2742.7724 <= x < -0.8995">
##FORMAT=<ID=AD,Number=.,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth (reads with MQ=255 or with bad mates are filtered)">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=MIN_DP,Number=1,Type=Integer,Description="Minimum DP observed within the GVCF block">
##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Normalized, Phred-scaled likelihoods for genotypes as defined in the VCF specification">
##FORMAT=<ID=SB,Number=4,Type=Integer,Description="Per-sample component statistics which comprise the Fisher\'s Exact Test to detect strand bias.">
##FORMAT=<ID=CV,Number=.,Type=String,Description="define which caller called the variant H: HaplotypeCaller F: Freebayes">
##FORMAT=<ID=AO,Number=A,Type=Integer,Description="Alternate allele observation count">
##FORMAT=<ID=GL,Number=G,Type=Float,Description="Genotype Likelihood, log10-scaled likelihoods of the data given the called genotype for each possible genotype generated from the reference and alternate alleles given the sample ploidy">
##FORMAT=<ID=QA,Number=A,Type=Integer,Description="Sum of quality of the alternate observations">
##FORMAT=<ID=QR,Number=1,Type=Integer,Description="Sum of quality of the reference observations">
##FORMAT=<ID=RO,Number=1,Type=Integer,Description="Reference allele observation count">
##INFO=<ID=AC,Number=A,Type=Integer,Description="Allele count in genotypes, for each ALT allele, in the same order as listed">
##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency, for each ALT allele, in the same order as listed">
##INFO=<ID=AN,Number=1,Type=Integer,Description="Total number of alleles in called genotypes">
##INFO=<ID=BaseQRankSum,Number=1,Type=Float,Description="Z-score from Wilcoxon rank sum test of Alt Vs. Ref base qualities">
##INFO=<ID=CCC,Number=1,Type=Integer,Description="Number of called chromosomes">
##INFO=<ID=ClippingRankSum,Number=1,Type=Float,Description="Z-score From Wilcoxon rank sum test of Alt vs. Ref number of hard clipped bases">
##INFO=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth; some reads may have been filtered">
##INFO=<ID=DS,Number=0,Type=Flag,Description="Were any of the samples downsampled?">
##INFO=<ID=END,Number=1,Type=Integer,Description="Stop position of the interval">
##INFO=<ID=FS,Number=1,Type=Float,Description="Phred-scaled p-value using Fisher\'s exact test to detect strand bias">
##INFO=<ID=GQ_MEAN,Number=1,Type=Float,Description="Mean of all GQ values">
##INFO=<ID=GQ_STDDEV,Number=1,Type=Float,Description="Standard deviation of all GQ values">
##INFO=<ID=HWP,Number=1,Type=Float,Description="P value from test of Hardy Weinberg Equilibrium">
##INFO=<ID=HaplotypeScore,Number=1,Type=Float,Description="Consistency of the site with at most two segregating haplotypes">
##INFO=<ID=InbreedingCoeff,Number=1,Type=Float,Description="Inbreeding coefficient as estimated from the genotype likelihoods per-sample when compared against the Hardy-Weinberg expectation">
##INFO=<ID=MLEAC,Number=A,Type=Integer,Description="Maximum likelihood expectation (MLE) for the allele counts (not necessarily the same as the AC), for each ALT allele, in the same order as listed">
##INFO=<ID=MLEAF,Number=A,Type=Float,Description="Maximum likelihood expectation (MLE) for the allele frequency (not necessarily the same as the AF), for each ALT allele, in the same order as listed">
##INFO=<ID=MQ,Number=1,Type=Float,Description="RMS Mapping Quality">
##INFO=<ID=MQ0,Number=1,Type=Integer,Description="Total Mapping Quality Zero Reads">
##INFO=<ID=MQRankSum,Number=1,Type=Float,Description="Z-score From Wilcoxon rank sum test of Alt vs. Ref read mapping qualities">
##INFO=<ID=NCC,Number=1,Type=Integer,Description="Number of no-called samples">
##INFO=<ID=NEGATIVE_TRAIN_SITE,Number=0,Type=Flag,Description="This variant was used to build the negative training set of bad variants">
##INFO=<ID=POSITIVE_TRAIN_SITE,Number=0,Type=Flag,Description="This variant was used to build the positive training set of good variants">
##INFO=<ID=QD,Number=1,Type=Float,Description="Variant Confidence/Quality by Depth">
##INFO=<ID=ReadPosRankSum,Number=1,Type=Float,Description="Z-score from Wilcoxon rank sum test of Alt vs. Ref read position bias">
##INFO=<ID=VQSLOD,Number=1,Type=Float,Description="Log odds ratio of being a true variant versus being false under the trained gaussian mixture model">
##INFO=<ID=culprit,Number=1,Type=String,Description="The annotation which was the worst performing in the Gaussian mixture model, likely the reason why the variant was filtered out">
##INFO=<ID=set,Number=1,Type=String,Description="Source VCF for the merged record in CombineVariants">
##INFO=<ID=AB,Number=A,Type=Float,Description="Allele balance at heterozygous sites: a number between 0 and 1 representing the ratio of reads showing the reference allele to all reads, considering only reads from individuals called as heterozygous">
##INFO=<ID=ABP,Number=A,Type=Float,Description="Allele balance probability at heterozygous sites: Phred-scaled upper-bounds estimate of the probability of observing the deviation between ABR and ABA given E(ABR/ABA) ~ 0.5, derived using Hoeffding\'s inequality">
##INFO=<ID=AO,Number=A,Type=Integer,Description="Alternate allele observations, with partial observations recorded fractionally">
##INFO=<ID=CIGAR,Number=A,Type=String,Description="The extended CIGAR representation of each alternate allele, with the exception that \'=\' is replaced by \'M\' to ease VCF parsing.  Note that INDEL alleles do not have the first matched base (which is provided by default, per the spec) referred to by the CIGAR.">
##INFO=<ID=DPB,Number=1,Type=Float,Description="Total read depth per bp at the locus; bases in reads overlapping / bases in haplotype">
##INFO=<ID=DPRA,Number=A,Type=Float,Description="Alternate allele depth ratio.  Ratio between depth in samples with each called alternate allele and those without.">
##INFO=<ID=EPP,Number=A,Type=Float,Description="End Placement Probability: Phred-scaled upper-bounds estimate of the probability of observing the deviation between EL and ER given E(EL/ER) ~ 0.5, derived using Hoeffding\'s inequality">
##INFO=<ID=EPPR,Number=1,Type=Float,Description="End Placement Probability for reference observations: Phred-scaled upper-bounds estimate of the probability of observing the deviation between EL and ER given E(EL/ER) ~ 0.5, derived using Hoeffding\'s inequality">
##INFO=<ID=GTI,Number=1,Type=Integer,Description="Number of genotyping iterations required to reach convergence or bailout.">
##INFO=<ID=LEN,Number=A,Type=Integer,Description="allele length">
##INFO=<ID=MEANALT,Number=A,Type=Float,Description="Mean number of unique non-reference allele observations per sample with the corresponding alternate alleles.">
##INFO=<ID=MQM,Number=A,Type=Float,Description="Mean mapping quality of observed alternate alleles">
##INFO=<ID=MQMR,Number=1,Type=Float,Description="Mean mapping quality of observed reference alleles">
##INFO=<ID=NS,Number=1,Type=Integer,Description="Number of samples with data">
##INFO=<ID=NUMALT,Number=1,Type=Integer,Description="Number of unique non-reference alleles in called genotypes at this position.">
##INFO=<ID=ODDS,Number=1,Type=Float,Description="The log odds ratio of the best genotype combination to the second-best.">
##INFO=<ID=PAIRED,Number=A,Type=Float,Description="Proportion of observed alternate alleles which are supported by properly paired read fragments">
##INFO=<ID=PAIREDR,Number=1,Type=Float,Description="Proportion of observed reference alleles which are supported by properly paired read fragments">
##INFO=<ID=PAO,Number=A,Type=Float,Description="Alternate allele observations, with partial observations recorded fractionally">
##INFO=<ID=PQA,Number=A,Type=Float,Description="Alternate allele quality sum in phred for partial observations">
##INFO=<ID=PQR,Number=1,Type=Float,Description="Reference allele quality sum in phred for partial observations">
##INFO=<ID=PRO,Number=1,Type=Float,Description="Reference allele observation count, with partial observations recorded fractionally">
##INFO=<ID=QA,Number=A,Type=Integer,Description="Alternate allele quality sum in phred">
##INFO=<ID=QR,Number=1,Type=Integer,Description="Reference allele quality sum in phred">
##INFO=<ID=RO,Number=1,Type=Integer,Description="Reference allele observation count, with partial observations recorded fractionally">
##INFO=<ID=RPP,Number=A,Type=Float,Description="Read Placement Probability: Phred-scaled upper-bounds estimate of the probability of observing the deviation between RPL and RPR given E(RPL/RPR) ~ 0.5, derived using Hoeffding\'s inequality">
##INFO=<ID=RPPR,Number=1,Type=Float,Description="Read Placement Probability for reference observations: Phred-scaled upper-bounds estimate of the probability of observing the deviation between RPL and RPR given E(RPL/RPR) ~ 0.5, derived using Hoeffding\'s inequality">
##INFO=<ID=RUN,Number=A,Type=Integer,Description="Run length: the number of consecutive repeats of the alternate allele in the reference genome">
##INFO=<ID=SAF,Number=A,Type=Integer,Description="Number of alternate observations on the forward strand">
##INFO=<ID=SAP,Number=A,Type=Float,Description="Strand balance probability for the alternate allele: Phred-scaled upper-bounds estimate of the probability of observing the deviation between SAF and SAR given E(SAF/SAR) ~ 0.5, derived using Hoeffding\'s inequality">
##INFO=<ID=SAR,Number=A,Type=Integer,Description="Number of alternate observations on the reverse strand">
##INFO=<ID=SRF,Number=1,Type=Integer,Description="Number of reference observations on the forward strand">
##INFO=<ID=SRP,Number=1,Type=Float,Description="Strand balance probability for the reference allele: Phred-scaled upper-bounds estimate of the probability of observing the deviation between SRF and SRR given E(SRF/SRR) ~ 0.5, derived using Hoeffding\'s inequality">
##INFO=<ID=SRR,Number=1,Type=Integer,Description="Number of reference observations on the reverse strand">
##INFO=<ID=TYPE,Number=A,Type=String,Description="The type of allele, either snp, mnp, ins, del, or complex.">
##INFO=<ID=technology.ILLUMINA,Number=A,Type=Float,Description="Fraction of observations supporting the alternate observed in reads from ILLUMINA">
##INFO=<ID=ANNO_INT,Number=A,Type=String,Description="If the variant is on a interested gene">
##INFO=<ID=ANNO_POPMAX,Number=1,Type=Float,Description="Max MAF of the FIRST Alternative allele in other popular DBs, including 1000G, ESP6500,cg69.">
##INFO=<ID=ANNO_CADD,Number=1,Type=Float,Description="CADD score of the position">
##INFO=<ID=ANNO_B_MAF,Number=1,Type=Float,Description="Batch MAF of the FIRST alternative allele">
##INFO=<ID=ANNO_InHo_MAF,Number=1,Type=Float,Description="In house MAF of the FIRST alternative allele">
##INFO=<ID=ANNO_InHoG_MAF,Number=1,Type=Float,Description="In house GATK MAF of the FIRST alternative allele">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT';
	return $header_string;
}

sub CNV_VCFHeader {
	my $header_string=
'##fileformat=VCFv4.1
##ALT=<ID=CNV,Description="Copy number variable region">
##FILTER=<ID=LowQual,Description="Low quality">
##FILTER=<ID=VQSRTrancheINDEL99.00to99.90,Description="Truth sensitivity tranche level for INDEL model at VQS Lod: -4.6047 <= x < 0.4197">
##FILTER=<ID=VQSRTrancheINDEL99.90to100.00+,Description="Truth sensitivity tranche level for INDEL model at VQS Lod < -28779.2906">
##FILTER=<ID=VQSRTrancheINDEL99.90to100.00,Description="Truth sensitivity tranche level for INDEL model at VQS Lod: -28779.2906 <= x < -4.6047">
##FILTER=<ID=VQSRTrancheSNP99.90to100.00+,Description="Truth sensitivity tranche level for SNP model at VQS Lod < -2742.7724">
##FILTER=<ID=VQSRTrancheSNP99.90to100.00,Description="Truth sensitivity tranche level for SNP model at VQS Lod: -2742.7724 <= x < -0.8995">
##FORMAT=<ID=AD,Number=.,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth (reads with MQ=255 or with bad mates are filtered)">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=MIN_DP,Number=1,Type=Integer,Description="Minimum DP observed within the GVCF block">
##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Normalized, Phred-scaled likelihoods for genotypes as defined in the VCF specification">
##FORMAT=<ID=SB,Number=4,Type=Integer,Description="Per-sample component statistics which comprise the Fisher\'s Exact Test to detect strand bias.">
##FORMAT=<ID=AO,Number=A,Type=Integer,Description="Alternate allele observation count">
##FORMAT=<ID=GL,Number=G,Type=Float,Description="Genotype Likelihood, log10-scaled likelihoods of the data given the called genotype for each possible genotype generated from the reference and alternate alleles given the sample ploidy">
##FORMAT=<ID=QA,Number=A,Type=Integer,Description="Sum of quality of the alternate observations">
##FORMAT=<ID=QR,Number=1,Type=Integer,Description="Sum of quality of the reference observations">
##FORMAT=<ID=RO,Number=1,Type=Integer,Description="Reference allele observation count">
##FORMAT=<ID=CN,Number=A,Type=Integer,Description="Copy number, 1 indicates deletion, 2 indicates duplication">
##FORMAT=<ID=CV,Number=.,Type=String,Description="define which caller called the variant H: HaplotypeCaller F: Freebayes">
##INFO=<ID=AC,Number=A,Type=Integer,Description="Allele count in genotypes, for each ALT allele, in the same order as listed">
##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency, for each ALT allele, in the same order as listed">
##INFO=<ID=AN,Number=1,Type=Integer,Description="Total number of alleles in called genotypes">
##INFO=<ID=BaseQRankSum,Number=1,Type=Float,Description="Z-score from Wilcoxon rank sum test of Alt Vs. Ref base qualities">
##INFO=<ID=CCC,Number=1,Type=Integer,Description="Number of called chromosomes">
##INFO=<ID=ClippingRankSum,Number=1,Type=Float,Description="Z-score From Wilcoxon rank sum test of Alt vs. Ref number of hard clipped bases">
##INFO=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth; some reads may have been filtered">
##INFO=<ID=DS,Number=0,Type=Flag,Description="Were any of the samples downsampled?">
##INFO=<ID=END,Number=1,Type=Integer,Description="Stop position of the interval">
##INFO=<ID=FS,Number=1,Type=Float,Description="Phred-scaled p-value using Fisher\'s exact test to detect strand bias">
##INFO=<ID=GQ_MEAN,Number=1,Type=Float,Description="Mean of all GQ values">
##INFO=<ID=GQ_STDDEV,Number=1,Type=Float,Description="Standard deviation of all GQ values">
##INFO=<ID=HWP,Number=1,Type=Float,Description="P value from test of Hardy Weinberg Equilibrium">
##INFO=<ID=HaplotypeScore,Number=1,Type=Float,Description="Consistency of the site with at most two segregating haplotypes">
##INFO=<ID=InbreedingCoeff,Number=1,Type=Float,Description="Inbreeding coefficient as estimated from the genotype likelihoods per-sample when compared against the Hardy-Weinberg expectation">
##INFO=<ID=MLEAC,Number=A,Type=Integer,Description="Maximum likelihood expectation (MLE) for the allele counts (not necessarily the same as the AC), for each ALT allele, in the same order as listed">
##INFO=<ID=MLEAF,Number=A,Type=Float,Description="Maximum likelihood expectation (MLE) for the allele frequency (not necessarily the same as the AF), for each ALT allele, in the same order as listed">
##INFO=<ID=MQ,Number=1,Type=Float,Description="RMS Mapping Quality">
##INFO=<ID=MQ0,Number=1,Type=Integer,Description="Total Mapping Quality Zero Reads">
##INFO=<ID=MQRankSum,Number=1,Type=Float,Description="Z-score From Wilcoxon rank sum test of Alt vs. Ref read mapping qualities">
##INFO=<ID=NCC,Number=1,Type=Integer,Description="Number of no-called samples">
##INFO=<ID=NEGATIVE_TRAIN_SITE,Number=0,Type=Flag,Description="This variant was used to build the negative training set of bad variants">
##INFO=<ID=POSITIVE_TRAIN_SITE,Number=0,Type=Flag,Description="This variant was used to build the positive training set of good variants">
##INFO=<ID=QD,Number=1,Type=Float,Description="Variant Confidence/Quality by Depth">
##INFO=<ID=ReadPosRankSum,Number=1,Type=Float,Description="Z-score from Wilcoxon rank sum test of Alt vs. Ref read position bias">
##INFO=<ID=VQSLOD,Number=1,Type=Float,Description="Log odds ratio of being a true variant versus being false under the trained gaussian mixture model">
##INFO=<ID=culprit,Number=1,Type=String,Description="The annotation which was the worst performing in the Gaussian mixture model, likely the reason why the variant was filtered out">
##INFO=<ID=set,Number=1,Type=String,Description="Source VCF for the merged record in CombineVariants">
##INFO=<ID=AB,Number=A,Type=Float,Description="Allele balance at heterozygous sites: a number between 0 and 1 representing the ratio of reads showing the reference allele to all reads, considering only reads from individuals called as heterozygous">
##INFO=<ID=ABP,Number=A,Type=Float,Description="Allele balance probability at heterozygous sites: Phred-scaled upper-bounds estimate of the probability of observing the deviation between ABR and ABA given E(ABR/ABA) ~ 0.5, derived using Hoeffding\'s inequality">
##INFO=<ID=AO,Number=A,Type=Integer,Description="Alternate allele observations, with partial observations recorded fractionally">
##INFO=<ID=CIGAR,Number=A,Type=String,Description="The extended CIGAR representation of each alternate allele, with the exception that \'=\' is replaced by \'M\' to ease VCF parsing.  Note that INDEL alleles do not have the first matched base (which is provided by default, per the spec) referred to by the CIGAR.">
##INFO=<ID=DPB,Number=1,Type=Float,Description="Total read depth per bp at the locus; bases in reads overlapping / bases in haplotype">
##INFO=<ID=DPRA,Number=A,Type=Float,Description="Alternate allele depth ratio.  Ratio between depth in samples with each called alternate allele and those without.">
##INFO=<ID=EPP,Number=A,Type=Float,Description="End Placement Probability: Phred-scaled upper-bounds estimate of the probability of observing the deviation between EL and ER given E(EL/ER) ~ 0.5, derived using Hoeffding\'s inequality">
##INFO=<ID=EPPR,Number=1,Type=Float,Description="End Placement Probability for reference observations: Phred-scaled upper-bounds estimate of the probability of observing the deviation between EL and ER given E(EL/ER) ~ 0.5, derived using Hoeffding\'s inequality">
##INFO=<ID=GTI,Number=1,Type=Integer,Description="Number of genotyping iterations required to reach convergence or bailout.">
##INFO=<ID=LEN,Number=A,Type=Integer,Description="allele length">
##INFO=<ID=MEANALT,Number=A,Type=Float,Description="Mean number of unique non-reference allele observations per sample with the corresponding alternate alleles.">
##INFO=<ID=MQM,Number=A,Type=Float,Description="Mean mapping quality of observed alternate alleles">
##INFO=<ID=MQMR,Number=1,Type=Float,Description="Mean mapping quality of observed reference alleles">
##INFO=<ID=NS,Number=1,Type=Integer,Description="Number of samples with data">
##INFO=<ID=NUMALT,Number=1,Type=Integer,Description="Number of unique non-reference alleles in called genotypes at this position.">
##INFO=<ID=ODDS,Number=1,Type=Float,Description="The log odds ratio of the best genotype combination to the second-best.">
##INFO=<ID=PAIRED,Number=A,Type=Float,Description="Proportion of observed alternate alleles which are supported by properly paired read fragments">
##INFO=<ID=PAIREDR,Number=1,Type=Float,Description="Proportion of observed reference alleles which are supported by properly paired read fragments">
##INFO=<ID=PAO,Number=A,Type=Float,Description="Alternate allele observations, with partial observations recorded fractionally">
##INFO=<ID=PQA,Number=A,Type=Float,Description="Alternate allele quality sum in phred for partial observations">
##INFO=<ID=PQR,Number=1,Type=Float,Description="Reference allele quality sum in phred for partial observations">
##INFO=<ID=PRO,Number=1,Type=Float,Description="Reference allele observation count, with partial observations recorded fractionally">
##INFO=<ID=QA,Number=A,Type=Integer,Description="Alternate allele quality sum in phred">
##INFO=<ID=QR,Number=1,Type=Integer,Description="Reference allele quality sum in phred">
##INFO=<ID=RO,Number=1,Type=Integer,Description="Reference allele observation count, with partial observations recorded fractionally">
##INFO=<ID=RPP,Number=A,Type=Float,Description="Read Placement Probability: Phred-scaled upper-bounds estimate of the probability of observing the deviation between RPL and RPR given E(RPL/RPR) ~ 0.5, derived using Hoeffding\'s inequality">
##INFO=<ID=RPPR,Number=1,Type=Float,Description="Read Placement Probability for reference observations: Phred-scaled upper-bounds estimate of the probability of observing the deviation between RPL and RPR given E(RPL/RPR) ~ 0.5, derived using Hoeffding\'s inequality">
##INFO=<ID=RUN,Number=A,Type=Integer,Description="Run length: the number of consecutive repeats of the alternate allele in the reference genome">
##INFO=<ID=SAF,Number=A,Type=Integer,Description="Number of alternate observations on the forward strand">
##INFO=<ID=SAP,Number=A,Type=Float,Description="Strand balance probability for the alternate allele: Phred-scaled upper-bounds estimate of the probability of observing the deviation between SAF and SAR given E(SAF/SAR) ~ 0.5, derived using Hoeffding\'s inequality">
##INFO=<ID=SAR,Number=A,Type=Integer,Description="Number of alternate observations on the reverse strand">
##INFO=<ID=SRF,Number=1,Type=Integer,Description="Number of reference observations on the forward strand">
##INFO=<ID=SRP,Number=1,Type=Float,Description="Strand balance probability for the reference allele: Phred-scaled upper-bounds estimate of the probability of observing the deviation between SRF and SRR given E(SRF/SRR) ~ 0.5, derived using Hoeffding\'s inequality">
##INFO=<ID=SRR,Number=1,Type=Integer,Description="Number of reference observations on the reverse strand">
##INFO=<ID=TYPE,Number=A,Type=String,Description="The type of allele, either snp, mnp, ins, del, or complex.">
##INFO=<ID=technology.ILLUMINA,Number=A,Type=Float,Description="Fraction of observations supporting the alternate observed in reads from ILLUMINA">
##INFO=<ID=ANNO_INT,Number=A,Type=String,Description="If the variant is on a interested gene">
##INFO=<ID=ANNO_POPMAX,Number=1,Type=Float,Description="Max MAF of the FIRST Alternative allele in other popular DBs, including 1000G, ESP6500,cg69.">
##INFO=<ID=ANNO_CADD,Number=1,Type=Float,Description="CADD score of the position">
##INFO=<ID=ANNO_B_MAF,Number=1,Type=Float,Description="Batch MAF of the FIRST alternative allele">
##INFO=<ID=ANNO_InHo_MAF,Number=1,Type=Float,Description="In house MAF of the FIRST alternative allele">
##INFO=<ID=ANNO_InHoG_MAF,Number=1,Type=Float,Description="In house GATK MAF of the FIRST alternative allele">
##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description="Imprecise structural variation">
##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant described in this record">
##INFO=<ID=SVLEN,Number=-1,Type=Integer,Description="Difference in length between REF and ALT alleles">
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">
##INFO=<ID=ANNO_CNV_BF,Number=1,Type=Float,Description="Bayes Factor, the higher the more confidence of the CNV call">
##INFO=<ID=ANNO_CNV_EXP_READS,Number=1,Type=Float,Description="Expected count of the CNV">
##INFO=<ID=ANNO_CNV_OBS_READS,Number=1,Type=Float,Description="Observed count of the CNV">
##INFO=<ID=ANNO_CNV_RATIO,Number=1,Type=Float,Description="Ratio between the observed count and the expected coun for the CNV">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT';
	return $header_string;
}

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

