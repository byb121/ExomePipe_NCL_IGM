#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;


# YX edit: Compound het V1/V2 is counted as homozygous of the first minor allele.
# Filter a need  to add a filter a 
# Read in a folder and take only vcf files
#
#
#define variables
my $vcf=""; 
my $DupSamples = "Yes"; # If there are duplicated samples included in VCFs, "Yes" mean the script will ignore those duplicated ones.

print "\n";
print "\n";
print "##################################################################################\n";
print "  YX edit: Compound het V1/V2 is counted as homozygous of the first minor allele, \n";
print "           but it does not affect MAFs only homo counts.                          \n";
print "##################################################################################\n";
print "\n";
print "\n";

my %var;
my %vCount;
my %SampleNames;

my $Results=GetOptions("vcf=s"=>\$vcf, "DupSamples=s"=>\$DupSamples);

#open each vcf file and add variants into %vars
my $exomecount = 0;
my $exomecount_this_file = 0;
my $file = $vcf;
open INPUT, $file or die "Cannot open $file\n";
my %skipCol;
my $lineCount=0;
print "Start to process file: $file\n";
line: foreach my $Line (<INPUT>){
	$lineCount++;
	if($Line =~m/^#CHR/){
		chomp $Line;
		my @linesplit = split (/\t/,$Line);
		for(my $i=9;$i<scalar(@linesplit);$i++) {
			if (! exists $SampleNames{$linesplit[$i]} ) {
				 $SampleNames{$linesplit[$i]} = 1;
			} else {
				$skipCol{$i} = 1;
			}
		}
		next line;
	} 
	
	if($Line=~/^#/){next line;} 
	
	chomp $Line;
	my @linesplit = split (/\t/,$Line);
	my $chr = $linesplit[0];
	$chr =~ s/chr//;
	my $pos_1 = $linesplit[1];
	my $pos_2 = $pos_1;
	my $ref = $linesplit[3];
	my $var = $linesplit[4];
	my $format = $linesplit[8];
	
	if ($linesplit[6] !~ m/PASS/) { next; } #vcf filter
	
	my $Sample_Call="";
	sample: for(my $i=9;$i<scalar(@linesplit);$i++) {
		if ( $DupSamples =~ m/^y/i ) {			
			if (exists $skipCol{$i} ) {
				next sample;
			} else {
				$Sample_Call=$Sample_Call."\t".$linesplit[$i];
			}
		} else {
			$Sample_Call=$Sample_Call."\t".$linesplit[$i];
		}
	}
	$Sample_Call =~ s/^\t//;
	
	my @vars=();
	if($var!~/,/){push(@vars,$var);}
	if($var=~/,/){
		@vars=split(/,/,$var);
	}
	
	if($ref=~/,/){print "Warnning: \n$Line\n";} #check for multiple ref alleles!
	foreach my $v (@vars){
		my $ref_now = $ref;
		my $v_now = $v;
		my $pos_1_now = $pos_1;
		my $pos_2_now = $pos_2;
		my $Sample_Call_Processed = AddGenoTypeToSampleCalls_CompondHet($v, $var, $format, $Sample_Call);
		my $out_var=$v;
		for (my $i=0;$i<scalar @vars;$i++) {
			if($v ne $vars[$i]) {
				$out_var=$out_var.",".$vars[$i];
			}
		}
		
		if ($ref_now =~ m/^[atgcATGC]$/) { # singlge base ref
			#insertion
			if(length($v_now)>length($ref_now)){
				#my @ref_elements=split("", $ref_now);
				my @var_chars=split(//, $v_now);
				if ($ref_now eq $var_chars[0]) {
					$ref_now="-";$v_now=~s/^.//;$pos_1_now++;$pos_2_now=$pos_1_now;
				}
			} #remove first 'ref' base if it is the same as first base of the variant
		} else { # multiple bases ref
			my @var_chars=split(//, $v_now);
			my @ref_chars=split(//, $ref_now);
			if(length($v_now) >= length($ref_now)){ #insertion
				Compare1:	for(my $i=0;$i<scalar @ref_chars;$i++) {
					if($ref_chars[$i] eq $var_chars[$i]) {
						if ($i != scalar @ref_chars-1) {
							my @temp_chars = split(//, $ref_now);
							shift @temp_chars;
							$ref_now=join '', @temp_chars;
						} else {
							$ref_now="-";
						}
						my @temp_chars = split(//, $v_now);
						shift @temp_chars;
						$v_now=join '', @temp_chars;
						$pos_1_now++;
					} else {
						last Compare1;
					}
				}
				$pos_2_now = $pos_1_now + length($ref_now) - 1;
			} else { #deletion
				Compare2:	for(my $i=0;$i<scalar @var_chars;$i++) {
					if($var_chars[$i] eq $ref_chars[$i]) {
						if ($i != scalar @var_chars -1) {
							my @temp_chars = split(//, $v_now);
							shift @temp_chars;
							$v_now=join '', @temp_chars;
						} else {
							$v_now="-";
						}
						my @temp_chars = split(//, $ref_now);
						shift @temp_chars;
						$ref_now=join '', @temp_chars;
						$pos_1_now++;
					} else {
						last Compare2;
					}
				}
				$pos_2_now = $pos_1_now + length($ref_now) - 1;
			}
		}
		
		#print OUT "$chr\t$pos_1_now\t$pos_2_now\t$ref_now\t$v_now\t$linesplit[0]\t$linesplit[1]\t$linesplit[2]\t".
           #"$linesplit[3]\t$out_var\t$linesplit[5]\t$linesplit[6]\t$linesplit[7]\t$linesplit[8]\t$Sample_Call_Processed\n";
		my @samplesCalls = split("\t", $Sample_Call_Processed);
		my @temp = split("\t", $Sample_Call);
		if ($exomecount_this_file == 0) {$exomecount_this_file = scalar(@temp);}
		
		for(my $i=1;$i<=scalar(@samplesCalls);$i=$i+2){# checking on samples
			if ($samplesCalls[$i] =~ /R\/R/ || $samplesCalls[$i] =~ /NA/) {
				next;
			} elsif ($samplesCalls[$i] =~ /R\/V/) {
				if ($samplesCalls[$i-1] =~ /0\/1/) {
					if (! exists $var{$chr}{$pos_1_now}{$ref_now}{$v_now}{'hets'}) { $var{$chr}{$pos_1_now}{$ref_now}{$v_now}{'hets'} = 1;}
					else {$var{$chr}{$pos_1_now}{$ref_now}{$v_now}{'hets'} = $var{$chr}{$pos_1_now}{$ref_now}{$v_now}{'hets'} + 1;}
				}
				if ($samplesCalls[$i-1] =~ /0\/1/) {
					if (! exists $vCount{$chr}{$pos_1_now}{$ref_now}{$v_now}) { $vCount{$chr}{$pos_1_now}{$ref_now}{$v_now} = 1;}
					else {$vCount{$chr}{$pos_1_now}{$ref_now}{$v_now} = $vCount{$chr}{$pos_1_now}{$ref_now}{$v_now} + 1;}
				}
			} elsif ($samplesCalls[$i] =~ /V\/V/) {
				if ($samplesCalls[$i-1] =~ /1\/1/) {
					if (! exists $var{$chr}{$pos_1_now}{$ref_now}{$v_now}{'homo'}) { $var{$chr}{$pos_1_now}{$ref_now}{$v_now}{'homo'} = 1;}
					else {$var{$chr}{$pos_1_now}{$ref_now}{$v_now}{'homo'} = $var{$chr}{$pos_1_now}{$ref_now}{$v_now}{'homo'} + 1;}
				}
				if ($samplesCalls[$i-1] =~ /1\/1/) {
					if (! exists $vCount{$chr}{$pos_1_now}{$ref_now}{$v_now}) { $vCount{$chr}{$pos_1_now}{$ref_now}{$v_now} = 2;}
					else {$vCount{$chr}{$pos_1_now}{$ref_now}{$v_now} = $vCount{$chr}{$pos_1_now}{$ref_now}{$v_now} + 2;}
				}
			} elsif ($samplesCalls[$i] =~ /V1\/V2/) {
				if ($samplesCalls[$i-1] =~ /1\/\d/) {
					if (! exists $var{$chr}{$pos_1_now}{$ref_now}{$v_now}{'homo'}) { $var{$chr}{$pos_1_now}{$ref_now}{$v_now}{'homo'} = 1;}
					else {$var{$chr}{$pos_1_now}{$ref_now}{$v_now}{'homo'} = $var{$chr}{$pos_1_now}{$ref_now}{$v_now}{'homo'} + 1;}
				}
				if ($samplesCalls[$i-1] =~ /1\/\d/) {
					if (! exists $vCount{$chr}{$pos_1_now}{$ref_now}{$v_now}) { $vCount{$chr}{$pos_1_now}{$ref_now}{$v_now} = 1;}
					else {$vCount{$chr}{$pos_1_now}{$ref_now}{$v_now} = $vCount{$chr}{$pos_1_now}{$ref_now}{$v_now} + 1;}
				}
				
			}
		}
	}
	if ($lineCount % 10000 == 0) {
		print "Processed $lineCount lines.\n";
	}
}	           
close INPUT;
$exomecount = $exomecount + $exomecount_this_file;
$exomecount_this_file = 0;


my $output_file = $vcf."_MAFs_".$exomecount."exomes.txt";
open(OUT, ">$output_file") or die "Cannot open file \"$output_file\" to write to!\n";
foreach my $c (sort keys %var){
	foreach my $p (sort {$a<=>$b} keys %{$var{$c}}){
		foreach my $r (sort keys %{$var{$c}{$p}}){
			foreach my $v (sort keys %{$var{$c}{$p}{$r}}){
			my $maf=$vCount{$c}{$p}{$r}{$v}/($exomecount*2);
				if (! exists $var{$c}{$p}{$r}{$v}{'hets'}) {$var{$c}{$p}{$r}{$v}{'hets'} = 0;}
				if (! exists $var{$c}{$p}{$r}{$v}{'homo'}) {$var{$c}{$p}{$r}{$v}{'homo'} = 0;}
				my $total = $var{$c}{$p}{$r}{$v}{'hets'} + $var{$c}{$p}{$r}{$v}{'homo'};	
				print OUT "$c\t$p\t$r\t$v\t$total\t$var{$c}{$p}{$r}{$v}{'hets'}\t$var{$c}{$p}{$r}{$v}{'homo'}\t$maf\n";
			}
		}
	}
}
close OUT;

exit;

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



