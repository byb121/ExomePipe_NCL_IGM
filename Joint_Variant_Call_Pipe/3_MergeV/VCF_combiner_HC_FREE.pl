#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;

print "\n";
print "####################################################################################\n";
print "#         output a new vcf file which have both HC and Freebayes calls             #\n";
print "#  insert a new format field for each sample call -  CV (Called in variant caller) #\n";
print "# H: Called in HC only, HF: Called in both callers, F: Called in  Freebayes only.  #\n";
print "#              Consensus rules may be vary from version to version                 #\n";
print "#                                Check with an expert for details                  #\n";
print "####################################################################################\n";


my $fai_file; #this is used to retreive chromosome names
my $input_FREE_vcf;
my $input_HC_vcf;
my $output_vcf;
my $output_MAF;

my $help;
usage() if ( @ARGV < 1 || ! GetOptions('help|?' => \$help, "fai=s"=>\$fai_file, "FREE=s"=>\$input_FREE_vcf, 
	'HC=s' => \$input_HC_vcf, 'output=s' => \$output_vcf, 'MAF=s' => \$output_MAF ) || defined $help );

#read in the fai file
my %chros;
my @chrosomes;
open FAI, "$fai_file" or die "Can not open the fasta index file: $fai_file\n";
while (my $line = <FAI> ) {
	chomp $line;
	my @eles = split(/\s+/, $line);
	if (! exists $chros{$eles[0]}) {
		$chros{$eles[0]} = $eles[1];
		push @chrosomes, $eles[0];
	} else {
		print "Error: duplicated chrosome name found.\n";
		exit 1;
	}
}

# need several line to play with the headers
open OUTPUT, ">$output_vcf" or die "Can not output file: $output_vcf\n";
print "Write vcf header to output\n";
my $header = Define_Header($input_FREE_vcf, $input_HC_vcf);
print OUTPUT $header;

foreach my $chr (@chrosomes) {
	#open both inputs and output
	#my $chr="MT"; # to test, will need to loop all chromosomes
	print "Start to processing files for $chr..\n";
	open my $FREE, "$input_FREE_vcf" or die "Can not reads the freebayes vcf file: $input_FREE_vcf\n";
	open my $HC, "$input_HC_vcf" or die "Can not open the HC vcf file: $input_HC_vcf\n";
	
	my $free_line = read_file_line($FREE);
	my $HC_line = read_file_line($HC);
	#print "wo zai zhe ------------------------ 1.\n";
	while ( $free_line or $HC_line) {
		if ($free_line and $HC_line) {
			#print "wo zai zhe ------------------------ 2.\n";
			my $free_chr;
			my $free_pos;
			my $HC_chr;
			my $HC_pos;
			if($free_line =~ m/^$chr\t/ && $HC_line =~ m/^$chr\t/) {
				#print "wo zai zhe ------------------------ 3.\n";
				my @f_eles = split("\t", $free_line);
				$free_chr = $f_eles[0];
				$free_pos = $f_eles[1];
				my @h_eles = split("\t", $HC_line);
				$HC_chr = $h_eles[0];
				$HC_pos = $h_eles[1];
				### has to be a keep going loop
				if ($free_pos < $HC_pos) {
					#print "wo zai zhe ------------------------ 4.\n";
					my $output_line = FREE_filter($free_line);
					if ($output_line) {
								#print OUTPUT "free:".$free_line."\n";
								#print OUTPUT "HC:".$HC_line."\n";
								#print OUTPUT "next free:".$output_line."\n";
						print OUTPUT $output_line."\n";
					}
					$free_line = read_file_line($FREE);
				} elsif ($free_pos > $HC_pos) {
					#print "wo zai zhe ------------------------ 5.\n";
					my $output_line = HC_filter($HC_line);
					if ($output_line) {
								#print OUTPUT "free:".$free_line."\n";
								#print OUTPUT "HC:".$HC_line."\n";
								#print OUTPUT "next HC:".$output_line."\n";
						print OUTPUT $output_line."\n";
					}
					$HC_line = read_file_line($HC);
				} else {
					#print "wo zai zhe ------------------------ 6.\n";
					my $output_line = Consensus_call($HC_line, $free_line);
					if ($output_line) {
								#print OUTPUT "free:".$free_line."\n";
								#print OUTPUT "HC:".$HC_line."\n";
								#print OUTPUT "Consensus:".$output_line."\n";
						print OUTPUT $output_line."\n";		
					}
					$HC_line = read_file_line($HC);
					$free_line = read_file_line($FREE);
				}
			} elsif ($free_line =~ m/^$chr\t/) {
				#print "wo zai zhe ------------------------- 7.\n";
				#print "HC line is not started with $chr: $HC_line\n";
				$HC_line = read_file_line($HC);
			} elsif ($HC_line =~ m/^$chr\t/) {
				#print "wo zai zhe ------------------------ 8.\n";
				#print "Free line is not started with $chr: $free_line\n";
				$free_line = read_file_line($FREE);
			} else {
				$HC_line = read_file_line($HC);
				$free_line = read_file_line($FREE);
			}
		} elsif ($free_line) {
			#print "wo zai zhe ------------------------ 9.\n";
			#print "HC line in the end: $HC_line\n";
			if($free_line =~ m/^$chr\t/) {
				#print "wo zai zhe ------------------------ 10.\n";
				my $output_line = FREE_filter($free_line);
				if ($output_line) {
								#print OUTPUT "free:".$free_line."\n";
								#print OUTPUT "free no HC:".$output_line."\n";
					print OUTPUT $output_line."\n";
				}
			}
			$free_line = read_file_line($FREE);
		} else {
			#print "wo zai zhe ------------------------ 11.\n";
			#print "free line in the end: $free_line\n";
			if ($HC_line =~ m/^$chr\t/) {
				#print "wo zai zhe ------------------------ 12.\n";
				my $output_line = HC_filter($HC_line);
				if ($output_line) {
								#print OUTPUT "HC:".$HC_line."\n";
								#print OUTPUT "HC no free:".$output_line."\n";
					print OUTPUT $output_line."\n";
				}
			}
			$HC_line = read_file_line($HC);
		}		
	}
	
	close $FREE;
	close $HC;
}
close OUTPUT;

print "Start to calculate batch MAFs from file $output_vcf\n";
BatchMAF($output_vcf, $output_MAF); # in order to remove common variants of this batch

print "All Done!\n";

exit;
#compare lines of the same chromosome

sub FREE_filter {
	my $line=shift;
	my $QUAL_threshold = 20;
	my $Depth_threshold = 5;
	my $AA_sup_count_threshold = 5;
	#my $GQ_threshold = 20;
	
	my @temp = split("\t", $line);
	if ($temp[5] >= $QUAL_threshold) {
		my $format = $temp[8];
		my @format_fields = split(":", $format);
		my $GT_index;
		my $GQ_index;
		my $DP_index;
		my $RO_index;
		my $AO_index;
		for (my $i=0;$i<scalar(@format_fields);$i++){
			if ($format_fields[$i] =~ m/GT/) {
				$GT_index = $i;
				#print $line."$GT_index\n";
			}
			if ($format_fields[$i] =~ m/GQ/) {
				$GQ_index = $i;
				#print $line."$GQ_index\n";
			}
			if ($format_fields[$i] =~ m/DP/) {
				$DP_index = $i;
				#print $line."$DP_index\n";
			}
			if ($format_fields[$i] =~ m/RO/) {
				$RO_index = $i;
				#print $line."$RO_index\n";
			}
			if ($format_fields[$i] =~ m/AO/) {
				$AO_index = $i;
				#print $line."$AO_index\n";
			}
		}
		
		if (defined $GT_index && defined $GQ_index && defined $DP_index && defined $RO_index
			&& defined $AO_index ) {
			my $output_string="";
			my $NA_count = 0;
			for(my $i=0;$i < scalar @temp; $i++){
				if ($i == 6) {
					$output_string = $output_string."PASS"."\t";
				} elsif ($i <=7) {
					$output_string = $output_string.$temp[$i]."\t";
				} elsif($i == 8) {
					$output_string = $output_string."GT:GQ:DP:AD:CV"."\t"; #change format
				} else {
					my $sample = $temp[$i];
					#print $sample."\n";
					if ($sample !~ m/^\./) {
						my @fields = split(":", $sample);
						my $GT = $fields[$GT_index];
						my $GQ = $fields[$GQ_index];
						my $DP = $fields[$DP_index];
						my $RO = $fields[$RO_index];
						my $AO = $fields[$AO_index];
						my $AD = "$RO,$AO";
						
						#Allele depth coreponding to the genotype for filtering
						my @bases = split ("/", $GT);
						my @alt_alleles = split(",", $AO);
						if ($GT !~ m/0\/0/) {
							if ($DP >= $Depth_threshold && $alt_alleles[$bases[1]-1] >= $AA_sup_count_threshold) {
								$output_string = $output_string."$GT:$GQ:$DP:$AD:F"."\t";
							} else {
								$NA_count += 1;
								$output_string = $output_string."./."."\t";
							}
						} else {
							if ($DP >= $Depth_threshold) {
								$output_string = $output_string."$GT:$GQ:$DP:$AD:F"."\t";
							} else {
								$NA_count += 1;
								$output_string = $output_string."./."."\t";
							}
						}
					} else {
						$NA_count += 1;
						$output_string = $output_string."./."."\t";
					}
				}
			}
			$output_string =~ s/\t$//;
			
			if ($NA_count < (scalar @temp - 9)) {
				return $output_string;
			} else {
				return 0;
			}
			
		} else {
			print "Error: Freebayes output not expected format line:\n$line\n;";
			exit;
		}
	} else {
		return 0;
	}
}

sub HC_filter {
	my $line=shift;
	#my $GQ_threshold = 20;
	
	my @temp = split("\t", $line);
	if ($temp[6] eq "PASS") {
		my $format = $temp[8];
		my @format_fields = split(":", $format);
		my $GT_index;
		my $GQ_index;
		my $DP_index;
		my $AD_index;
		for (my $i=0;$i<scalar(@format_fields);$i++){
			if ($format_fields[$i] =~ m/GT/) {
				$GT_index = $i;
				#print $line."$GT_index\n";
			}
			if ($format_fields[$i] =~ m/GQ/) {
				$GQ_index = $i;
				#print $line."$GQ_index\n";
			}
			if ($format_fields[$i] =~ m/DP/) {
				$DP_index = $i;
				#print $line."$DP_index\n";
			}
			if ($format_fields[$i] =~ m/AD/) {
				$AD_index = $i;
				#print $line."$AD_index\n";
			}
		}
		
		if (defined $GT_index && defined $GQ_index && defined $DP_index && defined $AD_index ) {
			my $output_string="";
			my $bad_count = 0;
			for(my $i=0;$i < scalar @temp; $i++){
				if ($i <=7) {
					$output_string = $output_string.$temp[$i]."\t";
				} elsif($i == 8) {
					$output_string = $output_string."GT:GQ:DP:AD:CV"."\t"; #change format
				} else {
					my $sample = $temp[$i];
					#print $sample."\n";
					if ($sample !~ m/^\./) {
						my @fields = split(":", $sample);
						my $GT = $fields[$GT_index];
						my $GQ = $fields[$GQ_index];
						my $DP = $fields[$DP_index];
						my $AD = $fields[$AD_index];
						my @bases = split ("/", $GT);
						my @alt_alleles = split(",", $AD);
						if ($DP =~ /^\./ || $alt_alleles[$bases[1]] == 0) {
							$bad_count += 1;
							$output_string = $output_string."./."."\t";
						} else {
							$output_string = $output_string."$GT:$GQ:$DP:$AD:H"."\t";
						}
					} else {
						$bad_count += 1;
						$output_string = $output_string."./."."\t";
					}
				}
			}
			$output_string =~ s/\t$//;
			if ($bad_count < (scalar @temp - 9)) {
				return $output_string;
				print "Error: HaplotypeCaller produces no valid call on this line:\n$line\n;";
			} else {
				return 0;
			}
		} elsif (defined $GT_index && defined $GQ_index) {
			my $output_string="";
			my $bad_count = 0;
			for(my $i=0;$i < scalar @temp; $i++){
				if ($i <=7) {
					$output_string = $output_string.$temp[$i]."\t";
				} elsif($i == 8) {
					$output_string = $output_string."GT:GQ:CV"."\t"; #change format
				} else {
					my $sample = $temp[$i];
					#print $sample."\n";
					if ($sample !~ m/^\./) {
						my @fields = split(":", $sample);
						my $GT = $fields[$GT_index];
						my $GQ = $fields[$GQ_index];
						$output_string = $output_string."$GT:$GQ:H"."\t";
					} else {
						$bad_count += 1;
						$output_string = $output_string."./."."\t";
					}
				}
			}
			$output_string =~ s/\t$//;
			if ($bad_count < (scalar @temp - 9)) {
				return $output_string;
			} else {
				print "Error: HaplotypeCaller produces no valid call on this line:\n$line\n;";
				return 0;
			}
		}else {
			print "Error: HaplotypeCaller outputs not expected format line:\n$line\n;";
			return 0;
		}
	} else {
		return 0;
	}
}

sub Consensus_call {
	my ($HC_line, $free_line)=@_;
	my $free_QUAL_threshold = 20; #use to filter freebayes output
	my $Depth_threshold = 5;
	my $AA_sup_count_threshold = 5;
	my @free_temp = split("\t", $free_line);
	my @HC_temp = split("\t", $HC_line);
	
	my @HC_format_fields = split(":", $HC_temp[8]); # to find out if HC line has very format field
	my $HC_GT_index;
	my $HC_GQ_index;
	my $HC_DP_index;
	my $HC_AD_index;
	for (my $i=0;$i<scalar(@HC_format_fields);$i++){
		if ($HC_format_fields[$i] =~ m/GT/) {
			$HC_GT_index = $i;
			#print $line."$HC_GT_index\n";
		}
		if ($HC_format_fields[$i] =~ m/GQ/) {
			$HC_GQ_index = $i;
			#print $line."$HC_GQ_index\n";
		}
		if ($HC_format_fields[$i] =~ m/DP/) {
			$HC_DP_index = $i;
			#print $line."$HC_DP_index\n";
		}
		if ($HC_format_fields[$i] =~ m/AD/) {
			$HC_AD_index = $i;
			#print $line."$HC_AD_index\n";
		}
	}
	
	
	if ($free_temp[5] >= $free_QUAL_threshold) { #if freebayes line is statify the qulaity filtering
		if ($HC_temp[6] eq "PASS") {
			my $free_ref = $free_temp[3];
			my $HC_ref = $HC_temp[3];
			if ( $free_ref ne $HC_ref) { # if the reference base are different, then consider they are two different variants
				my $temp_output1 = HC_filter($HC_line);
				my $temp_output2 = FREE_filter($free_line);
				if ($temp_output1 && $temp_output2) {
					return $temp_output1."\n".$temp_output2;
				} elsif ($temp_output1) {
					return $temp_output1;
				} else {
					return $temp_output2;
				}
			} else { #same position, same reference
			
				if (defined $HC_GT_index && defined $HC_GQ_index && defined $HC_DP_index 
				&& defined $HC_AD_index ) {
					#record each variant
					my @HC_vars = split(",", $HC_temp[4]);
					my @free_vars = split(",", $free_temp[4]);
					my %vars;
					my $order = 0; # all alternative allele indexed from 1
					foreach my $v (@HC_vars) {
						if (! exists $vars{$v}) {
							$order += 1;
							$vars{$v} = $order;
						}
					}
					foreach my $v (@free_vars) {
						if (! exists $vars{$v}) {
							$order += 1;
							$vars{$v} = $order;
						}
					}
					
					#decide free vcf format index
					my @free_format_fields = split(":", $free_temp[8]);
					my $free_GT_index;
					my $free_GQ_index;
					my $free_DP_index;
					my $free_RO_index;
					my $free_AO_index;
					for (my $i=0;$i<scalar(@free_format_fields);$i++){
						if ($free_format_fields[$i] =~ m/GT/) {
							$free_GT_index = $i;
							#print $line."$free_GT_index\n";
						}
						if ($free_format_fields[$i] =~ m/GQ/) {
							$free_GQ_index = $i;
							#print $line."$free_GQ_index\n";
						}
						if ($free_format_fields[$i] =~ m/DP/) {
							$free_DP_index = $i;
							#print $line."$free_DP_index\n";
						}
						if ($free_format_fields[$i] =~ m/RO/) {
							$free_RO_index = $i;
							#print $line."$free_RO_index\n";
						}
						if ($free_format_fields[$i] =~ m/AO/) {
							$free_AO_index = $i;
							#print $line."$free_AO_index\n";
						}
					}
					
					#consensus call
					#current rules:
					# 1. genotype of variants over-ride genotype of reference
					# 2. if both called reference or (HC called ref and Freebayes has fewer alt evidence)
					#  then trust HC output
					# 3. take higher GQ, if same then take HC output (Freebayes tends to give higher scores even with
					#  only few reads support the varints, thus filtering before merge)
					my $output_string="";
					my $NA_count=0;
					for(my $i=0;$i < scalar @HC_temp; $i++){
						if ($i == 4) {
							my $var_string="";
							foreach my $key ( sort { $vars{$a} <=> $vars{$b} } keys (%vars) ) {
								$var_string=$var_string.$key.",";
							}
							$var_string =~ s/\,$//;
							$output_string = $output_string.$var_string."\t";
						} elsif($i == 6) {
							$output_string = $output_string."PASS"."\t";
						} elsif($i <=7) {
							#variant call quality is using HC caller results
							$output_string = $output_string.$HC_temp[$i]."\t";
						} elsif($i == 8) {
							$output_string = $output_string."GT:GQ:DP:AD:CV"."\t"; #change format
						} else {
							my $HC_sample = $HC_temp[$i];
							my $free_sample = $free_temp[$i];
							my @HC_fields = split(":", $HC_sample);
							my @free_fields = split(":", $free_sample);
							
							#determine freebayes variant call's alt allele count
													
							
							if ($HC_sample !~ m/^\./ && $free_sample !~ m/^\./) {
								my @free_bases = split ("/", $free_fields[$free_GT_index]);
								my @free_alt_alleles = split(",", $free_fields[$free_AO_index]);
								
								my ($GT, $GQ, $DP, $AD, $CV);
								$GT = $HC_fields[$HC_GT_index];
								my @bases = split ("/", $GT);
								my ($left, $right);
								if ($bases[0] == 0) { $left = 0;} else {
									$left = $vars{$HC_vars[$bases[0]-1]};
								}
								if ($bases[1] == 0) { $right = 0;} else {
									$right = $vars{$HC_vars[$bases[1]-1]};
								}
								my $GT_HC = $left."/".$right;
								
								$GT = $free_fields[$free_GT_index];
								@bases = split ("/", $GT);
								if ($bases[0] == 0) { $left = 0; } else {
									$left = $vars{$free_vars[$bases[0]-1]};
								} 
								if ($bases[1] == 0) { $right = 0; } else {
									$right = $vars{$free_vars[$bases[1]-1]};
								}
								my $GT_free = $left."/".$right;
								
								if ($HC_fields[$HC_GT_index] !~ m/0\/0/ && $free_fields[$free_GT_index] !~ m/0\/0/
									&& $free_alt_alleles[$free_bases[1]-1] >= $AA_sup_count_threshold
									&& $free_fields[$free_DP_index] >= $Depth_threshold) {
									if ($HC_fields[$HC_GQ_index] >= $free_fields[$free_GQ_index]) {
										if ($GT_HC eq $GT_free) {
											$CV = "HF";
										} else {
											$CV = "H";
										}
										$GQ = $HC_fields[$HC_GQ_index];
										$DP = $HC_fields[$HC_DP_index];
										#Allele depth is corrupted when merge, so when use HC GT, variants in free only
										# will be assigned as 0
										#Arranging AD to vars
										my %var_p_HC;
										my @AD_HC = split(",", $HC_fields[$HC_AD_index]);
										for (my $i=1;$i < scalar @AD_HC; $i++ ) {
											$var_p_HC{$HC_vars[$i-1]} = $AD_HC[$i];
										}
										$AD = $AD_HC[0].",";
										foreach my $key ( sort { $vars{$a} <=> $vars{$b} } keys (%vars) ) {
											if (exists $var_p_HC{$key}) {
												$AD = $AD.$var_p_HC{$key}.",";
											} else {
												$AD = $AD."0".",";
											}
										}
										$AD =~ s/\,$//;
										$output_string = $output_string."$GT_HC:$GQ:$DP:$AD:$CV"."\t";
									} else {
										if ($GT_HC eq $GT_free) {
											$CV = "HF";
										} else {
											$CV = "F";
										}
										$GQ = $free_fields[$free_GQ_index];
										$DP = $free_fields[$free_DP_index];
										#Allele depth is corrupted when merge, so when use HC GT, variants in free only
										# will be assigned as 0
										#Arranging AD to vars
										my %var_p_free;
										my @AD_free = split(",", $free_fields[$free_AO_index]);
										for (my $i=0;$i < scalar @AD_free; $i++ ) {
											$var_p_free{$free_vars[$i]} = $AD_free[$i];
										}
										$AD = $free_fields[$free_RO_index].",";
										foreach my $key ( sort { $vars{$a} <=> $vars{$b} } keys (%vars) ) {
											if (exists $var_p_free{$key}) {
												$AD = $AD.$var_p_free{$key}.",";
											} else {
												$AD = $AD."0".",";
											}
										}
										$AD =~ s/\,$//;
										$output_string = $output_string."$GT_free:$GQ:$DP:$AD:$CV"."\t";
									}
								} elsif ($HC_fields[$HC_GT_index] =~ m/0\/0/ && $free_fields[$free_GT_index] !~ m/0\/0/ 
									&& $free_alt_alleles[$free_bases[1]-1] >= $AA_sup_count_threshold
									&& $free_fields[$free_DP_index] >= $Depth_threshold) {
									$CV = "F";
									$GQ = $free_fields[$free_GQ_index];
									$DP = $free_fields[$free_DP_index];
									#Allele depth is corrupted when merge, so when use HC GT, variants in free only
									# will be assigned as 0
									#Arranging AD to vars
									my %var_p_free;
									my @AD_free = split(",", $free_fields[$free_AO_index]);
									for (my $i=0;$i < scalar @AD_free; $i++ ) {
										$var_p_free{$free_vars[$i]} = $AD_free[$i];
									}
									$AD = $free_fields[$free_RO_index].",";
									foreach my $key ( sort { $vars{$a} <=> $vars{$b} } keys (%vars) ) {
										if (exists $var_p_free{$key}) {
											$AD = $AD.$var_p_free{$key}.",";
										} else {
											$AD = $AD."0".",";
										}
									}
									$AD =~ s/\,$//;
									$output_string = $output_string."$GT_free:$GQ:$DP:$AD:$CV"."\t";
								} elsif ($HC_fields[$HC_GT_index] !~ m/0\/0/) {
									if ($GT_HC eq $GT_free) {
										$CV = "HF";
									} else {
										$CV = "H";
									}
									$GQ = $HC_fields[$HC_GQ_index];
									$DP = $HC_fields[$HC_DP_index];
									my %var_p_HC;
									my @AD_HC = split(",", $HC_fields[$HC_AD_index]);
									for (my $i=1;$i < scalar @AD_HC; $i++ ) {
										$var_p_HC{$HC_vars[$i-1]} = $AD_HC[$i];
									}
									$AD = $AD_HC[0].",";
									foreach my $key ( sort { $vars{$a} <=> $vars{$b} } keys (%vars) ) {
										if (exists $var_p_HC{$key}) {
											$AD = $AD.$var_p_HC{$key}.",";
										} else {
											$AD = $AD."0".",";
										}
									}
									$AD =~ s/\,$//;
									$output_string = $output_string."$GT_HC:$GQ:$DP:$AD:$CV"."\t";
								} else { #if both called reference then trust HC output
									if ($GT_HC eq $GT_free) {
										$CV = "HF";
									} else {
										$CV = "H";
									}
									$GQ = $HC_fields[$HC_GQ_index];
									$DP = $HC_fields[$HC_DP_index];
									my %var_p_HC;
									my @AD_HC = split(",", $HC_fields[$HC_AD_index]);
									for (my $i=1;$i < scalar @AD_HC; $i++ ) {
										$var_p_HC{$HC_vars[$i-1]} = $AD_HC[$i];
									}
									$AD = $AD_HC[0].",";
									foreach my $key ( sort { $vars{$a} <=> $vars{$b} } keys (%vars) ) {
										if (exists $var_p_HC{$key}) {
											$AD = $AD.$var_p_HC{$key}.",";
										} else {
											$AD = $AD."0".",";
										}
									}
									$AD =~ s/\,$//;
									$output_string = $output_string."$GT_HC:$GQ:$DP:$AD:$CV"."\t";
								}
							} elsif ($HC_sample !~ m/^\./) {
								my @HC_fields = split(":", $HC_sample);
								my ($GT, $GQ, $DP, $AD, $CV);
								$GT = $HC_fields[$HC_GT_index];
								my @bases = split ("/", $GT);
								my ($left, $right);
								if ($bases[0] == 0) { $left = 0;} else {
									$left = $vars{$HC_vars[$bases[0]-1]};
								}
								if ($bases[1] == 0) { $right = 0;} else {
									$right = $vars{$HC_vars[$bases[1]-1]};
								}
								my $GT_HC = $left."/".$right;
								$CV = "H";
								$GQ = $HC_fields[$HC_GQ_index];
								$DP = $HC_fields[$HC_DP_index];
								my %var_p_HC;
								my @AD_HC = split(",", $HC_fields[$HC_AD_index]);
								for (my $i=1;$i < scalar @AD_HC; $i++ ) {
									$var_p_HC{$HC_vars[$i-1]} = $AD_HC[$i];
								}
								$AD = $AD_HC[0].",";
								foreach my $key ( sort { $vars{$a} <=> $vars{$b} } keys (%vars) ) {
									if (exists $var_p_HC{$key}) {
										$AD = $AD.$var_p_HC{$key}.",";
									} else {
										$AD = $AD."0".",";
									}
								}
								$AD =~ s/\,$//;
								$output_string = $output_string."$GT_HC:$GQ:$DP:$AD:$CV"."\t";
							} elsif ($free_sample !~ m/^\./) {
								my @free_bases = split ("/", $free_fields[$free_GT_index]);
								my @free_alt_alleles = split(",", $free_fields[$free_AO_index]);
								
								if ($free_fields[$free_GT_index] =~ m/0\/0/ || 
									($free_alt_alleles[$free_bases[1]-1] >= $AA_sup_count_threshold
									&& $free_fields[$free_DP_index] >= $Depth_threshold) ) {
									my ($GT, $GQ, $DP, $AD, $CV);
									my @free_fields = split(":", $free_sample);
									$GT = $free_fields[$free_GT_index];
									my @bases = split ("/", $GT);
									my ($left, $right);
									if ($bases[0] == 0) { $left = 0; } else {
										$left = $vars{$free_vars[$bases[0]-1]};
									} 
									if ($bases[1] == 0) { $right = 0; } else {
										$right = $vars{$free_vars[$bases[1]-1]};
									}
									my $GT_free = $left."/".$right;
									$CV = "F";
									$GQ = $free_fields[$free_GQ_index];
									$DP = $free_fields[$free_DP_index];
									#Allele depth is corrupted when merge, so when use HC GT, variants in free only
									# will be assigned as 0
									#Arranging AD to vars
									my %var_p_free;
									my @AD_free = split(",", $free_fields[$free_AO_index]);
									for (my $i=0;$i < scalar @AD_free; $i++ ) {
										$var_p_free{$free_vars[$i]} = $AD_free[$i];
									}
									$AD = $free_fields[$free_RO_index].",";
									foreach my $key ( sort { $vars{$a} <=> $vars{$b} } keys (%vars) ) {
										if (exists $var_p_free{$key}) {
											$AD = $AD.$var_p_free{$key}.",";
										} else {
											$AD = $AD."0".",";
										}
									}
									$AD =~ s/\,$//;
									$output_string = $output_string."$GT_free:$GQ:$DP:$AD:$CV"."\t";
								} else {
									$NA_count += 1;
									$output_string = $output_string."./."."\t";
								}
							} else {
								$NA_count += 1;
								$output_string = $output_string."./."."\t";
							}
						}
					}
					#add NA counts
					$output_string =~ s/\t$//;
					if ($NA_count < (scalar @HC_temp - 9)) {
						return $output_string;
					} else {
						return 0;
					}
				} else {
					return FREE_filter($free_line);
				}
			}
		} else {
			return FREE_filter($free_line);
		}
	} else {
		return HC_filter($HC_line);
	}
}

sub usage {
    print "Unknown option: @_\n" if ( @_ );
    print "\nusage: VCF_combiner_HC_FREE.pl \n";
    print "--fai samtools index output fai file of the reference;\n";
    print "--FREE freebayes output.\n"; 
    print "--HC recalibrated HC output;\n";
    print "--MAF output file of the batch MAFs;\n";
    print "--output combined vcf output.\n";
    return(1);
}

sub read_file_line { #copied here: http://stackoverflow.com/questions/2498937/how-can-i-walk-through-two-files-simultaneously-in-perl
# not understanding what should be followed after shift, although I know it should be a file handler.
  my $fh = shift;
  if ($fh and my $line = <$fh>) {
    chomp $line;
    return $line;
  }
  return;
}

sub Define_Header {
	my ($input_FREE_vcf, $input_HC_vcf) = @_;
	my @output;
	my %ID;
	my $sample_line;
	open HC, "$input_HC_vcf" or die "Can not open the HC vcf file: $input_HC_vcf\n";
	while (my $line = <HC>) {
		chomp $line;
		if ($line =~ m/^\#/) {
			if ($line !~ m/^\#\#GATKCommandLine/ &&  $line !~ m/^\#\#contig/ &&  $line !~ m/^\#CHROM/) {
				my @eles = split(",", $line);
				if ( ! exists $ID{$eles[0]}) {
					$ID{$eles[0]} = 1;
					push @output, $line."\n";
				}
			} elsif ($line =~ m/^\#CHROM/) {
				$sample_line = $line;
			}
		} else {
			last;
		}
	}
	close HC;
	
	open FREE, "$input_FREE_vcf" or die "Can not reads the freebayes vcf file: $input_FREE_vcf\n";
	while (my $line = <FREE>) {
		chomp $line;
		if ($line =~ m/^\#/) {
			if ($line !~ m/^\#\#GATKCommandLine/ && $line !~ m/^\#\#fileformat/ &&  $line !~ m/^\#\#contig/ &&  $line !~ m/^\#CHROM/) {
				my @eles = split(",", $line);
				if ( ! exists $ID{$eles[0]}) {
					$ID{$eles[0]} = 1;
					push @output, $line."\n";
				}
			}
		} else {
			last;
		}
	}
	close FREE;
	push @output, '##FORMAT=<ID=CV,Number=.,Type=String,Description="define which caller called the variant H: HaplotypeCaller F: Freebayes">'."\n";
	push @output, $sample_line."\n";
	
	my $return_string="";
	foreach my $line (@output) {
		$return_string = $return_string.$line;
	}
	return $return_string;
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

sub BatchMAF {
	my ($vcf, $output_file)=@_; 
	my $DupSamples = "No"; # If there are duplicated samples included in VCFs, "Yes" mean the script will ignore those duplicated ones.
	
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
	
	
	#my $output_file = $vcf."_MAFs_".$exomecount."exomes.txt";
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
}







