#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;

###############################################
# Emit vary tables of coverage summary on exons of a gene list
# include average coverage on exons of each sample
# more than a certain times covered percentages of each sample
###############################################

my $Output_prefix="prefix_test";
my $SamplePath="/users/a5907529/lustre/Kidney_20130211/A1569_Fastq/";
#my $SampleNames="Sample_1,Sample_2,Sample_3,Sample_4,Sample_5,Sample_6,Sample_7,Sample_8,Sample_D09254,Sample_D09255,Sample_D115953,Sample_D116785,Sample_D118160,Sample_D156602,Sample_D21608,Sample_D21683,Sample_D44659,Sample_D53464,Sample_D53686,Sample_D53687,Sample_D62086,Sample_D62088,Sample_D74357,Sample_D79574";
my $SampleNames="Sample_1,Sample_2";
my $CovFileName="test.txt";

my $Results=GetOptions("Output_prefix=s"=>\$Output_prefix, "SamplePath=s"=>\$SamplePath, "SampleNames=s"=>\$SampleNames, "CovFileName=s"=>\$CovFileName);

my @Names = split(",", $SampleNames);
$SamplePath =~ s/\/$//;
my $output_file_mean_coverage = $SamplePath."/".$Output_prefix."mean_coverage_on_exons_of_genes.txt";
my $output_file_1x_coverage = $SamplePath."/".$Output_prefix."1x_coverage_on_exons_of_genes.txt";
my $output_file_10x_coverage = $SamplePath."/".$Output_prefix."10x_coverage_on_exons_of_genes.txt";
my $output_file_20x_coverage = $SamplePath."/".$Output_prefix."20x_coverage_on_exons_of_genes.txt";
my $output_file_40x_coverage = $SamplePath."/".$Output_prefix."40x_coverage_on_exons_of_genes.txt";

print AverageCoverageOnExonsOfSamples($SamplePath, \@Names, $CovFileName, $output_file_mean_coverage);
print FractionOfLeastCoverageOnExonsOfSamples($SamplePath,  \@Names, 1, $CovFileName, $output_file_1x_coverage);
print FractionOfLeastCoverageOnExonsOfSamples($SamplePath,  \@Names, 10, $CovFileName, $output_file_10x_coverage);
print FractionOfLeastCoverageOnExonsOfSamples($SamplePath,  \@Names, 20, $CovFileName, $output_file_20x_coverage);
print FractionOfLeastCoverageOnExonsOfSamples($SamplePath,  \@Names, 40, $CovFileName, $output_file_40x_coverage);

exit;

sub AverageCoverageOnExonsOfSamples{
	
	my ($SamplePath, $names_ref, $CovFileName, $output_file) = @_;
	
	my @Names = @{$names_ref};
	
	my %mean_count_on_exons;
	
	foreach my $sample (@Names) {
		my $file = $SamplePath."/".$sample.$CovFileName;
		open INPUT, $file or die "coverage_summary_on_gene_list_on_file_list.pl: Cannot open ".$file."\n";
		while (my $line =<INPUT>) {
			if ($line =~ m/^chr/) {
				chomp $line;
				my @words = split(/\t/,$line);
				my $temp_name=$words[3];
				if ( exists $mean_count_on_exons{$temp_name}{$sample} ) {
					$mean_count_on_exons{$temp_name}{$sample} = $mean_count_on_exons{$temp_name}{$sample} + $words[4]*$words[7];
				} else {
					$mean_count_on_exons{$temp_name}{$sample} = $words[4]*$words[7];
				}
			}
		}
		close INPUT;
		print "coverage_summary_on_gene_list_on_file_list.pl: $file\n";
	}
	
	open OUTPUT, ">$output_file" or die "coverage_summary_on_gene_list_on_file_list.pl: Cannot open ".$output_file." to output.\n";
	print OUTPUT "Region_Name";
	foreach my $sample (@Names) {
		print OUTPUT "\t".$sample;
	}
	print OUTPUT "\n";
	foreach my $exon (sort keys %mean_count_on_exons) {
		print OUTPUT $exon;
		foreach my $sample (@Names) {
				print OUTPUT "\t".sprintf( "%.1f",$mean_count_on_exons{$exon}{$sample});
		}
		print OUTPUT "\n";
	}
	close OUTPUT;
	
	return "Success.\n";
}

sub FractionOfLeastCoverageOnExonsOfSamples {
	
	my ($SamplePath, $names_ref, $least_coverage, $CovFileName, $output_file ) = @_;
	
	my @Names = @{$names_ref};
	
	my %fraction_on_exons;
	
	foreach my $sample (@Names) {
		my $file = $SamplePath."/".$sample.$CovFileName;
		open INPUT, $file or die "coverage_summary_on_gene_list_on_file_list.pl: Cannot open ".$file."\n";
		while (my $line =<INPUT>) {
			if ($line =~ m/^chr/) {
				chomp $line;
				my @words = split(/\t/,$line);
				my $temp_name=$words[3];
				if ( $words[4] >= $least_coverage ){
					if ( exists $fraction_on_exons{$temp_name}{$sample} ) {
						$fraction_on_exons{$temp_name}{$sample} = $fraction_on_exons{$temp_name}{$sample} + $words[7];
					} else {
						$fraction_on_exons{$temp_name}{$sample} = $words[7];
					}
				} else {
					if ( !exists $fraction_on_exons{$temp_name}{$sample} ) {
						$fraction_on_exons{$temp_name}{$sample} = 0
					}
				}
			}
		}
		close INPUT;
		print "coverage_summary_on_gene_list_on_file_list.pl: $file\n";
	}
	
	open OUTPUT, ">$output_file" or die "coverage_summary_on_gene_list_on_file_list.pl: Cannot open ".$output_file." to output.\n";
	print OUTPUT "Region_Name";
        foreach my $sample (@Names) {
                print OUTPUT "\t".$sample;
        }
        print OUTPUT "\n";
	foreach my $exon (sort keys %fraction_on_exons) {
		print OUTPUT $exon;
		foreach my $sample (@Names) {
				print OUTPUT "\t".sprintf( "%.1f", $fraction_on_exons{$exon}{$sample}*100 );
		}
		print OUTPUT "\n";
	}
	close OUTPUT;
	return "Success.\n";
}


