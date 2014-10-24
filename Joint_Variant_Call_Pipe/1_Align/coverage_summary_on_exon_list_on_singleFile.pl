#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;

###############################################
# Emit vary tables of coverage summary on exons of a gene list
# include average coverage on exons of each sample
# more than a certain times covered percentages of each sample
###############################################

my $CovFileName;
my $output_file;

my $Results=GetOptions("CovFileName=s"=>\$CovFileName, "output=s"=>\$output_file);

my $output_file_mean_coverage_ref = AverageCoverageOnExonsOfSamples($CovFileName);
my $output_file_1x_coverage_ref = FractionOfLeastCoverageOnExonsOfSamples(1, $CovFileName);
my $output_file_5x_coverage_ref = FractionOfLeastCoverageOnExonsOfSamples(5, $CovFileName);
my $output_file_10x_coverage_ref = FractionOfLeastCoverageOnExonsOfSamples(10, $CovFileName);
my $output_file_20x_coverage_ref = FractionOfLeastCoverageOnExonsOfSamples(20, $CovFileName);
my $output_file_40x_coverage_ref = FractionOfLeastCoverageOnExonsOfSamples(40, $CovFileName);

my %output_file_mean_coverage = %$output_file_mean_coverage_ref;
my %output_file_1x_coverage = %$output_file_1x_coverage_ref;
my %output_file_5x_coverage = %$output_file_5x_coverage_ref;
my %output_file_10x_coverage = %$output_file_10x_coverage_ref;
my %output_file_20x_coverage = %$output_file_20x_coverage_ref;
my %output_file_40x_coverage = %$output_file_40x_coverage_ref;

my @keys = keys %output_file_mean_coverage;
my $numberOfExons = scalar @keys;
my $total_mean_cov = 0;
my $numberOfExonsOver1Xcoverage = 0;
my $numberOfExonsOver5Xcoverage = 0;
my $numberOfExonsOver10Xcoverage = 0;
my $numberOfExonsOver20Xcoverage = 0;
my $numberOfExonsOver40Xcoverage = 0;


open OUTPUT, ">$output_file" or die "Cannot open the file to ouput $output_file\n";
print OUTPUT "Exon\t1X_Cov_frac\t5X_Cov_frac\t15X_Cov_frac\t20X_Cov_frac\t40X_Cov_frac\n";
foreach my $exon (sort keys %output_file_mean_coverage) {
	if (exists $output_file_1x_coverage{$exon} && exists $output_file_5x_coverage{$exon}
		&& exists $output_file_10x_coverage{$exon} && exists $output_file_20x_coverage{$exon}
		&& exists $output_file_40x_coverage{$exon}) {
		$total_mean_cov += $output_file_mean_coverage{$exon};
		if ($output_file_1x_coverage{$exon} >= 0.9) {
			$numberOfExonsOver1Xcoverage += 1;
		}
		if ($output_file_5x_coverage{$exon} >= 0.9) {
			$numberOfExonsOver5Xcoverage += 1;
		}
		if ($output_file_10x_coverage{$exon} >= 0.9) {
			$numberOfExonsOver10Xcoverage += 1;
		}
		if ($output_file_20x_coverage{$exon} >= 0.9) {
			$numberOfExonsOver20Xcoverage += 1;
		}
		if ($output_file_40x_coverage{$exon} >= 0.9) {
			$numberOfExonsOver40Xcoverage += 1;
		}
		print OUTPUT $exon."\t".$output_file_mean_coverage{$exon}."\t".$output_file_1x_coverage{$exon}
			."\t".$output_file_5x_coverage{$exon}."\t".$output_file_10x_coverage{$exon}
			."\t".$output_file_20x_coverage{$exon}."\t".$output_file_40x_coverage{$exon}."\n";	
	} else {
		print "Something is wrong, the exon does not exists for the other coverage calculation. exit\n";
		exit;
	}
}

print OUTPUT "Number_of_exons\t$numberOfExons\n";
print OUTPUT "Total_mean\t".$total_mean_cov/$numberOfExons."\n";
print OUTPUT "Number_of_exons_at_least_90Percent_of_exon_covered_at_least_1X\t$numberOfExonsOver1Xcoverage\n";
print OUTPUT "Number_of_exons_at_least_90Percent_of_exon_covered_at_least_5X\t$numberOfExonsOver5Xcoverage\n";
print OUTPUT "Number_of_exons_at_least_90Percent_of_exon_covered_at_least_10X\t$numberOfExonsOver10Xcoverage\n";
print OUTPUT "Number_of_exons_at_least_90Percent_of_exon_covered_at_least_20X\t$numberOfExonsOver20Xcoverage\n";
print OUTPUT "Number_of_exons_at_least_90Percent_of_exon_covered_at_least_40X\t$numberOfExonsOver40Xcoverage\n";

close OUTPUT;
print "Done.\n";
exit;

sub AverageCoverageOnExonsOfSamples{
	my ($CovFileName) = @_;
	my %mean_count_on_exons;
	my $file = $CovFileName;
	print "reading: $file\n";
	open INPUT, $file or die "coverage_summary_on_gene_list_on_file_list.pl: Cannot open ".$file."\n";
	while (my $line =<INPUT>) {
		if ($line !~ m/^all\t/) {
			chomp $line;
			my @words = split(/\t/,$line);
			if ( exists $mean_count_on_exons{$words[3]} ) {
				$mean_count_on_exons{$words[3]} = $mean_count_on_exons{$words[3]} + $words[4]*$words[7];
			} else {
				$mean_count_on_exons{$words[3]} = $words[4]*$words[7];
			}
		}
	}
	close INPUT;
	return \%mean_count_on_exons;
}

sub FractionOfLeastCoverageOnExonsOfSamples {
	my ($least_coverage, $CovFileName) = @_;
	my %fraction_on_exons;
	my $file =  $CovFileName;
	print "reading: $file\n";
	open INPUT, $file or die "coverage_summary_on_gene_list_on_file_list.pl: Cannot open ".$file."\n";
	while (my $line =<INPUT>) {
		if ($line !~ m/^all\t/) {
			chomp $line;
			my @words = split(/\t/,$line);
			if ( $words[4] >= $least_coverage ){
				if ( exists $fraction_on_exons{$words[3]} ) {
					$fraction_on_exons{$words[3]} = $fraction_on_exons{$words[3]} + $words[7];
				} else {
					$fraction_on_exons{$words[3]} = $words[7];
				}
			} else {
				if ( !exists $fraction_on_exons{$words[3]} ) {
					$fraction_on_exons{$words[3]} = 0
				}
			}
		}
	}
	close INPUT;
	return \%fraction_on_exons;
}


