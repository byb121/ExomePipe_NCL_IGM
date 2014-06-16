#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;

my $fastq1;
my $fastq2;
my $Results=GetOptions("1=s"=>\$fastq1, "2=s"=>\$fastq2);
my $output_file1 = $fastq1.".polyN_removed.txt";
my $output_file2 = $fastq2.".polyN_removed.txt";

my $i = 1;
my %shorts;
my ($id1, $id2, $reads1, $reads2, $quality_score1, $quality_score2);

open FASTQ_1, "$fastq1" or die "Cannot open fastq file $fastq1\n";
open FASTQ_2, "$fastq2" or die "Cannot open fastq file $fastq2\n";
open OUTPUT_1, ">$output_file1" or die "Cannot open the file $output_file1\n";
open OUTPUT_2, ">$output_file2" or die "Cannot open the file $output_file2\n";
while (my $line1=<FASTQ_1>) {
	
	my $line2 = <FASTQ_2>;
	
	if ($i ==1 && $line1 =~ m/^\@/) {
		chomp $line1;
		chomp $line2;
		$id1 = $line1;
		$id2 = $line2;
		$i += 1;
		next;
	} elsif ($i ==2) {
		chomp $line1;
		$reads1 = $line1;
		while ($reads1 =~ m/N$/) {
			$reads1 =~ s/N$//;
		}
		chomp $line2;
		$reads2 = $line2;
		while ($reads2 =~ m/N$/) {
			$reads2 =~ s/N$//;
		}
		$i += 1;
		next;
	} elsif ($i ==3 ) {
		$i += 1;
		next;
	} elsif ($i == 4) {
		chomp $line1;
		chomp $line2;
		if (length($reads1) >= 20 && length($reads2) >= 20) {
			$quality_score1 = substr ($line1, 0, length($reads1));
			$quality_score2 = substr ($line2, 0, length($reads2));	
			print OUTPUT_1 $id1."\n".$reads1."\n"."+\n".$quality_score1."\n";
			print OUTPUT_2 $id2."\n".$reads2."\n"."+\n".$quality_score2."\n";
		}
		
		$i = 1;
		$id1 = "";
		$id2 = "";
		$reads1 = "";
		$reads2 = "";
		$quality_score1 = "";
		$quality_score2 = "";
	} else {
		print "Fatal error: either fastq format is incorrct or software failure.\n";
	}
}

close (FASTQ_1);
close (FASTQ_2);
close (OUTPUT_1);
close (OUTPUT_2);

exit;
