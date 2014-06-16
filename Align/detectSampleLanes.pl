#!/usr/bin/perl
use strict;
use warnings;

my ($sample_path, $sample_id)  = @ARGV;
$sample_path =~ s/\/$//;
my $folder=$sample_path."/".$sample_id;
my $lanes="";
opendir (DIR, $folder) or die "Cann't find the directory $folder";
while (my $file = readdir(DIR)) {
	next if ($file =~ m/^\./);
	
	if($file =~ m/$sample_id\_L00(\d+)\_R1.+\.fastq$/) {
	
	#/*********************** if QC ed with Trimming.sh
	#if($file =~ m/$sample_id\_L00(\d+)\_1.+\.fq$/) {
	#***********************/ 
		$lanes=$1." ".$lanes;
	}
}
close DIR;
chomp $lanes;
print $lanes."\n";
exit;

	
