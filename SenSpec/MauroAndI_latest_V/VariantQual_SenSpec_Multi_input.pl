#!/usr/bin/perl;
use strict;
use warnings;
use Getopt::Long;
use Cwd qw(abs_path);	
# use IPC::System::Simple qw(system);
	
my ($help, $vcf, $ref);
my $currentdir=`pwd`;
my $CurrentScript=abs_path($0);
$CurrentScript=~/(.*)\/VariantQual/;
my $scriptdir=$1; #="~/scripts/darren";#### Path to R script ###
chomp $currentdir;

my $home_dir = $ENV{"HOME"};
chomp $home_dir;

#my $commonTail=".vcf";
#my $samplePath="temp";
#my $sampleNames="Sample_1";

my $commonTail;
my $samplePath;
my $sampleNames;
my $commonHeader="";

usage() if ( @ARGV < 1 || ! GetOptions('help|?' => \$help, 'sampleNames=s' => \$sampleNames, 'refSNPs=s' => \$ref, 'commonTail=s' => \$commonTail, 'samplePath=s' => \$samplePath,  'commonHeader=s' => \$commonHeader) || defined $help );

unless (defined $sampleNames) {
	die "You have not supplied sample names using --sampleNames\n";
}

unless (defined $ref) {
	die "You have not supplied a reference list of SNPs using --refSNPs\n";
}

unless (defined $samplePath) {
	die "You have not supplied the path that have all vcf files of samples using --samplePath\n";
}

unless (defined $commonTail) {
	die "You have not supplied common tail of vcf file names using --commonTail\n";
}

$samplePath =~ s/^\~/$home_dir/;
$samplePath = abs_path($samplePath);
$samplePath =~ s/\/$//;

if ($sampleNames !~ /\,/) {
	$vcf = $samplePath."/".$commonHeader.$sampleNames.$commonTail;
	print "processing file: $vcf\n";
	$vcf = abs_path($vcf);
	my $input=Filter($vcf, $ref);
	r($input,$scriptdir);
} else {
	my @samples=split(",", $sampleNames);
	
	my $output = $samplePath."/".$samples[0]."_to_".$samples[$#samples]."_SenSpec.txt";
	my $sample_string="";
	my $sens_string="";
	my $spec_string="";
	
	foreach my $sample (@samples) {
		$vcf = $samplePath."/".$commonHeader.$sample.$commonTail;
		print "processing file: $vcf\n";
		$vcf = abs_path($vcf);
		my $input=Filter($vcf, $ref);
		r($input,$scriptdir);
		my $routput=$vcf.".SensSpec.txt";
		open Routput, "$routput" or die "Cannot open the file $routput";
		my @lines = <Routput>;
		my @temp1 = split("\t", $lines[0]);
		my @temp2 = split("\t", $lines[1]);
		
		if ($sample_string eq "") {
			$sample_string = $sample;
		} else {
			$sample_string = $sample_string."\t".$sample;
		}
		
		if ($sens_string eq "") {
			chomp $temp1[1];
			$sens_string = $temp1[1];
		} else {
			chomp $temp1[1];
			$sens_string = $sens_string."\t".$temp1[1];
		}
		
		if ($spec_string eq "") {
			chomp $temp2[1];
			$spec_string = $temp2[1];
		} else {
			chomp $temp2[1];
			$spec_string = $spec_string."\t".$temp2[1];
		}
		close Routput;
	}
	
	open OUTPUT, ">$output" or die "Cannot open the file $output to output results.\n";
	chomp $sample_string;
	chomp 	$sens_string;
	chomp 	$spec_string;
	print OUTPUT $sample_string."\n";
	print OUTPUT $sens_string."\n";
	print OUTPUT $spec_string."\n";
	close OUTPUT;
}
exit;

####SUBROUTINES#################################################################################

sub Filter{ 
	my $File=$_[0];
	my $ref=$_[1];
	my %VariantHash=();
	my @Array=();
	
	open FILE, "$File";
	while(<FILE>){
		chomp $_;
		if($_=~/\#/){
			next;
		}else{
			my @SplitLine=split(/\t/, $_);
			my $Chr=$SplitLine[0];
			unless($SplitLine[0]=~/chr\S+/){
				$Chr="chr".$SplitLine[0];
			}
			unless($SplitLine[6]=~/^PASS$/ || $SplitLine[6]=~/^\.$/ ){
				next;
			}
			my $Pos=$SplitLine[1];
			my $Ref=uc($SplitLine[3]);
			my $Variant=uc($SplitLine[4]);
			my @SplitInfo=split(';', $SplitLine[7]);
			my @SplitGT=split(':', $SplitLine[9]);
			my $GT=$SplitGT[0];
			my $Match=$Chr."_".$Pos."_".$Ref;
			if($Variant=~/\S+\,/){
				my @SplitVar=split(',', $Variant);
				foreach my $Var(@SplitVar){
					if($GT eq "1/1"){
						$VariantHash{$Match}="3";
					}elsif($GT eq "0/1"){
						$VariantHash{$Match}="2";
					}
				}
			}else{
				if($GT eq "1/1"){
					$VariantHash{$Match}="3";
				}elsif($GT eq "0/1"){
					$VariantHash{$Match}="2";
				}					
			}

		}
	}	
	close FILE;
	my $rinput="$File.Input.txt";
	open I, ">$rinput";
	print I "rsNumber\tChange\tChromosome\tPosition\tFrequency\tSample\n";
	
	open REF, "$ref";
	while(<REF>){
		chomp $_;
		if($_=~/rsNumber\t/){
			next;
		}else{
			my @Split=split('\t', $_);
			my @splitrefvar=split('/', $Split[1]);
			my $VariantMatch=$Split[2]."_".$Split[3]."_".$splitrefvar[0];
			if(exists $VariantHash{$VariantMatch}){
				print I $Split[0]."\t".$Split[1]."\t".$Split[2]."\t".$Split[3]."\t".$Split[4]."\t".$VariantHash{$VariantMatch}."\n";
			}else{		
				print I $Split[0]."\t".$Split[1]."\t".$Split[2]."\t".$Split[3]."\t".$Split[4]."\t"."1"."\n";
			}
		}
	}
	
	close REF;
	close I;
	return($rinput);
}

sub r{
	my ($in,$scriptdir)=@_;
	my $invokeR = "R --vanilla --silent --no-save --args";
	my $scriptR = "$scriptdir/VariantQual.r > $in.R.log";
	my $Garbage=system("$invokeR $in < $scriptR ");
	system("rm $in $in.R.log");
	return(1);
}

sub usage {
    print "Unknown option: @_\n" if ( @_ );
    print "\nusage: VariantQual_SenSpec_Multi_inout.pl [--sampleNames Sample_1,Sample_2,...][--samplePath Path of the folder of vcfs] [--commonTail common tails of file names exclude the sample name part]  [-refSNPs reference list of SNPs] [--commonHeader if the files has common headers before sample names][-help|-?]\n\n";
	return(1);
}

