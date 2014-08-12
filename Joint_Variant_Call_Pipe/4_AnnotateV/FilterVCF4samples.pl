#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;

print "\n";
print "##########################################################################################\n";
print "#  output a new vcf file fitlered by named samples (multiple samples separated with (,)  #\n";
print "#               reference calls will be removed. Filter values are ignored.              #\n";
print "##########################################################################################\n";


my ($vcf, $sample_string, $output);

my $help;
usage() if ( @ARGV < 1 || ! GetOptions('help|?' => \$help, "vcf=s"=>\$vcf, "samples=s"=>\$sample_string,
'output=s' => \$output ) || defined $help );

my @samples = split(",", $sample_string);
my %index;

open VCF, "$vcf" or die "Cannot open the file $vcf, check if it exists\n";
open OUT, ">$output" or die "Cannot open the file $output\n";

my $line_count=0;
Line: while ( my $line = <VCF> ) {
	chomp $line;
	$line_count += 1;
	if($line_count%10000 == 0) {
		print "processed $line_count lines.\n";
	}
	if ($line =~ m/^\#/) {
        if ($line =~ m/^\#CHROM\t/) {
            my $temp_string="";
            my @temp = split("\t", $line);
            for (my $i = 0; $i < 9; $i++) {
                $temp_string=$temp_string."$temp[$i]\t";
            }
            foreach my $sample (@samples) {
                for (my $i = 9; $i < scalar @temp; $i++) {
                    if ($temp[$i] =~ m/^$sample$/) {
                        print "Found sample $sample on column $i.\n";
                        $index{$i} = 1;
                        $temp_string=$temp_string."$sample\t";
                    }
                }
            }
    
            if (scalar keys(%index) != scalar @samples) {
                print "Error: input samples were not found in the vcf.\n";
                exit;
            } else {
                $temp_string =~ s/\t$//;
                print OUT $temp_string."\n";
            }
        } else {
            print OUT $line."\n";
        }
    } else {
        my $na_count = 0;
        my $ref_count = 0;
		my @temp = split("\t", $line);
		my $format = $temp[8];
		my @format_fields = split(":", $format);
		my $GT_index;
		for (my $i=0;$i<scalar(@format_fields);$i++){
			if ($format_fields[$i] =~ m/GT/) {
				$GT_index = $i;
				#print $line."$GT_index\n";
			}
		}
		
        my $temp_string = "";
		if (defined $GT_index) {
            for (my $i=0;$i<9;$i++) {
                $temp_string=$temp_string."$temp[$i]\t";
            }
			for(my $i=9;$i< scalar @temp; $i++){
                if (exists $index{$i}) {
                    my $sample = $temp[$i];
                    #print "use column $i\n";
                    $temp_string=$temp_string."$sample\t";
                    if ($sample !~ m/^\./) {
                        my @sample_fields = split(":", $sample);
                        if ($sample_fields[$GT_index] =~ m/^0\/0/) {
                            $ref_count += 1;
                        }
                    } else {
                        $na_count += 1;
                    }
                }
            }
		} else {
            #print "skip line without GT tag\n";
            next Line;
        }
        
        if (scalar @samples > ($na_count+$ref_count) ) {
		        #my $sum=$na_count+$ref_count;
		        #print "$line_count: na plus ref: $sum\n";
		        $temp_string =~ s/\t$//;
			print OUT $temp_string."\n";
        	}
	}
}

close VCF;
close OUT;
print "Done!\n";

exit;

sub usage {
    print "Unknown option: @_\n" if ( @_ );
    print "\nusage: FilterVCF4samples.pl: \n";
    print "--vcf input vcf file;\n";
    print "--samples sample name, or multiple names separated with (,).\n";
    print "--output output file name.\n";
    return(1);
}

