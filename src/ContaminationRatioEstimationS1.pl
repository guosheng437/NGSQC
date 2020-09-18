######################################################
# Sheng Guo, PhD (guosheng@crownbio.com)
# July 2020
# Copyright (C) 2020 Sheng Guo, Crown Bioscience Inc.
######################################################
#!/usr/bin/perl -w
use strict;
use Cwd;

######################################################
#--global parameters--#
######################################################
my $MINREADS = 1000;
my $LARGE_SMALL_CUTOFF = 0.20;
my $MINIMAL_FREQUENCY = 0.01; #10/1000 as the error cut-off
my $PATH2DEPTH = '../data/testedcelllines_Table2/.';  #This is the path to the fold storing all .depth files for tested cell lines
#my $PATH2DEPTH = '../data/testedcelllines_Table3/.';  #This is the path to the fold storing all .depth files for tested cell lines
my $PWD = getcwd();


my %READ2CUT;
#--read in binomial threshold--#
open FH, "binomial_threshold.txt" or die;
while(<FH>){
	s/^\s+|\s+$//g;
	next if(/^#|^total/);
	my @s = split(/\s+/);
	$READ2CUT{$s[0]}=$s[1];
	
}close FH;


chdir $PATH2DEPTH or die;
my @file = glob("*.depth");

print "cellLine\t#informativeSNPs\tCoarseEstimationOfcontamination/heterogeneityRatio(%)\n";
foreach(@file){
	my $sample = $_;
	#print $_, "\n";
	
	my $outfile = $_;
	$outfile =~ s/depth/SNPratio/g;
	open OUT, ">$outfile" or die;
	
	open FH, $_ or die;
	my @ratios;
	while(<FH>){
		s/^\s+|\s+$//g;
		next if(/^\s*$|^pos|NA/);
		my @s = split(/\s+/);
		shift(@s);
		shift(@s);
		my $sum = getSum(@s);
		next if($sum < $MINREADS);
		
		#--remove reads less than MINIMAL_FREQUENCY
		@s = removeReads(@s);
		$sum = getSum(@s);
		
		#--determine large and small nucleotides
		my @idx_large = largeNucleotideIndx(@s);
		my @idx_small = smallNucleotideIndx(@s);
		
		next if(@idx_small==0);		
		
		#--4 scenarios
		#--case 1: large A and small T
		if(@idx_large==1 && @idx_small == 1){
			my $large_count = $s[ $idx_large[0] ]; 
			my $small_count = $s[ $idx_small[0] ]; 
			my $r = 1.5 * $small_count/($small_count + $large_count);
			push(@ratios, $r);
			next;
		}
		
		#--case 2: large A, small T and small G
		if(@idx_large==1 && @idx_small == 2){
			my $large_count = $s[ $idx_large[0] ]; 
			my $small_count = $s[ $idx_small[0] ] + $s[ $idx_small[1] ]; 
			my $r = $small_count/($small_count + $large_count);
			push(@ratios, $r);
			next;
		}		
		
		#--case 3: large A and T, small G
		if(@idx_large==2 && @idx_small == 1){
			my $large_count = $s[ $idx_large[0] ] + $s[ $idx_large[1] ]; 
			my $small_count = $s[ $idx_small[0] ]; 
			my $r = 1.5 * $small_count/($small_count + $large_count);
			push(@ratios, $r);
			next;
		}		
		
		#--case 4: large A and T, small G and C
		if(@idx_large==2 && @idx_small == 2){
			my $large_count = $s[ $idx_large[0] ] + $s[ $idx_large[1] ]; 
			my $small_count = $s[ $idx_small[0] ] + $s[ $idx_small[1] ]; 
			my $r =  $small_count/($small_count + $large_count);
			push(@ratios, $r);
			next;
		}	
		
	}close FH;
	
	#--get median ratios--#
	@ratios = sort{$a<=>$b}(@ratios);
	my $median = sprintf("%.2f", 100*$ratios[ @ratios/2  ]); 
	print "$sample\t", scalar @ratios, "\t$median\n";

	foreach(@ratios){
	 printf OUT ("%.4f\n",$_);
	}
	close OUT;
	
}

sub largeNucleotideIndx{
	my @t = @_;
	my $sum = getSum(@t);
	my @idx;
	for(my $i=0; $i<@t; $i++){
		if($t[$i] >= $sum* $LARGE_SMALL_CUTOFF){
			push(@idx, $i);
		}
	}
	return @idx;
}

sub smallNucleotideIndx{
	my @t = @_;
	my $sum = getSum(@t);
	my @idx;
	for(my $i=0; $i<@t; $i++){
		if($t[$i]>0 && $t[$i] < $sum* $LARGE_SMALL_CUTOFF){
			push(@idx, $i);
		}
	}
	return @idx;
}

sub removeReads{
	my @t = @_;
	my $sum = getSum(@t);
	for(my $i=0; $i<@t; $i++){
		if($t[$i] <= $READ2CUT{$sum}){
			$t[$i] = 0;
		}
	}
	return @t;
}

sub getSum{
	my $x = 0;
	foreach(@_){
		$x = $x + $_;
	}
	return $x;
}

