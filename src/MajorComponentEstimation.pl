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
my @ATCG = ('A','T','C','G');
my $MINIMAL_FREQUENCY = 0.2;
my $MINIMAL_FREQUENCY_REF = 0.1;
my $MINREADS = 500;
my $MIN_IDENTITY = 0.8;
my $PATH2MURREF = '../data/mutfre/.'; #This is the path to the fold storing all .mutfre files for reference cell lines
my $PATH2DEPTH = '../data/testedcelllines_Table2/.';  #This is the path to the fold storing all .depth files for tested cell lines
#my $PATH2DEPTH = '../data/testedcelllines_Table3/.';  #This is the path to the fold storing all .depth files for tested cell lines
my $PWD = getcwd();

#--read in reference cell line genotype--#
chdir $PATH2MURREF or die;

my @file = glob("*.mutfre.txt");
my %cellLine_genotype;

foreach(@file){
	my $cellLine = $_;
	$cellLine =~ s/\.mutfre\.txt//g;	
	
	open FH, $_ or die;
	my %snp2genotype;
	while(<FH>){
		s/^\s+|\s+$//g;
		next if(/^\s*$|^pos|NA/);
		my @s = split(/\s+/);
		my $snp = shift @s;
		shift @s;
		my @genotype;
		for(my $i=0;$i<@s;$i++){
			if($s[$i] > $MINIMAL_FREQUENCY_REF){
				push(@genotype, $ATCG[$i]);
			}
		}
		$snp2genotype{$snp} = \@genotype;
	}
	close FH;
	
	$cellLine_genotype{$cellLine} = \%snp2genotype;
}


chdir $PWD or die;
chdir $PATH2DEPTH or die;
@file = glob("*.depth");

print "cellLine\t#match\t#nonmatch\t#total\t#identity\n";
foreach(@file){
	my $sample = $_;
	print "$_\t";
	open FH, $_ or die;	
	my %snp2genotype;
	
	while(<FH>){
		s/^\s+|\s+$//g;
		next if(/^\s*$|^pos|NA/);
		my @s = split(/\s+/);
		my $snp = shift(@s);
		shift(@s);
		my $sum = getSum(@s);
		next if($sum < $MINREADS);
		
		#--remove reads less than MINIMAL_FREQUENCY
		@s = removeReads(@s);

		#--get genotype and check if match--#
		my @genotype;
		for(my $i=0;$i<@s;$i++){
			if($s[$i] > $MINIMAL_FREQUENCY){
				push(@genotype, $ATCG[$i]);
			}
		}
		$snp2genotype{$snp} = \@genotype;		
	}close FH;
	
	#--loop through all cell lines to get best matches--#
	foreach(keys %cellLine_genotype){
		my $cellLine = $_;
		my @results = Identity(\%snp2genotype, $cellLine_genotype{$cellLine});
		my $match = $results[0];
		my $nonmatch = $results[1];
		my $total = $match + $nonmatch;
		my $identity = sprintf("%.4f", $results[2]);
		if($identity>$MIN_IDENTITY){
			print "\t$cellLine\t$match\t$nonmatch\t$total\t$identity";
		}
	}
	print "\n";
}


sub Identity{
	my %rh_snp2genotype = %{ shift @_ };
	my %rh_snp2genotype_ref = %{ shift @_ };
	my $match = 0;
	my $nonmatch = 0;
	
	foreach(keys %rh_snp2genotype){
		my $snp = $_;
		my @genotype = @{ $rh_snp2genotype{$snp} };
		
		next if(! exists $rh_snp2genotype_ref{$snp});
		my @genotype_ref = @{ $rh_snp2genotype_ref{$snp} };
		if(scalar @genotype != scalar @genotype_ref){
			$nonmatch++;
			next;
		}
		
		my $c=0;
		for(my $i=0; $i<@genotype; $i++){
			$c++ if($genotype[$i] ne $genotype_ref[$i]);
		}
		if($c==0){
			$match++;
		}else{
			$nonmatch++;
		}	
	}
	
	my $identity = $match/($match+$nonmatch);
	my @returned = ($match, $nonmatch, $identity);
	return @returned;
}

sub removeReads{
	my @t = @_;
	my $sum = getSum(@t);
	for(my $i=0; $i<@t; $i++){
		if($t[$i] < $sum* $MINIMAL_FREQUENCY){
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