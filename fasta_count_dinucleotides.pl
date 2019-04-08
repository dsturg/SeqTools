#!/usr/bin/perl

# Prints out perc. GC in a multi-fasta file.


# Usage: perl fasta_count_dinucleotides.pl [filename] (repeatmask [mask/unmask])
#

use strict;
#use FileHandle;
#use POSIX qw(ceil floor);

my $filename = $ARGV[0];
my $repeatmask = $ARGV[1];

chomp $filename;

my $count = 0;

my $DNA;
my @DNA;
my $base;
my $errors = 0;
my $totalcheck = 0;
my $count_of_A = 0;
my $count_of_C = 0;
my $count_of_G = 0;
my $count_of_T = 0;
my $count_of_N = 0;
my $masked = 0;
my $gc = 0;
my $percgc = 0;
my $nucleotides = 0;
my $freqCG = 0;
######################
#@@@@@@@@@@@@@@@@@@@@@
######################
my $datafile;
my @datafile;
my $seq;
my %seqhash;
my $newseq = 0;

my $id;
my $idminus1;

# This is added to a
#my $

unless ($repeatmask =~ /mask|unmask/) {
	die "Invalid option for repeatmask: [mask/unmask]";
}

open(DATAFILE, "$filename");

while ($datafile = <DATAFILE>) {

#! This section changes depending on identifier structure
#! The most generic one is on the bottom

	#if ($datafile =~ /^>lcl\|(NM_[0-9]+).+/) {
	#if ($datafile =~ /^>(FBgn[0-9]+).+/) {
	#if ($datafile =~ /^>([a-zA-Z0-9\.\_]+)/) {
	#if ($datafile =~ /^>(\S+)/) {
	if ($datafile =~ /^>(.+)/) {
	
		$newseq += 1;
		$idminus1 = $id;
		$id = $1; 
		#print "\n*** $promoter **\n";
		if ($newseq > 1) {
			$seq = join "", @datafile;
			$seq =~ s/\n//g;
			$seq =~ s/\s//g;
			if (exists ($seqhash{$idminus1})) {
				print "Warning!  The id is non-unique\n";
				quit();
			}
			$seqhash{$idminus1} = $seq;			
			@datafile = ();
		}

	} else {
	
	push (@datafile, $datafile);

	}

}

if (exists ($seqhash{$id})) {
	print "Warning!  The id is non-unique\n";
	quit();
}
$seq = join "", @datafile;
$seq =~ s/\n//g;
$seq =~ s/\s//g;
$seqhash{$id} = $seq;			

my $k;

my @dinucs = qw(AA AC AG AT CA CC CG CT GA GC GG GT TA TC TG TT);


print "ID\tA\tC\tG\tT\tN\tnonACGT\tmasked\tNUCS\tPercGC\t";
print join("\t", @dinucs), "\tfreqCG\n";

foreach $k (sort keys %seqhash) {
    $count_of_A = 0;
	$count_of_C = 0;
	$count_of_G = 0;
	$count_of_T = 0;
	$count_of_N = 0;
	$errors = 0;
	$totalcheck = 0;

     
	if ($repeatmask eq "mask") {
		$DNA = $seqhash{$k};
	} else {
		$DNA = uc $seqhash{$k};
	}
	@DNA = split('',$DNA);
	foreach $base(@DNA) {
		if ($base eq 'A') {
			++$count_of_A;
		} elsif ($base eq 'C') {
			++$count_of_C;
		} elsif ($base eq 'G') {
			++$count_of_G;
		} elsif ($base eq 'T') {
			++$count_of_T;
		} elsif ($base eq 'N') {
			++$count_of_N;
			++$errors;
		} elsif ($base =~ /[acgt]/) {
			++$masked;
		} else {
			++$errors;
		}
	}
		

	$count = length($DNA);
	$nucleotides = $count_of_A + $count_of_C + $count_of_G + $count_of_T;

	# Print the results
	print $k, "\t";
	print $count_of_A , "\t"; 
	print  $count_of_C , "\t"; 
	print $count_of_G , "\t"; 
	print $count_of_T , "\t"; 
	print $count_of_N , "\t"; 
	print $errors , "\t"; 
	print $masked , "\t"; 
	print $nucleotides, "\t";
	$totalcheck = $count_of_C + $count_of_G + $count_of_T + $count_of_A + $errors + $masked;
	#if ($totalcheck ne $count) {
	#	print "Error in counting\n";
	#}
	$percgc = 0;
	if ($nucleotides > 0) {
		$percgc = ($count_of_C + $count_of_G) / $nucleotides * 100;
	} 
	
	$percgc = sprintf("%.4f", $percgc);
	print $percgc, "\t";
	
	#######
	#my %count;
	#my $sequence = $DNA;
	#print "setting sequence ewul\n";
	#print "length\t", length $DNA;
	#my $sequence = 'AACGTACTGACGTACTGGTAAATGGTACGA';
	#my $sequence ='AACGTACTGA CGTACTGGTA AATGGTACGA';

	my $dinuc;
	my @results = ();
	my $dinuc_count;

	##########
	# Note that it comes out the same in EMBOSS COMPSEQ
	##########

	##########
	# This solution comes from 
	# Perl for bioinformatics 2nd edition (Arun Jagota)
	##########
	
	my $i;
	my $len;
	my %dimer_counts;
	my $let1;
	my $let2;
	my $dinucs_sum = 0;

	my @res=split(//,$DNA); $len=@res;
	foreach $i (0..$len - 2) { $dimer_counts{$res[$i]}{$res[$i+1]}++; }

	foreach $let1 ("A","C","G","T") {
		foreach $let2 ("A","C","G","T") {
			if ($dimer_counts{$let1}{$let2}) {
				#print "#($let1,$let2)=",$dimer_counts{$let1}{$let2},"\n";
				#print "$let1$let2",$dimer_counts{$let1}{$let2},"\t";
				push(@results,$dimer_counts{$let1}{$let2});
				$dinucs_sum += $dimer_counts{$let1}{$let2};
				#print "foo\n";
				#print "$let1$let2\t$dimer_counts{$let1}{$let2}\n";
				#print $dinucs_sum, "\n";
			} else {
				push(@results,0);
				#print "*$let1$let2\n";

			}
		}
	}
	$DNA = "";
	#print "CG is ", $dimer_counts{"C"}{"G"}, "\n";
	#print "dinucs sum is ", $dinucs_sum, "\n";
	if ($dinucs_sum > 0) {
		$freqCG = $percgc = sprintf("%.4f", $dimer_counts{"C"}{"G"} / $dinucs_sum * 100);
	} else {
		$freqCG = 0;
	}
	print join("\t", @results), "\t$freqCG\n";



	
########

	
	
}
######################
#@@@@@@@@@@@@@@@@@@@@@
######################


