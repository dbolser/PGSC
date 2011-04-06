#!/usr/bin/perl -w

## A script to check the 'depth' of alignment at every position on
## every reference(hit) and query sequence in a Blast or a SSAHA2
## tabular report. Alignments mostly within regions above the 'depth
## threshold' are discarded.



use strict;
use Getopt::Long;

## DEFAULT: Regions of alignments with a 'depth' this deep (or deeper)
## are deemed 'repetative'.
my $depthThreshold = 5;

## DEFAULT: The fraction of the alignment that must have a depth BELOW
## the depth threshold to be acceptable.
my $coverageThreshold = 0.8;

## DEFAULT: Do we measure alignment depth by position on the QUERY,
## HIT or BOTH?
my $sequenceFlag = "BOTH";



my $file;
my $format = 'ssaha';

my $verbose = 0;

## Parse the command line for options
GetOptions(
	   "infile|i=s"   => \$file,
	   "format|f=s"   => \$format,
	   "depth|d=i"    => \$depthThreshold,
	   "coverage|c=f" => \$coverageThreshold,
	   "type|t=s"     => \$sequenceFlag,
	   "verbose|v+"   => \$verbose,
	  )
  or die "could not parse command line for options\n";

$file = $ARGV[0] unless $file;

die "give me an infile!\n"
  unless $file;

die "no such format!\n"
  unless $format =~ /ssaha|blasttable/i;

die "depth should be a positive integer greater than 1\n"
  unless $depthThreshold > 1;

die "coverage should be a fraction in the range 0.10 to 0.90\n"
  unless $coverageThreshold > 0 && $coverageThreshold < 1;

die "invalid type!\n"
  unless $sequenceFlag =~ /query|hit|both/i;



warn "using the following settings:\n";
warn "infile : $file\n";
warn "format : $format\n";
warn "depth threshold : $depthThreshold\n";
warn "coverage threshold : $coverageThreshold\n";
warn "on sequence : $sequenceFlag\n";
warn "\n";





## Go time

## Calculate the positional 'depth' accross this sequence

my %depth;

warn "parsing alignment report\n";

open DATA, '<', $file
  or die "failed to open input file : $!\n";

while(<DATA>){
  my ($queryName, $hitName,
      $qs, $qe, $hs, $he);
  
  ## blast table
  if($format =~ /blasttable/i){
    ($queryName, $hitName, undef, undef, undef, undef,
     $qs, $qe, $hs, $he, undef, undef,
    ) = split/\t/;
  }
  
  ## ssaha2 output
  if($format =~ /ssaha/i){
    (undef, $queryName, $hitName,
     $qs, $qe, $hs, $he, undef, undef, undef, undef,
    ) = split/\s+/;
  }
  
  ## Checking
  unless(defined($queryName) && defined($hitName) &&
	 defined($qs) && defined($qe) &&
	 defined($hs) && defined($he)){
    warn $_;
    die "is this format corect : $format?\n";
  }
  
  ## Note: sorting is required for the following code
  ($qs, $qe) = sort {$a<=>$b} ($qs, $qe);
  ($hs, $he) = sort {$a<=>$b} ($hs, $he);
  
  # Explicitly create a 'depth vector' for each sequence (query, hit
  # or both)
  
  if ($sequenceFlag =~ /query|both/i){
    for my $i ( $qs .. $qe ){
      $depth{$queryName}[$i]++
    }
  }
  if ($sequenceFlag =~ /hit|both/i){
    for my $i ( $hs .. $he ){
      $depth{$hitName  }[$i]++
    }
  }
}

warn "created depth vector for ", scalar( keys %depth ), " sequencs\n";
warn "where sequence type is : $sequenceFlag\n";





warn "applying the filter\n";

open DATA, '<', $file
  or die "failed to open input file : $!\n";

while(<DATA>){
  my ($queryName, $hitName,
      $qs, $qe, $hs, $he);
  
  ## blast table
  if($format =~ /blasttable/i){
    ($queryName, $hitName, undef, undef, undef, undef,
     $qs, $qe, $hs, $he, undef, undef,
    ) = split/\t/;
  }
  
  ## ssaha2 output
  if($format =~ /ssaha/i){
    (undef, $queryName, $hitName,
     $qs, $qe, $hs, $he, undef, undef, undef, undef,
    ) = split/\s+/;
  }
  
  ## Note: sorting is required for the following code
  ($hs, $he) = sort {$a<=>$b} ($hs, $he);
  ($qs, $qe) = sort {$a<=>$b} ($qs, $qe);
  
  
  
  # Filter alignments by depth
  
  if ($sequenceFlag =~ /query|both/i){
    my $ql = 0;
    my $qt = 0;
    for my $i ($qs .. $qe){
      $ql++;
      $qt++ if $depth{$queryName}[$i] >= $depthThreshold;
    }
    #next if $qt>0; # Simulate $coverateThrshold = 0;
    next if $qt / $ql > ( 1 - $coverageThreshold );
  }
  
  if ($sequenceFlag =~ /hit|both/i){
    my $hl = 0;
    my $ht = 0;
    for my $i ($hs .. $he){
      $hl++;
      $ht++ if $depth{$hitName  }[$i] >= $depthThreshold;
    }
    #next if $ht>0; # Simulate $coverateThrshold = 0;
    next if $ht / $hl > ( 1 - $coverageThreshold );
  }
  
  # We passed the filter...
  print "$_";
}

warn "OK\n";
