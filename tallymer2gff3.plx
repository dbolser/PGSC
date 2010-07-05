#!/usr/bin/perl -w

use strict;
use Getopt::Long;

## Define some paramters
my $k = 19;
my $merge_thresh = $k;
my $min_len_thresh = 50;

my $seq_file;



## Get options
GetOptions
  (
   "k=i"    => \$k,
   "min=i"  => \$min_len_thresh,
   "file=s" => \$seq_file,
  )
  or die "failed to parse command line\n";

die "pass a sequence file\n"
  unless -s $seq_file;



## Read in sequence names
warn "reading in sequence names from $seq_file\n";
my @seq_names;

open NAMES, '<', $seq_file
  or die "failed to open $seq_file\n";

while(<NAMES>){
  next unless /^>/; chomp;
  push @seq_names, substr($_, 1);
}

warn "read ", scalar(@seq_names), "\n";
warn "OK\n";



## Parse the output of the tallymer search

my ($p_seq, $p_pos, $i) = (-1, -1, -1);

while(<>){
  # ignore commets
  next if /^#/;
  
  chomp;
  my ($seq, $pos, $x, undef) = split;
  
  ## are we beginning a new sequence?
  if($seq != $p_seq){
    #warn $seq, "\n";
    
    # we should dump the current region
    &print_region( $p_seq, $p_pos, $i )
      if $i >= $min_len_thresh;
    
    # and reset our counters
    $i = 0;
    $p_pos = $pos;
    $p_seq = $seq;
    next;
  }
  
  ## if not, are we continuing a region?
  if($pos - $p_pos <= $merge_thresh){
    $i += $pos - $p_pos;
    $p_pos = $pos;
    next;
  }
  
  ## if not, we must be leaving a region...
  
  # we should dump the current region
  &print_region( $p_seq, $p_pos, $i )
    if $i >= $min_len_thresh;
  
  # and reset our counters
  $i = 0;
  $p_pos = $pos;
}

warn "OK\n";



sub print_region () {
  $p_seq = shift;
  $p_pos = shift;
  $i = shift;
  
  our $j;
  
  print
    join("\t",
	 $seq_names[$p_seq],
	 'tallymer',
	 'repeat_region',
	 $p_pos - $i,
	 $p_pos +  0,
	 '1',
	 '+',
	 '.',
	 join(';',
	      'ID='. sprintf("%08d", $j++),
	      'dbxref=SO:0000657'
	     ),
	), "\n";
  
  return 1;
}
