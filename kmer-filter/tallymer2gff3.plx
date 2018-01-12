#!/usr/bin/perl -w

use strict;
use Getopt::Long;

## Define some default paramters (configurable on the cli)

## Size of the kmer we're using to mine repeats
my $k;

## Minimum length to consider a region a repeat. Regions smaller than
## this will be ignored. Set to k + 1 by default.
my $minimum_len;

## Maximum gap between two high 'kmer coverage' regions before they
## are considered separate features
my $maximum_gap;

## Minimum 'kmer coverage' to consider repetative. This is safe to
## leave at zero.
my $minimum_cov = 0;

## Sequence file (used to simply get sequene names). This is expected
## to contain sequence names in fasta (header) format.
my $seq_file = '';



GetOptions
  (
   "k=i" => \$k,
   "min=i" => \$minimum_len,
   "gap=i" => \$maximum_gap,
   "cov=i" => \$minimum_cov,
   "seq-file=s" => \$seq_file,
  )
  or die "failed to parse command line\n";



my @seq_names =
    seq_names_from_file ( $seq_file );

die "pass a (non-empty) sequence (header) file to extract sequence names\n"
    unless @seq_names;



## Define some sensible defaults for the above parameters, if not
## given on the cli

## This should match what was run by tallymer (and is usually reported
## in the header of the tmer file)
$k = 21
  unless defined $k;

$minimum_len = $k
  unless defined $minimum_len;

## Changing this from $k + 1 to $k * 3 reduced the number of repeat
## features by a factor of about 15 but only increased the number of
## masked bases by about 5%, which seemed to be a good tradeoff...
$maximum_gap = $k * 3
  unless defined $maximum_gap;





## Parse the output of the tallymer search

my $r_len = -1; # Rolling length
my $r_cov = -1; # Rolling coverage

my $p_seq = -1; # Previous seq id
my $p_pos = -1; # Previous position

while(<>){
  # ignore comments and empty lines
  next if /^#/; next if /^$/;
  
  chomp;
  my ($seq, $pos, $cov) = split;
  
  ## By default, we trust the 'minocc' coverage, however, we can go
  ## larger here (makes me why bother setting minocc)...
  next if $cov < $minimum_cov;

  ## are we beginning a new sequence?
  if($seq != $p_seq){
    
    ## we should dump the previous region (unless we just began...)
    &print_region( $p_seq, $p_pos, $r_len, $r_cov, $k, $minimum_len )
        if $p_seq > -1;
    
    ## and reset our counters
    $r_len = 0;
    $r_cov = $cov;
    $p_seq = $seq;
    $p_pos = $pos;
    next;
  }
  
  ## if not, are we continuing a region?
  if( $maximum_gap > $pos - $p_pos - $k - 1 ){
    $r_len += $pos - $p_pos;
    $r_cov += $cov;
    $p_pos  = $pos;
    next;
  }
  
  ## if not, we must be leaving a region...
  
  ## we should dump the previous region
  &print_region( $p_seq, $p_pos, $r_len, $r_cov, $k, $minimum_len );
  
  ## and reset our counters
  $r_len = 0;
  $r_cov = $cov;
  $p_pos = $pos;
}



## And finally...
my ($num_repeats, $num_bases) =
    &print_region( $p_seq, $p_pos, $r_len, $r_cov, $k, $minimum_len );

warn "found ", $num_repeats || 0, " repeats ".
  "covering ", $num_bases || 0, " bases\n";

warn "OK\n";





sub seq_names_from_file {
    my $seq_file = shift;
    my @seq_names;

    warn "reading sequence names from seq-file '$seq_file'\n";

    open NAMES, '<', $seq_file
        or warn "failed to open file '$seq_file'\n";

    while(<NAMES>){
        next unless /^>/;
        chomp;
        push @seq_names, substr($_, 1);
    }

    warn "read ", scalar(@seq_names), " sequence names\n";

    return @seq_names;
}





sub print_region () {
  my $seq = shift;
  my $pos = shift;
  my $len = shift;
  my $cov = shift;

  my $k = shift;
  my $minimum_length = shift;

  our $num_repeats;
  our $num_bases;

  return ($num_repeats, $num_bases)
      if $len + $k < $minimum_len;

  ## This is just the first thing I could think of, it's the average
  ## kmer coverage across the repeat region
  my $score = $cov / ($len || 1);

  print
    join("\t",
	 $seq_names[$p_seq],
	 'tallymer',
	 'repeat_region',
	 $pos + 1 - $len,
	 $pos + $k,
	 sprintf("%d", $score),
	 '+',
	 '.',
	 join(';',
	      'ID='. sprintf("%08d", $num_repeats++),
	      'dbxref=SO:0000657',
              "length=". ($len + $k),
	     ),
	), "\n";

  $num_bases += $len + $k;

  return ($num_repeats, $num_bases);
}
