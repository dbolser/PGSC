#!/usr/bin/perl -w

use strict;
use Getopt::Long;

use Bio::SeqIO;
use Bio::FeatureIO;



## Set options

my $verbose = 0;

my $feature_file;
my $sequence_file;

my $feature_format  = 'GFF';
my $sequence_format = 'fasta';

my $feature_to_mask = 'repeat_region';

my $seq_mask_character = 'X';



GetOptions
  (
   "verbose"    => \$verbose,
   
   "feature_file|f=s"  => \$feature_file,
   "sequence_file|s=s" => \$sequence_file,
   
   "feature_format|ff=s"  => \$feature_format,
   "sequence_format|sf=s" => \$sequence_format,
   
   "feature_to_mask|m=s"  => \$feature_to_mask,
   
   "seq_mask_character|c=s" => \$seq_mask_character,
  )
  or die "failed to parse command line options\n";



## Check options

## A value should be passed
die usage() unless
  $sequence_file &&
  $feature_file;

## The files should exist
die "problem with feature file '$feature_file' : $!\n"
  unless -s $feature_file;
die "problem with sequence file '$sequence_file' : $!\n"
  unless -s $sequence_file;

## The formats should be valid
$feature_format = lc($feature_format);
die "ERROR: feature format '$feature_format' is not supported!\n\n"
  unless eval( "require Bio::FeatureIO::$feature_format" );

$sequence_format = lc($sequence_format);
die "ERROR: sequence format '$sequence_format' is not supported!\n\n"
  unless eval( "require Bio::SeqIO::$sequence_format" );

## Erm...
die "1: what are you trying to do?\n"
  unless $feature_to_mask;

die "2: what are you trying to do?\n"
  unless length($seq_mask_character) == 1;



=head1 NAME

 bp_repeat_mask_sequence.pl - mask sequence features

=head1 DESCRIPTION

 Takes an input sequence file and a feature file, and returns the
 sequence with 'repeat_region' features masked out (replaced with
 X's). This is useful for downstream processing of the sequence file.

 The masked sequence is written to STDOUT

=head1 USAGE

 bp_repeat_mask_sequence.pl <options>

 Options:

    -f
    --feature_file        The file from which the sequence features will
                          be read (for subsequent masking).

    -s
    --sequence_file       The sequence file (to be  masked).

    --ff
    --feature_format      The format of the feature file
                          (the default is GFF).

    --sf
    --sequence_format     The format of the sequence file
                          (the default is fasta).

    -m
    --feature_to_mask     The type of feature to mask
                          (the default is 'repeat_region').

    -c
    --seq_mask_character  The 'mask' character to use in the sequence.
                          (the default is 'X').

    -v
    --verbose             Generate some debugging output

=cut





## Set up the BioPerl objects

my $gff_reader =
  Bio::FeatureIO->new( -file => $feature_file,
		       -format => $feature_format
		     );

my $seq_reader =
  Bio::SeqIO->new( -file => $sequence_file,
		   -format => $sequence_format,
		 );

my $seq_writer =
  Bio::SeqIO->new( -fh => \*STDOUT,
		   -format => $sequence_format,
		 );



## Run

warn "hashing features to mask\n";

my (%repeats, $c);

while ( my $feature = $gff_reader->next_feature() ) {
  if($verbose>0){
    print
      join("\t", #$feature,
	   $feature->seq_id,
	   $feature->type->name,
	   $feature->start,
	   $feature->end,
	  ), "\n";
  }
  
  if($feature->type->name eq $feature_to_mask){
    $c++;
    push @{$repeats{ $feature->seq_id }},
      [$feature->start,
       $feature->end];
  }
}

warn "read $c '$feature_to_mask' features for ",
  scalar keys(%repeats), " sequences\n";



warn "masking sequences\n";

while(my $seq = $seq_reader->next_seq){
  my $id = $seq->id;
  my $sequence = $seq->seq;
  
  print $id, "\n"
    if $verbose > 0;
  
  ## Do the masking
  for my $region (@{$repeats{ $id }}){
    my ($start, $end) = @$region;
    print "$start\t$end\n"
      if $verbose > 1;
    
    substr($sequence, $start, $end - $start,
	   $seq_mask_character x ($end - $start)
	  );
  }
  
  $seq->seq($sequence);
  
  $seq_writer->write_seq($seq);
}

warn "done\n";



sub usage{
  `perldoc -T ./$0`
}
