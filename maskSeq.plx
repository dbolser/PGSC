#!/usr/bin/perl -w

use strict;
use Getopt::Long;

use Bio::SeqIO;
use Bio::FeatureIO;



## Set options

my $verbose = 0;
my $seq_file;
my $gff_file;

GetOptions
  (
   "verbose" => \$verbose,
   "seq=s" => \$seq_file,
   "gff=s" => \$gff_file,
  )
  or die "failed to parse command line options\n";

die "fail $gff_file : $!\n"
  unless -s $gff_file;



## Set up the BioPerl objects

my $seq_reader =
  Bio::SeqIO->new( -file => $seq_file,
		   -format => 'fasta'
		 );

my $seq_writer =
  Bio::SeqIO->new( -fh => \*STDOUT,
		   -format => 'fasta',
		   -width => 100
		 );

my $gff_reader =
  Bio::FeatureIO->new( -file => $gff_file,
		       -format => 'GFF',
		     );

#warn $seq_reader->width, "\n"; exit;



## Run

my (%repeats, $c);

while ( my $feature = $gff_reader->next_feature() ) {
  if($verbose>1){
    print
      join("\t", #$feature,
	   $feature->seq_id,
	   $feature->type->name,
	   $feature->start,
	   $feature->end,
	  ), "\n";
  }
  
  if($feature->type->name eq 'repeat_region'){
    $c++;
    push @{$repeats{ $feature->seq_id }},
      [$feature->start,
       $feature->end];
  }
  
  # Debugging
  #last if $c > 100;
}

warn "read $c repeat_region features for ",
  scalar keys(%repeats), " sequences\n";



##

while(my $seq = $seq_reader->next_seq){
  my $id = $seq->id;
  my $sequence = $seq->seq;
  
  print $id, "\n"
    if $verbose > 0;
  
  print length($sequence), "\n"
    if $verbose > 0;
  
  for my $region (@{$repeats{ $id }}){
    my ($start, $end) = @$region;
    print "$start\t$end\n"
      if $verbose > 1;
    
    substr($sequence, $start, $end - $start, 'X' x ($end - $start));
  }
  
  print length($sequence), "\n"
    if $verbose > 0;
  
  $seq->seq($sequence);
  
  $seq_writer->write_seq($seq);
  
  # Debugging;
  #last;
}

warn "OK\n";

