#!/usr/bin/perl -w

use strict;
use Bio::SeqIO;

my $verbose = 0;



## read in gff

warn "reading GFF\n";

my %gff;

open GFF, '<', 'repeat.gff3'
  or die "fail\n";

while(<GFF>){
  my ($seqid, undef, undef, $start, $end,
      undef, undef, undef, $attrs) = split;
  
  push @{$gff{$seqid}}, [$start, $end, $attrs];
}

warn "OK\n";



## Do the fasta

my $seqio = Bio::SeqIO->
  new( -file => 'Assembly/PGSC0003DM/PGSC0003DMS.fa',
       -format => 'fasta' )
  or die "double fail\n";


while(my $sobj = $seqio->next_seq){
  my $seqid = $sobj->id;
  
  unless(defined($gff{$seqid})){
    warn "no repeats for $seqid\n";
    next;
  }
  
  my @repeats = @{$gff{$seqid}};
  
  my $seq = $sobj->seq;
  
  for(@repeats){
    my ($start, $end, $attrs) = @$_;
    
    warn join("\t", $start, $end, $attrs), "\n"
      if $verbose > 0;
    
    my %attrs = split(/=|;/, $attrs);
    
    print ">$seqid-". $attrs{"ID"}.
      "/$start-$end (". ($end-$start+1). ")\n";
    
    print substr($seq, $start, $end-$start+1), "\n";
  }
  
  #exit;
}

warn "OK\n";
