#!/usr/bin/perl -w

use strict;

=head1 DESCRIPTION

This script produces GFF from BES alignments (SSAHA2 format). It uses
a specific regexp to match one end to another, and uses a BES sequence
specific data table to construct the GFF. The script removes very
short or very long pairs, so they will never be seen!

=cut

use Getopt::Long;

## Set the source field, if desired
my $source_tag = 'dundee';

## Set to 1 to enable debugging output
my $verbose = 0;

GetOptions( "source=s"  => \$source_tag,
	    "v"         => \$verbose,
	  )
  or die "failed to communicate!\n";



## We expect to be passed the (filtered) SSAHA2 alignment file(s)

die "pass ssaha file\n"
  unless @ARGV;





## Read in the BES specif sequence header table so that we can match
## the sequence identifier in the alignment file to the more
## informative sequence header (the header is parsed into a table
## elsewhere).

## At this stage, we also construct the BES pairs for a BAC (clone)

warn "\nreading sequence data table\n";

my $seq_data_file =
  "Data/RH_BAC_ends_curated.20070603.header.tab";

open SEQ, '<', $seq_data_file
  or die "failed to open $seq_data_file : $!\n";

my (%sequence_pairs, $seq_count);

while(<SEQ>){
  my ($gi, $gb, $gb_ver, $clone, $clone_id,
      $direction, $length) = split/\t/;
  
  $seq_count++;
  
  ## No direction should be duplicated!
  die "$clone\n" if
    exists($sequence_pairs{$clone}{$direction});
  
  ## Store the 'putative pair'
  $sequence_pairs{$clone}{$direction} = "$clone\_$direction";
}

warn "got $seq_count BES (and ",
  scalar keys %sequence_pairs, " BACs)\n";

warn "\nNote, not all BACs have BES for both directions!\n";





## Read in the alignment results

warn "\nparsing ssaha results\n";

my (%hits, $hit_counter);

while(<>){
  my ($score, $query_name, $hit_name, $qs, $qe, $hs, $he,
      $strand, $align_len, $ident, $query_len) = split/\s+/;
  
  $hit_counter++;
  
  ## Store all the hits for this query together
  push
    @{$hits{$query_name}},
      [$score, $query_name, $hit_name, $qs, $qe, $hs, $he,
       $strand, $align_len, $ident, $query_len];
}

warn "got $hit_counter hits for ",
  scalar keys %hits, " BES\n";







## Produce the GFF

warn "\ndumping GFF\n";

print "##gff-version 3\n";
print "##species RH\n";
print "##genome-build BGI 3\n";

my %fail;

foreach my $clone (keys %sequence_pairs){
  
  ## Really debugging
  #next unless $clone eq 'LuSp113E12';
  
  ## Debugging
  print "C:\t$clone\n"
    if $verbose > 0;
  
  ## Get the hits for both sequences
  my $query_id1 = $sequence_pairs{$clone}{'F'} || 'seq missing';
  my $query_id2 = $sequence_pairs{$clone}{'R'} || 'seq missing';
  
  my @hits1 = exists($hits{$query_id1}) ? @{$hits{$query_id1}} : ();
  my @hits2 = exists($hits{$query_id2}) ? @{$hits{$query_id2}} : ();
  
  ## Debugging
  print "G:\t$query_id1\t", scalar(@hits1), "\n"
    if $verbose > 0;
  print "G:\t$query_id2\t", scalar(@hits2), "\n"
    if $verbose > 0;
  
  
  
  ## Failure accounting

  ## Hits for both ends? If not, we give up now (we don't want to
  ## 'pollute' the GFF with these kinds of problem)
  
  if(scalar @hits1 == 0){
    if($query_id1 eq 'seq missing'){
      $fail{'missing seq'}++;
    }
    else{
      $fail{'no align'}++;
    }
  }
  if(scalar @hits2 == 0){
    if($query_id2 eq 'seq missing'){
      $fail{'missing seq'}++;
    }
    else{
      $fail{'no align'}++;
    }
  }
  next unless @hits1 && @hits2;
  
  
  ## Try to pair the hits (best hits first)
  
  ## NOTE: if the best hits for each end don't form a 'good pair', we
  ## give up on the clone. About 3% of the clones could be 'recovered'
  ## by hunting for a lower scoring 'good pair', however, their
  ## quality is dubious, and they are complex to handle! We skip them
  ## here.
  
  my $hit1 = (sort hits_by_score @hits1)[0]; # Forward
  my $hit2 = (sort hits_by_score @hits2)[0]; # Reverse
  
  my ($sc1, $qn1, $hn1, $qs1, $qe1,
      $hs1, $he1, $st1, $al1, $nt1, $ql1) = @$hit1;
  my ($sc2, $qn2, $hn2, $qs2, $qe2,
      $hs2, $he2, $st2, $al2, $nt2, $ql2) = @$hit2;
  
  ## Debugging
  print join("\t", '1:', @$hit1), "\n"
    if $verbose > 0;
  print join("\t", '2:', @$hit2), "\n"
    if $verbose > 0;
  
  
  
  ## Check they hit the same scaffold. If not they must span a pair of
  ## scaffolds, and are handled differently (below).
  
  if ($hn1 eq $hn2){
    
    ## Check orientation 1/2
    if ($st1 eq 'F'){
      if ($st2 ne 'C'){
	## Debugging
	print "W:\tincorrect orientation (parallel)\n"
	  if $verbose > 0;
	$fail{'bad ori, para'}++;
	next;
      }
      
      ## OK. Check their distance
      my $bes_distance = $he2 - $hs1 + 1;
      print "I:\tdistance is $bes_distance\n"
	if $verbose > 0;
      
      if ($bes_distance < 5_000){
	if ($bes_distance < 0){
	  ## Debugging
	  print "W:\tincorrect orientation, (away)\n"
	    if $verbose > 0;
	  $fail{'bad ori, away'}++;
	}
	else{
	  ## Debugging
	  print "W:\ttoo short\n"
	    if $verbose > 0;
	  $fail{'too short'}++;
	}
	next;
      }
      
      if ($bes_distance > 500_000){
	## Debugging
	print "W:\ttoo long!\n"
	  if $verbose > 0;
	$fail{'too long'}++;
	next;
      }
      
      ## OK! Print three lines of GFF
      print join("\t", $hn1, $source_tag, 'BAC', $hs1, $he2,  '.', '+', '.', "ID=$clone;Name=$clone"), "\n";
      print join("\t", $hn1, $source_tag, 'BES', $hs1, $he1, $sc1, '+', '.', "ID=$query_id1;Parent=$clone"), "\n";
      print join("\t", $hn1, $source_tag, 'BES', $hs2, $he2, $sc2, '-', '.', "ID=$query_id2;Parent=$clone"), "\n";
      
      ## DONE!
      $fail{'GOOD PAIR!'}++;
      next;
    }
    
    
    
    ## Check orientation 2/2
    if ($st1 eq 'C'){
      if ($st2 ne 'F'){
	## Debugging
	print "W:\tincorrect orientation (parallel)\n"
	  if $verbose > 0;
	$fail{'bad ori, para'}++;
	next;
      }
      
      ## OK. Check their distance
      my $bes_distance = $he1 - $hs2 + 1;
      print "I:\tdistance is $bes_distance\n"
	if $verbose > 0;
      
      if ($bes_distance < 5_000){
	if ($bes_distance < 0){
	  ## Debugging
	  print "W:\tincorrect orientation, (away)\n"
	    if $verbose > 0;
	  $fail{'bad ori, away'}++;
	}
	else{
	  ## Debugging
	  print "W:\ttoo short\n"
	    if $verbose > 0;
	  $fail{'too short'}++;
	}
	next;
      }
      
      if ($bes_distance > 500_000){
	## Debugging
	print "W:\ttoo long!\n"
	  if $verbose > 0;
	$fail{'too long'}++;
	next;
      }
      
      ## OK! Print three lines of GFF
      print join("\t", $hn1, $source_tag, 'BAC', $hs2, $he1,  '.', '-', '.', "ID=$clone;Name=$clone"), "\n";
      print join("\t", $hn1, $source_tag, 'BES', $hs1, $he1, $sc1, '-', '.', "ID=$query_id1;Parent=$clone"), "\n";
      print join("\t", $hn1, $source_tag, 'BES', $hs2, $he2, $sc2, '+', '.', "ID=$query_id2;Parent=$clone"), "\n";
      
      ## DONE!
      $fail{'GOOD PAIR!'}++;
      next;
    }
  }
  
  
  
  ## Nothing paired, falling through to single end printer
  
  ## Debugging
  print "W:\tnot on the same scaffold\n"
    if $verbose > 0;
  
  ## The only remaining unpaired ends are those that span
  ## superscaffolds!
  
  ## Debugging
  print join("\t", '1:', @$hit1), "\n" if $verbose > 0;
  print join("\t", '2:', @$hit2), "\n" if $verbose > 0;
  
  ## ONE (Forward)
  if(0){}
  elsif($st1 eq 'F'){
    print join("\t", $hn1, "$source_tag link", 'BAC', $hs1, $he1+1000,  '.', '+', '.', "ID=$clone\_F;Name=$clone\_F;Note=Other end matches $hn2"), "\n";
    print join("\t", $hn1, "$source_tag link", 'BES', $hs1, $he1,      $sc1, '+', '.', "Parent=$clone\_F"), "\n";
  }
  elsif($st1 eq 'C'){
    print join("\t", $hn1, "$source_tag link", 'BAC', $hs1-1000, $he1,  '.', '-', '.', "ID=$clone\_F;Name=$clone\_F;Note=Other end matches $hn2"), "\n";
    print join("\t", $hn1, "$source_tag link", 'BES', $hs1,      $he1, $sc1, '-', '.', "Parent=$clone\_F"), "\n";
  }
  else{die}
  
  ## TWO (Reverse)
  if(0){}
  elsif($st2 eq 'F'){
    print join("\t", $hn2, "$source_tag link", 'BAC', $hs2, $he2+1000,  '.', '-', '.', "ID=$clone\_R;Name=$clone\_R;Note=Other end matches $hn1"), "\n";
    print join("\t", $hn2, "$source_tag link", 'BES', $hs2, $he2,      $sc2, '+', '.', "Parent=$clone\_R"), "\n";
  }
  elsif($st2 eq 'C'){
    print join("\t", $hn2, "$source_tag link", 'BAC', $hs2-1000, $he2,  '.', '+', '.', "ID=$clone\_R;Name=$clone\_R;Note=Other end matches $hn1"), "\n";
    print join("\t", $hn2, "$source_tag link", 'BES', $hs2, $he2,      $sc2, '-', '.', "Parent=$clone\_R"), "\n";
  }
  else{die}
}



## output failures
warn "$_\t$fail{$_}\n" for keys %fail;

warn "OK\n";





sub hits_by_score {
  # score
  $b->[0] <=> $a->[0] ||
    # ident
    $b->[9] <=> $a->[9];
}

