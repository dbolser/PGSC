#!/usr/bin/perl -w

use strict;


# -output unique nonunique nonuniquemulti total relative

my %types =
  ("# distribution of unique mers"                                                   => 'uniq',
   "# distribution of non unique mers (counting each non unique mer only once)"      => 'nun1',
   "# distribution of non unique mers (counting each non unique mer more than once)" => 'nunN',
   "# distribution of all mers (counting each non unique mer only once)"             => 'all1',
   "# distribution of all mers (counting each non unique mer more than once)"        => 'allN',
  );

warn "there are ", scalar keys %types, " types\n";



my %res;
#my ($mink, $maxk) = ( 1e9, -1 );

my $type;

while(<>){
  chomp;
  
  if(/^# /){
      if(/^# distribution of /){
          $type = $types{$_} or die;
          warn "'$type'\n";
      }
      next;
  }
  
  my ($k, $a, $b) = split(/\s/, $_);
  
  $a = "$a\t$b" if defined($b);
  
  #$mink = $k if $k < $mink;
  #$maxk = $k if $k > $maxk;
  
  $res{$k}{$type} = $a;
}


for my $k ( sort keys %res ){
  print
    join("\t", $k,
	 ( exists($res{$k}{'uniq'}) ? $res{$k}{'uniq'} : (0, 0) ),
	 ( exists($res{$k}{'nun1'}) ? $res{$k}{'nun1'} : (0, 0) ),
	 ( exists($res{$k}{'nunN'}) ? $res{$k}{'nunN'} : (0, 0) ),
	 ( exists($res{$k}{'all1'}) ? $res{$k}{'all1'} :     0  ),
	 ( exists($res{$k}{'allN'}) ? $res{$k}{'allN'} :     0  ),
	), "\n";
}
