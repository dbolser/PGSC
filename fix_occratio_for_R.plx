#!/usr/bin/perl -w

use strict;

my %types =
  ("# distribution of unique mers"                                                   => 'uniq',
   "# distribution of non unique mers (counting each non unique mer only once)"      => 'nun1',
   "# distribution of non unique mers (counting each non unique mer more than once)" => 'nunN',
   "# distribution of all mers (counting each non unique mer only once)"             => 'all1',
   "# distribution of all mers (counting each non unique mer more than once)"        => 'allN',
  );

warn "there are ", scalar keys %types, " types\n";



my %res;
my ($mink, $maxk) = ( 1e9, -1 );


my $type;

while(<>){
  chomp;
  
  if(/^#/){
    $type = $types{$_} or die;
    #warn "'$type'\n";
    next;
  }
  
  my ($k, $a, $b) = split(/\s/, $_);
  
  $a = "$a\t$b" if defined($b);
  
  $mink = $k if $k < $mink;
  $maxk = $k if $k > $maxk;
  
  $res{$type}{$k} = $a;
}


for my $k ( ($mink .. $maxk) ){
  print
    join("\t", $k,
	 ( defined($res{'uniq'}) ? $res{'uniq'}{$k} : (0, 0) ),
	 ( defined($res{'nun1'}) ? $res{'nun1'}{$k} : (0, 0) ),
	 ( defined($res{'nunN'}) ? $res{'nunN'}{$k} : (0, 0) ),
	 ( defined($res{'all1'}) ? $res{'all1'}{$k} :     0  ),
	 ( defined($res{'allN'}) ? $res{'allN'}{$k} :     0  ),
	), "\n";
}


