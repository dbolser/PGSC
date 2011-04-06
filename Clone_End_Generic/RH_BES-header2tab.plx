#!/usr/bin/perl -w

die "pass fasta\n"
  unless @ARGV;

my ($gi, $len);

while(<>){
  if(/^>/){
    if (defined($gi)){
      print "$gi\t$len\n";
      $len = 0;
    }
    chomp;
    
    die "fail : '$_'\n"
      unless /^>(RH\d{3}[A-P]\d{2})\_(F|R) gi\|(\d+)\|gb\|(E[IR]\d+)\.(1)\|\4\|(\S+)$/;
    
    $gi = join("\t", $3, $4, $5, $1, $6, $2);
  }
  else{
    chomp;
    $len += length($_)
  }
}

print "$gi\t$len\n";

warn "OK\n";
