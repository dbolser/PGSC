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
    
    die "fail 1 : '$_'\n"
      unless /^>gi\|(\d+)\|gb\|(GS\d+)\.(1)\| (LuSp\d+[A-P]\d+)_(\d+)_(TV|TP) /;
    
    $gi = join("\t", $1, $2, $3, $4, $5, $6);
  }
  else{
    chomp;
    $len += length($_)
  }
}

print "$gi\t$len\n";

warn "OK\n";
