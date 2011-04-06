#!/usr/bin/perl -w

die "pass fasta\n"
  unless @ARGV;

my ($gi, $len, %sanity_check);

while(<>){
  if(/^>/){
    if (defined($gi)){
      ## Print information after finding the next sequence header
      print "$gi\t$len\n";
      $len = 0;
    }
    chomp;
    
    die "fail 1 : '$_'\n"
      unless /^>gi\|(\d+)\|gb\|(FI\d+)\.(\d)\|\2 (\d+)_(TR|TF) /;
    
    ## GI / GB / GB_VER / CLONE / DIRECTION
    $gi = join("\t", $1, $2, $3, $4, $5);
    die if $sanity_check{$2}++;
  }
  else{
    chomp;
    $len += length($_)
  }
}

## Don't forget to print the last piece of info...
print "$gi\t$len\n";

warn "OK\n";
