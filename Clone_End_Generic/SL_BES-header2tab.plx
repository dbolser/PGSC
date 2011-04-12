#!/usr/bin/perl -w

use strict;

## We use this to assess duplicate IDs and sequences
use Digest::MD5 qw(md5_base64);

## We use this to parse the fasta file
use Bio::SeqIO;

die "pass bacends_combined_screened_and_trimmed.seq\n"
  unless @ARGV;

my $obj =
  Bio::SeqIO->new( -file => '<'. $ARGV[0],
		   -format => 'fasta' );



my $i = 0;

while ( my $seq = $obj->next_seq() ){
  $i++;
  
  my $seq_id = $seq->id;
  
  # Perhaps not the most efficient regexp, but hopfully clear...
  die "fail : '$seq_id'\n" unless
    $seq_id =~ /^(LE_HBa)  (\d{4}[A-P]\d{2})_(  T7| SP6)_(\d+)(?:_(\d+))?$/x ||
    $seq_id =~ /^(SL_MboI) (\d{4}[A-P]\d{2})_(  T7| SP6)_(\d+)(?:_(\d+))?$/x ||
    $seq_id =~ /^(SL_EcoRI)(\d{4}[A-P]\d{2})_(  T7| SP6)_(\d+)           $/x ||
    $seq_id =~ /^(SL_FOS)  (\d{4}[A-P]\d{2})_(pIBF|pIBR)_(\d+)           $/x ||
    $seq_id =~ /^(SL_s)    (\d{4}[A-P]\d{2})_(  T7| SP6)_(\d+)           $/x ||
    $seq_id =~ /^(SL_MT)   (\d{4}[A-P]\d{2})_(  T7| SP6)_(\d+)           $/x ||
    $seq_id =~ /^(LpenCOS) (\d{4}[A-P]\d{2})_(M13F|M13R)_(\d+)           $/x ||
    $seq_id =~ /^(LpenCOS) (\d{4}[A-P]\d{2})_(  T7| SP6)_(\d+)           $/x ||
    $seq_id =~ /^(LpenBAC) (\d{4}[A-P]\d{2})_(  T7| SP6)_(\d+)(?:_(\d+))?$/x;
  
  my ($libraryLongName, $bacName, $primer, $readSerial, $readSerialX) =
    ($1, $2, $3, $4, $5);
  
  my $libraryName = $libraryLongName;
  
  ## Do I look bothered by efficency?
  $libraryName  =~ s/^LE_HBa  $/HBa/x;
  $libraryName  =~ s/^SL_MboI $/SLm/x;
  $libraryName  =~ s/^SL_EcoRI$/SLe/x;
  $libraryName  =~ s/^SL_FOS  $/SLf/x;
  $libraryName  =~ s/^SL_s    $/SLs/x;
  $libraryName  =~ s/^SL_MT   $/SLt/x;
  $libraryName  =~ s/^LpenCOS $/LPc/x;
  $libraryName  =~ s/^LpenBAC $/LPb/x;
  
  die "failed to recover library : '$seq_id'\n"
    if $libraryName eq $libraryLongName;
  
  
  
  $primer =~ s/^  T7$/F/x;
  $primer =~ s/^ SP6$/R/x;
  
  $primer =~ s/^pIBF$/F/;
  $primer =~ s/^pIBR$/R/;
  
  $primer =~ s/^M13F$/F/;
  $primer =~ s/^M13R$/R/;
  
  die "failed to get read direction : '$seq_id'\n"
    if $primer ne 'F' and $primer ne 'R';
  
  
  
  print
    join("\t",
	 ## The original seq_id is needed to match the blast or SSAHA2
	 ## query results
	 $seq_id,
	 
	 ## This name and the direction is all we need
	 $libraryName.$bacName,
	 $primer,
	 
	 #$libraryName,
	 #$libraryLongName,
	 #$readSerial,
	 
	 ## Couple of 'stats'
	 $seq->length,
	 md5_base64($seq->seq),
	 
	), "\n";
}

warn "processed $i sequences\n";
warn "OK\n"
