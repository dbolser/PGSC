#!/usr/bin/perl -w

use strict;
use Digest::MD5 qw(md5_base64);

use Bio::SeqIO;

my $obj =
  Bio::SeqIO->new( -file => $ARGV[0],
		   -format => 'fasta' );

my $i = 0;

while ( my $seq = $obj->next_seq() ){
  $i++;
  
  my $seqid = $seq->id;
  
  # Perhaps not the most efficient regexp, but hopfully clear...
  die "fail : '$seqid'\n" unless
    $seqid =~ /^(LE_HBa)  (\d{4}[A-P]\d{2})_(  T7| SP6)_(\d+)(?:_(\d+))?$/x ||
    $seqid =~ /^(SL_MboI) (\d{4}[A-P]\d{2})_(  T7| SP6)_(\d+)(?:_(\d+))?$/x ||
    $seqid =~ /^(SL_EcoRI)(\d{4}[A-P]\d{2})_(  T7| SP6)_(\d+)           $/x ||
    $seqid =~ /^(SL_FOS)  (\d{4}[A-P]\d{2})_(pIBF|pIBR)_(\d+)           $/x ||
    $seqid =~ /^(SL_s)    (\d{4}[A-P]\d{2})_(  T7| SP6)_(\d+)           $/x ||
    $seqid =~ /^(SL_MT)   (\d{4}[A-P]\d{2})_(  T7| SP6)_(\d+)           $/x ||
    $seqid =~ /^(LpenCOS) (\d{4}[A-P]\d{2})_(M13F|M13R)_(\d+)           $/x ||
    $seqid =~ /^(LpenCOS) (\d{4}[A-P]\d{2})_(  T7| SP6)_(\d+)           $/x ||
    $seqid =~ /^(LpenBAC) (\d{4}[A-P]\d{2})_(  T7| SP6)_(\d+)(?:_(\d+))?$/x;
  
  my ($libraryLongName, $bacName, $primer, $readSerial, $readSerialX) =
    ($1, $2, $3, $4, $5);
  
  my $libraryName = $libraryLongName;
  
  $libraryName  =~ s/LE_HBa  /HBa/x;
  $libraryName  =~ s/SL_EcoRI/SLe/x;
  $libraryName  =~ s/SL_MboI /SLm/x;
  $libraryName  =~ s/SL_FOS  /SLf/x;
  $libraryName  =~ s/LpenCOS /LPc/x;
  $libraryName  =~ s/LpenBAC /LPb/x;
  
  $primer =~ s/  T7/F/x;
  $primer =~ s/ SP6/R/x;
  
  $primer =~ s/pIBF/F/;
  $primer =~ s/pIBR/R/;
  
  $primer =~ s/M13F/F/;
  $primer =~ s/M13R/R/;
  
  print
    join("\t",
	 $seqid,
	 $libraryName.$bacName,
	 #$libraryName,
	 #$libraryLongName,
	 #$readSerial,
	 $primer,
	 $seq->length,
	 #$seq->seq,
	 md5_base64($seq->seq),
	), "\n";
}

warn "processed $i sequences\n";
warn "OK\n"
