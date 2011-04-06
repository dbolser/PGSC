#!/usr/bin/perl -w

## Convert superscaffold feature GFF into pseudomolecule feature GFF:
## LIFT OVER.

## In this special case, we also pair up newly mated PE.

use strict;

use Data::Dumper;

use Getopt::Long;

## Query feature data
use Bio::DB::SeqFeature::Store;

## Module that provides mapping access to an AGP file
use AGP;

## Module that provides mapping access to the superscaffold 'chimera'
## GFF. Note, the location of the GFF file is encoded within the
## module.
use Chimera;



## HANDLE COMMAND LINE OPTIONS

my $agp_file;
my $mates_file;

my $feature_type;
my $library_id;

my $verbose = 0;

GetOptions( 'verbose+'       => \$verbose,
	    'library_id=s'   => \$library_id,
	    'feature_type=s' => \$feature_type,
	    'agp_file=s'     => \$agp_file,
	    'mates_file=s'   => \$mates_file,
	  )
  or die "failure to communicate\n";



## SANITY CHECK COMMAND LINE OPTIONS

die "please pass an agp file\n"
  unless defined($agp_file) && -s $agp_file;

die "please pass a mates file\n"
  unless defined($mates_file) && -s $mates_file;

die "please pass a feature type\n"
  unless defined($feature_type);

die "please pass a library_id\n"
  unless defined($library_id);


my ($feat_type, $feat_source) =
  split(/:/, $feature_type);

die "feature type ('$feature_type') needs to be in the format 'type:source'\n"
  unless defined($feat_type) && defined($feat_source);

die "feature type 'source' ('$feat_source') needs to be a 'link' type\n"
  unless $feat_source =~ /^(.+) link$/;
$feat_source = $1;





##
## BEGIN PROCESSING
##

## Feed in the appropriate AGP file to create an AGP object
my $agp = AGP->new( file => $agp_file );

warn "The AGP file ('$agp_file') has ",
  scalar $agp->components, " components\n\n";

my $chimera = Chimera->
  new( file => 'GFF/chimeric_superscaffold_v3_split_report.gff' );



## Process the mates file...
my (%mates, $curr_lib);

open M, '<', $mates_file
  or die "failed to open file '$mates_file' : $!\n";

while(<M>){
  chomp;
  if(/^library/){
    my (undef, $lib, $min, $max, $lib_regexp) = split/\t/;
    $curr_lib = $lib;
    
    $mates{$curr_lib}{min} = $min;
    $mates{$curr_lib}{max} = $max;
    $mates{$curr_lib}{lib_regexp} = qr/$lib_regexp/;
  }
  if(/^pair/){
    my (undef, $forward_regexp, $reverse_regexp) = split/\t/;
    
    $mates{$curr_lib}{forward_regexp} = qr/$forward_regexp/;
    $mates{$curr_lib}{reverse_regexp} = qr/$reverse_regexp/;
  }
}

warn "Read ", scalar keys %mates, " mates libs\n";
warn "OK\n\n";

my $library = $mates{$library_id}
  or die "lib: '$library_id' not found!\n";



## Connect to the feature database to query superscaffold features

warn "Connecting to superscaffold feature db\n";

my $db = Bio::DB::SeqFeature::Store->
  new( -adaptor => 'DBI::mysql',
       -dsn     => 'www-potato:mysql.compbio.dundee.ac.uk',
       -user    => 'www-potato',
       -pass    => 'abc123',
       -verbose => $verbose,
       
       ## We use namespaces to keep things nice and clean
       -namespace => 'gb_pot_qa',
     );

warn "OK\n\n";





## Begin analysis

warn "Collecting PE features ('$feature_type') for AGP components\n";

my (%clones, %hack);
my $i = 0;

for my $scaff ($agp->components){
  
  ## SUPERSCAFFOLD CHIMERA HACKERY. Revert to the old 'full'
  ## superscaffold id for feature lookup (so that we can use the new
  ## 'split' superscaffold id's in our AGP without worrying)...
  
  my $scaf2 = $scaff;
  $scaf2 =~ s/PGSC0003DMB10\d0/PGSC0003DMB0000/;
  
  ## Sanity check
  die "failure of sanity : '$scaf2'\n"
    unless $scaf2 =~ /^PGSC0003DMB0000/;
  
  ## Prevent double lookup for split superscaffolds (this is needed
  ## because we don't hash by scaff)
  next if $hack{$scaf2}++;
  
  ## Debugging
  #next unless
  #  $scaf2 eq 'PGSC0003DMB000000177' ||
  #  $scaf2 eq 'PGSC0003DMB000000268';
  
  
  
  ## Get a list of all the currently unpaired PEs features in this
  ## superscaffold using the appropriate GFF feature 'type' code.
  ## These feature types all end in the keyword 'link'.
  
  my @clone_ends = $db->
    features(
	     -seqid => $scaf2,
	     -type  => $feature_type,
	    );
  
  warn "collected ", scalar @clone_ends, " clone ends from $scaf2\n"
    if $verbose > 0;
  
  
  
  ## For each clone_end feature on this superscaffold
  
  for my $clone_end (@clone_ends){
    warn Dumper $clone_end if $verbose > 2;
    
    ## Lookup the parent id from the clone_end (the parent clone)
    my $clone_id = ($clone_end->get_tag_values('parent_id'))[0];
    
    ## Check it and grab the 'base' clone id (the 'link' type clones
    ## have a directional clone id (forward or reverse).
    die "failed to recognize clone id '$clone_id'\n"
      unless $clone_id =~ $library->{lib_regexp};
    my $base_clone_id = $1;
    warn "$base_clone_id\n" if $verbose > 1;
    
    
    
    ## Now we grab the actual clone feature to get its 'source',
    ## 'type', and 'Notes' attributes (used to rebuild the GFF)
    my @clone = $db->
      get_feature_by_name( -name => $clone_id );
    
    ## One parent only please!
    warn Dumper scalar @clone if @clone != 1 && $verbose > 0;
    die "found more than one parent for '$clone_id'!\n", if @clone != 1;
    
    warn Dumper @clone if $verbose > 2;

    my $p_source = $clone[0]->source;
    my $p_type   = $clone[0]->type;
    
    ## Dear Jesus! (More code is less good).
    ($p_type, undef) = split(/:/, $p_type);
    
    ## Add the required clone attributes to the clone end feature
    $clone_end->add_tag_value( p_source => $p_source );
    $clone_end->add_tag_value( p_type   => $p_type   );
    
    die "FUCK!\n" if !$clone[0] ->has_tag('Note');
    die "FUCK!\n" if  $clone_end->has_tag('Note');
    
    $clone_end->add_tag_value( 'Note' => ($clone[0]->get_tag_values('Note'))[0] );
    
    warn Dumper $clone_end if $verbose > 2;
    
    
    
    ## Get the feature positions from 'full' superscaffolds back into
    ## 'split' superscaffold coordinate space (or do nothing if this
    ## isn't a 'split' superscaffold). See $scaf2 above.
    
    my $clone_end_mapped = $chimera->map( $clone_end );
    
    warn Dumper $clone_end_mapped if $verbose > 2;
    
    $clone_end->seq_id( $clone_end_mapped->seq_id );
    $clone_end->start ( $clone_end_mapped->start  );
    $clone_end->end   ( $clone_end_mapped->end    );
    $clone_end->strand( $clone_end_mapped->strand );
    
    ## Reject the features from parts of the 'full' superscaffold that
    ## are not in our AGP (built using 'split' superscaffolds).
    next unless $scaff eq $clone_end->seq_id;
    
    ## Store each clone end feature by it's parent clone
    $i++;
    push @{ $clones{$base_clone_id} }, $clone_end;
  }
}

warn "found ", scalar keys %clones, " clones and $i clone ends\n";
warn "OK\n\n";







warn "producing GFF...\n";

for my $clone_id (keys %clones){
  
  ## DEBUGGING
  #next unless $clone_id eq 'GBSKQZK02J130R';
  
  my @clone_ends = @{ $clones{$clone_id} };
  
  warn join("\t", $clone_id, scalar(@clone_ends)), "\n"
    if $verbose > 1;
  
  
  
  ##
  ## Paired within the AGP object?
  ##
  
  if (@clone_ends == 2){
    my $ce1 = $clone_ends[0];
    my $ce2 = $clone_ends[1];
    
    warn Dumper $ce1 if $verbose > 1;
    warn Dumper $ce2 if $verbose > 1;
    
    
    
    ## Map the scaffold coordinates onto the AGP object
    my $ce1_mapped = $agp->map( $ce1 );
    my $ce2_mapped = $agp->map( $ce2 );
    
    ## Because of chimerisms, some clone ends can map outside the
    ## given superscaffold coordinates. The mapping object
    ## 'conveniently' reports this...
    next if $ce1_mapped->purge_gaps;
    next if $ce2_mapped->purge_gaps;
    
    
    
    ## Did they pair within the same AGP object?
    
    if( $ce1_mapped->seq_id ne
	$ce2_mapped->seq_id ){
      warn "spanning ultrascaffolds\n" if $verbose > 1;
      
      ## DID NOT PAIR: PRINT EACH CLONE END SEPARATELY
      
      ## Forward or reverse read?
      my $par1 = ($ce1->get_tag_values('parent_id'))[0];
      my $dir1;
      
      if(    $par1 =~ $library->{forward_regexp} ){ $dir1 = 'F' }
      elsif( $par1 =~ $library->{reverse_regexp} ){ $dir1 = 'R' }
      else{die "Waaaaaaaaa\n"}
      
      print
	join("\t",
	     $ce1_mapped->seq_id, "$feat_source link",
	     ($ce1->get_tag_values('p_type'))[0],
	     $ce1_mapped->start + ($ce1_mapped->strand > 0 ? 0 : -10000),
	     $ce1_mapped->end   + ($ce1_mapped->strand > 0 ? +10000 : 0),
	     0,
	     ## Note, a forward read in the + strand indicates a + clone
	     ## Hwvr, a reverse read in the + strand indicates a - clone
	     $dir1 eq 'F' ?
	      $ce1_mapped->strand > 0 ? '+' : '-':
	      $ce1_mapped->strand > 0 ? '-' : '+',
	     '.',
	     join(';',
		  "ID=$par1", "Name=$clone_id", "Note=Other end matches ". $ce2->seq_id. " ". $ce2_mapped->seq_id,
		 ),
	    ), "\n";
      
      print
	join("\t",
	     $ce1_mapped->seq_id, "$feat_source link", $feat_type,
	     $ce1_mapped->start,
	     $ce1_mapped->end,
	     $ce1->score,
	     $ce1_mapped->strand > 0 ? '+' : '-',
	     '.',
	     join(';',
		  "Parent=$par1", "Note=On ". $ce1->seq_id,
		 ),
	    ), "\n";
      
      ## Forward or reverse read?
      my $par2 = ($ce2->get_tag_values('parent_id'))[0];
      my $dir2;

      if(    $par2 =~ $library->{forward_regexp} ){ $dir2 = 'F' }
      elsif( $par2 =~ $library->{reverse_regexp} ){ $dir2 = 'R' }
      else{die "Waaaaaaaaa\n"}
      
      print
	join("\t",
	     $ce2_mapped->seq_id, "$feat_source link",
	     ($ce2->get_tag_values('p_type'))[0],
	     $ce2_mapped->start + ($ce2_mapped->strand > 0 ? 0 : -10000),
	     $ce2_mapped->end   + ($ce2_mapped->strand > 0 ? +10000 : 0),
	     0,
	     ## Note, a forward read in the + strand indicates a + clone
	     ## Hwvr, a reverse read in the + strand indicates a - clone
	     $dir2 eq 'F' ?
	       $ce2_mapped->strand > 0 ? '+' : '-':
	       $ce2_mapped->strand > 0 ? '-' : '+',
	     '.',
	     join(';',
		  "ID=$par2", "Name=$clone_id", "Note=Other end matches ". $ce1->seq_id. " ". $ce1_mapped->seq_id,
		 ),
	    ), "\n";
      
      print
	join("\t",
	     $ce2_mapped->seq_id, "$feat_source link", $feat_type,
	     $ce2_mapped->start,
	     $ce2_mapped->end,
	     $ce2->score,
	     $ce2_mapped->strand > 0 ? '+' : '-',
	     '.',
	     join(';',
		  "Parent=$par2", "Note=On ". $ce2->seq_id,
		 ),
	    ), "\n";
    }
    
    
    
    ##
    ## PAIRED!
    ##
    
    else{
      
      ## Calculate clone size (sidestepping pair 'correctness')
      my ($st, $en) =
	range($ce1_mapped->start, $ce1_mapped->end,
	      $ce2_mapped->start, $ce2_mapped->end);
      my $insert_size = $en - $st;
      
      
      
      ## Get into the dirty details
      
      ## Forward or reverse read?
      my $par1 = ($ce1->get_tag_values('parent_id'))[0];
      my $dir1;
      
      if(    $par1 =~ $library->{forward_regexp} ){ $dir1 = 'F' }
      elsif( $par1 =~ $library->{reverse_regexp} ){ $dir1 = 'R' }
      else{die "Waaaaaaaaa\n"}
      
      ## Note, a forward read in the + strand indicates a + clone
      ## Hwvr, a reverse read in the + strand indicates a - clone
      my $str1 =
	$dir1 eq 'F' ?
	  $ce1_mapped->strand > 0 ? '+' : '-':
	  $ce1_mapped->strand > 0 ? '-' : '+';
	      
      ## Forward or reverse read?
      my $par2 = ($ce2->get_tag_values('parent_id'))[0];
      my $dir2;
      
      if(    $par2 =~ $library->{forward_regexp} ){ $dir2 = 'F' }
      elsif( $par2 =~ $library->{reverse_regexp} ){ $dir2 = 'R' }
      else{die "Waaaaaaaaa\n"}
      
      ## Note, a forward read in the + strand indicates a + clone
      ## Hwvr, a reverse read in the + strand indicates a - clone
      my $str2 =
	$dir2 eq 'F' ?
	  $ce1_mapped->strand > 0 ? '+' : '-':
	  $ce1_mapped->strand > 0 ? '-' : '+';
      
      ## Sanity check
      die "WTSF?\n" if $dir1 eq $dir2;
      
      
      
      my $ori = 'b'; ## b for bad
      my $str = 0;
      
      if($dir1 eq 'F' && $dir2 eq 'R'){
	if($ce1_mapped->start < $ce2_mapped->start){
	  if($str1 eq '+' && $str2 eq '-'){
	    $ori = 'g';
	    $str = '+';
	  }
	}
	else{
	  if($str2 eq '+' && $str1 eq '-'){
	    $ori = 'g';
	    $str = '-';
	  }
	}
      }
      else{
	if($ce2_mapped->start < $ce1_mapped->start){
	  if($str2 eq '+' && $str1 eq '-'){
	    $ori = 'g';
	    $str = '+';
	  }
	}
	else{
	  if($str1 eq '+' && $str2 eq '-'){
	    $ori = 'g';
	    $str = '-';
	  }
	}
      }
      
      my $happy;
      $happy = 'h';
      $happy = 's' if $insert_size < $library->{min};
      $happy = 'l' if $insert_size > $library->{max};
      $happy = 'f' if $insert_size > $library->{max} * 4;
      
      
      
      print
	join("\t",
	     $ce1_mapped->seq_id, "$feat_source $happy$ori",
	     ($ce1->get_tag_values('p_type'))[0],
	     $st, $en, 0, $str, '.',
	     join(';',
		  "ID=$clone_id", "Name=$clone_id",
		 ),
	    ), "\n";
      
      print
	join("\t",
	     $ce1_mapped->seq_id, "$feat_source", $feat_type,
	     $ce1_mapped->start,
	     $ce1_mapped->end,
	     $ce1->score, $str1, '.',
	     join(';',
		  "ID=$par1", "Parent=$clone_id", "Note=On ". $ce1->seq_id,
		 ),
	    ), "\n";
      
      print
	join("\t",
	     $ce2_mapped->seq_id, "$feat_source", $feat_type,
	     $ce2_mapped->start,
	     $ce2_mapped->end,
	     $ce2->score, $str2, '.',
	     join(';',
		  "ID=$par2", "Parent=$clone_id", "Note=On ". $ce2->seq_id,
		 ),
	    ), "\n";
      #exit;
    }
  }
  
  
  
  ##
  ## clone_end to nowhere! (PE not found anywhere in AGP)
  ##
  
  elsif (@clone_ends == 1){
    my $ce0 = $clone_ends[0];
    
    warn Dumper $ce0 if $verbose > 1;
    
    ## Map the scaffold coordinates onto the AGP object
    my $ce0_mapped = $agp->map( $ce0 );
    
    ## Because of chimerisms, some clone ends can map outside the
    ## given superscaffold coordinates. The mapping object
    ## 'conveniently' reports this...
    next if $ce0_mapped->purge_gaps;
    
    ## Forward or reverse read?
    my $par0 = ($ce0->get_tag_values('parent_id'))[0];
    my $dir0;
    
    if(    $par0 =~ $library->{forward_regexp} ){ $dir0 = 'F' }
    elsif( $par0 =~ $library->{reverse_regexp} ){ $dir0 = 'R' }
    else{die "Waaaaaaaaa\n"}
    
    
    
    print
      join("\t",
	   $ce0_mapped->seq_id,
	   "$feat_source link",
	   ($ce0->get_tag_values('p_type'))[0],
	   $ce0_mapped->start + ($ce0_mapped->strand > 0 ? 0 : -10000),
	   $ce0_mapped->end   + ($ce0_mapped->strand > 0 ? +10000 : 0),
	   0,
	   ## Note, a forward read in the + strand indicates a + clone
	   ## Hwvr, a reverse read in the + strand indicates a - clone
	   $dir0 eq 'F' ?
	     $ce0_mapped->strand > 0 ? '+' : '-':
	     $ce0_mapped->strand > 0 ? '-' : '+',
	   '.',
	   join(';',
		"ID=$par0", "Name=$clone_id", #"Note=". $ce0->desc,
	       ),
	  ), "\n";
    
    print
      join("\t",
	   $ce0_mapped->seq_id, "$feat_source link", $feat_type,
	   $ce0_mapped->start,
	   $ce0_mapped->end,
	   $ce0->score,
	   $ce0_mapped->strand > 0 ? '+' : '-',
	   '.',
	   join(';',
		"Parent=$par0", "Note=On ". $ce0->seq_id,
	       ),
	  ), "\n";
    #exit;
  }
  
  
  
  ##
  ## Unesplainable!
  ##
  
  else{
    warn "This clone '$clone_id' has ", scalar(@clone_ends), " clone ends!\n";
    warn "continuing\n";
  }
}

warn "OK\n\n";



sub range {
  my $min = $_[+0];
  my $max = $_[-1];
  
  for(@_[1..($#_-1)]){
    $min = $_ if $_ < $min;
    $max = $_ if $_ > $max;
  }
  
  return ($min, $max);
}


