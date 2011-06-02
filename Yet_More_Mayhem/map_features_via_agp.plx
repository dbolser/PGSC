#!/usr/bin/perl -w

## This script assumes we connect to a specific feature database (see
## below), and that the AGPs 'objects' are sequences with features in
## that database.

## IMPORTANT NOTE, the features are on 'full' superscaffolds, we
## assume that the AGP deals with the 'split' superscaffolds (although
## it should work fine if they are not split). For that reason, we use
## the 'chimera' module to map between the full and the split
## sequences.

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

my $feature_type;

my $verbose = 0;

GetOptions( 'verbose+'       => \$verbose,
	    'feature_type=s' => \$feature_type,
	    'agp_file=s'     => \$agp_file,
	  )
  or die "failure to communcate\n";

die "please pass a feature type to map! (-f)\n"
  unless defined($feature_type);

die "please pass an agp file to map over (-a)\n"
  unless -s $agp_file;



## GO TIME

## Initalise the AGP object
my $agp = AGP->new();

$agp->load_agp_file($agp_file);

warn "the AGP ('$agp_file') has ",
  scalar $agp->components, " components\n\n";



##
## Connect to the feature database
##

warn "connecting to feature db\n";

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



## Query features from the database

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
  #  $scaf2 eq 'PGSC0003DMB000000029';
  
  ## Get a list of all the features of the given type on this scaffold
  my @features = $db->
    features( -seqid => $scaf2,
	      -type  => $feature_type,
	    );
  
  warn "collected ", scalar @features, " features from $scaf2\n"
    if $verbose > 0;
  
  
  
  ## Map the AGP component coordinates onto the AGP object coordinates
  
  foreach my $feature (@features){
    
    ## Debugging
    warn Dumper $feature
      if $verbose > 1;
    
    my $mapped_feature = $agp->map( $feature );
    
    ## Debugging
    warn Dumper $mapped_feature
      if $verbose > 1;
    
    ## Because of chimerisms, some features can map outside the given
    ## superscaffold coordinates. The mapping object 'conveniently'
    ## reports this...
    next if $mapped_feature->purge_gaps;
    
    
    
    ## Rebuild the attribues string. I'm sure this is harder than it should be
    my @attrs;
    for ($feature->get_all_tags){
      next if /load_id/;
      push @attrs,
	"$_=". join(",", $feature->get_tag_values($_));
    }
    
    my $id = ($feature->get_tag_values('load_id'))[0];
    
    my $attrs =
      join(";",
	   "ID=$id",
	   "Name=". ($feature->display_name || $id),
	   @attrs
	  );
    
    
    
    ## Dump the GFF
    
    print
      join("\t",
	   $mapped_feature-> seq_id,
	   $feature->        source_tag,
	   $feature->        primary_tag,
	   $mapped_feature-> start,
	   $mapped_feature-> end,
	   $feature->        score || 0,
	   $mapped_feature-> strand,
	   $feature->        phase || '.',
	   $attrs,
	  ), "\n";
  }
}

warn "OK\n\n";
