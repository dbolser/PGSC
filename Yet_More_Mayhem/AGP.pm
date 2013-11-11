package AGP;

use Data::Dumper;

use Moose;
use MooseX::FileAttribute;

use Bio::Location::Simple;
use Bio::Coordinate::Pair;
use Bio::Coordinate::Collection;

## Could use Bio::AGP::LowLevel to parse the AGP!

has 'mapper' =>
  (
   is => 'ro',
   isa => 'Bio::Coordinate::Collection',

   ## Delegate
   handles => [ qw( map add_mapper swap ) ],

   ## Build it
   default => sub{ Bio::Coordinate::Collection->new() },
  );

## Just so that we can pass a file at create time...
has_file 'agp_file' =>
  (
   must_exist => 1,
   init_arg => 'file',
   trigger => \&load_agp_file,
  );


sub load_agp_file {
  my $self = shift;
  my $file = shift;

  confess "pass an AGP file plz\n"
    unless -s $file;

  open C, '<', $file
    or die "failed to open file '$file' : $!\n";

  while(<C>){
    ## Ignore comments or blank lines
    next if /^#/;
    next if /^\s*$/;
    chomp;
    
    ## Parse the AGP
    my ($obj, $obj_beg, $obj_end,
	$comp_idx, $comp_type, $comp_id,
	$comp_beg, $comp_end, $comp_ori) = split/\t/;
    
    next unless $comp_type eq 'W';
    
    $comp_ori = +1 if $comp_ori eq 'unknown';

    ## Create a Bio::Coordinate::Pair (map) to store the mapping
    ## between the feature and its reference sequence

    my $scaff = Bio::Location::Simple->
      new( -seq_id => $comp_id,
	   -start  => $comp_beg,
	   -end    => $comp_end,
	   -strand => +1,
	 );
    #print Dumper $scaff;
    
    my $scaff_on_chr = Bio::Location::Simple->
      new( -seq_id => $obj,
	   -start  => $obj_beg,
	   -end    => $obj_end,
	   -strand => $comp_ori,
	 );
    #print Dumper $scaff_on_chr;
    
    my $map = Bio::Coordinate::Pair->
      new( -in  => $scaff,
	   -out => $scaff_on_chr,
	 );
    #print Dumper $map;
    
    $self->add_mapper( $map );
  }

  return $self->components;
}

sub components {
  my $self = shift;
  
  my @components;
  
  ## $self->mapper is a Bio::Coordinate::Colleciton, composed of
  ## several 'mappers'. Each 'mapper' is a Bio::Coordinate::Pair.
  push @components,
    $_->in->seq_id for $self->mapper->mappers;
  
  return @components;
}

1;
