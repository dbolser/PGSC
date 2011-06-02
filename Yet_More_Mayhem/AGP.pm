package AGP;

use Data::Dumper;

use Moose;
use MooseX::FileAttribute;

use Bio::Location::Simple;
use Bio::Coordinate::Pair;
use Bio::Coordinate::Collection;


has 'agp' =>
  ( is => 'ro',
    isa => 'Bio::Coordinate::Collection',
    
    ## Build it
    default => sub{ Bio::Coordinate::Collection->new() },
    
    ## Delegate
    handles => [ qw( map add_mapper swap ) ],
  );


## just so that we can pass a single file at create time
has_file 'agp_file' =>
  ( must_exist => 1,
    init_arg => 'file',
    trigger => \&load_agp_file,
  );


sub load_agp_file {
  my $self = shift;
  my $file = shift;
  
  confess "pass an AGP file plz\n"
    unless -s $file;
  
  open AGP, '<', $file
    or die "failed to open file '$file' : $!\n";
  
  while(<AGP>){
    next if /^#/;
    next if /^\s*$/;
    chomp;
    
    my ($obj, $obj_beg, $obj_end,
	$comp_idx, $comp_type, $comp_id,
	$comp_beg, $comp_end, $comp_ori) = split/\t/;
    
    next unless $comp_type eq 'W';
    
    $comp_ori = +1 if $comp_ori eq 'unknown';
    
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
  
  return 1;
}

sub components {
  my $self = shift;
  
  my @components;
  
  ## $self->agp is a Bio::Coordinate::Colleciton, composed of several
  ## 'mappers'. Each 'mapper' is a Bio::Coordinate::Pair.
  push @components,
    $_->in->seq_id for $self->agp->mappers;
  
  return @components;
}

1;
