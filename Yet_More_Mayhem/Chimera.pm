package Chimera;

use Data::Dumper;

use Moose;
use MooseX::FileAttribute;

use Bio::Location::Simple;
use Bio::Coordinate::Pair;
use Bio::Coordinate::Collection;


has 'chimera' =>
  (
   is => 'ro',
   isa => 'Bio::Coordinate::Collection',
   
   ## Build it
   default => sub{ Bio::Coordinate::Collection->new() },
   
   ## Delegate
   handles => [ qw( map add_mapper swap ) ],
  );


## just so that we can pass a single file at create time
has_file 'chimera_gff_file' =>
  ( must_exist => 1,
    init_arg => 'file',
    trigger => \&load_chimera_gff,
  );





=head1 NAME

Chimera - Perl module for mapping over the chimeric superscaffolds

=head1 SYNOPSIS

  use Chimera;

  $feature = FooFoo::old_to_new_scaff($feature);

=head1 DESCRIPTION

This module helps manage mapping features from old superscaffolds onto
the new superscaffolds.

=head2 Methods

=over 4

=item * old_to_new_scaff

Does the mapping from old to new scaffolds.

=back

=head1 AUTHOR

Dan B (dan.bolser@gmail.com)

=cut







sub load_chimera_gff {
  my $self = shift;
  my $file = shift;
  
  confess "pass a GFF file plz\n"
    unless -s $file;
  
  open C, '<', $file
    or die "failed to open file '$file' : $!\n";
  
  while(<C>){
    chomp;
    next if /^#/;
    next if /^\s*$/;
    
    ## Grab the GFF fields we need
    my ($seq_id, undef, undef, $st, $en,
	undef, $strand, undef, $attrs) = split/\t/;
    
    ## Parse the attributes field
    my %attrs = split(/=|;/, $attrs);
    
    ## Sanity checks for the chimera GFF
    die "$_\n" unless $strand eq '+';
    die "$_\n" unless $attrs{'ID'} =~ /^PGSC0003DMB\d{9}$/;
    die "$_\n" unless $attrs{'Note'} =~ /^Chr\. \d\d$/;
    
    
    
    ## NB, the coordinates (extent) of the 'new' sequence are not the
    ## same as the coordinates of the new sequence in the 'old'
    ## sequence i.e. "old sequence 50 to 100" becomes "new sequence 1
    ## to 51".
    
    my $sx = 1;
    my $ex = $en - $st + 1;
    
    
    
    ## Create a Bio::Coordinate::Pair (map) to store the mapping
    ## between the old and the new sequence
    
    my $old_scaff = Bio::Location::Simple->
      new( -seq_id => $seq_id,
	   -start  => $st,
	   -end    => $en,
	   -strand => +1,
	 );
    #print Dumper $old_scaff;
    
    my $new_scaff = Bio::Location::Simple->
      new( -seq_id => $attrs{'ID'},
	   -start  => $sx,
	   -end    => $ex,
	   -strand => +1,
	 );
    #print Dumper $new_scaff;
    
    
    my $map = Bio::Coordinate::Pair->
      new( -in  => $old_scaff,
	   -out => $new_scaff
	 );
    #print dumper $map;
    
    $self->add_mapper( $map );
  }
  
  return 1;
}

1;

