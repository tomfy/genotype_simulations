#!/usr/bin/perl -w
use strict;
use warnings;
use List::Util qw(min max sum);
# use Graphics::GnuplotIF qw(GnuplotIF);
use Getopt::Long;
use Math::GSL::RNG  qw( :all );
use Math::GSL::Randist  qw( :all );
# use GenotypeSim qw ( :all );
use Genotype;

# fasta2paupnexus.pl : read in fasta format genotypes
# construct Genotype objects and output nexus
# files for use as paup input.

{                               # main

   # defaults:

   my $rng_type = $gsl_rng_mt19937;
   my $rng_seed = undef;
   my $sample_size = undef;
   my $nexus_character_set = '012'; # anything else gives aA, aa, AA.
   my $nexus_filename = undef;

   GetOptions(
              'rng_type=s' => \$rng_type, # mt -> 'Mersenne Twister'
              'seed=i' => \$rng_seed,
              'sample_size=i' => \$sample_size,
              'character_set=s' => \$nexus_character_set,
              'output_nexus_filename=s' => \$nexus_filename,
             );

   # ############## set up the random number generator ##########################################
   if ($rng_type eq 'sys') {
      $rng_type = $gsl_rng_default;
   } elsif ($rng_type eq 'mt') {
      $rng_type = $gsl_rng_mt19937;
   } elsif ($rng_type eq 'lux') {
      $rng_type = $gsl_rng_ranlxd2;
   }
   print STDERR "RNG type: $rng_type \n";
   my $the_rng = (defined $rng_seed)? Math::GSL::RNG->new($rng_type, $rng_seed) : Math::GSL::RNG->new($rng_type) ; # RNG
   # ############################################################################################


   my @genotype_objects = (); 
   while (<>) {
      if (/^>(\S+)\s+(\S+)\s+(\S+)/) {
         my ($id, $gen, $pedigree) = ($1, $2, $3);
         my $gstring = <>;
         my $gobj = Genotype->new_from_012string($the_rng, $gstring, $gen, $id, $pedigree);
         push @genotype_objects, $gobj;
      }
   }
   if ( defined $sample_size and ($sample_size < scalar @genotype_objects) ) {
       @genotype_objects = $the_rng->shuffle(@genotype_objects);
      @genotype_objects = @genotype_objects[0..$sample_size-1];
   }
   print genotypes_as_nexus(\@genotype_objects);


}


################################################################

sub genotypes_as_nexus{
   my $genotype_objects = shift; # array ref of Genotype objects.

   my $ntaxa = scalar @$genotype_objects;
   my $taxa_block = "begin taxa; \n";;
   $taxa_block .= " dimensions ntax=$ntaxa \n";
   $taxa_block .= " taxlabels \n";
   my $characters_block = "begin characters; \n" . " dimensions nchar=" . "NCHARS; \n";
   $characters_block .= " format missing=? gap=- matchchar=. datatype=dna; \n";
   $characters_block .= " options gapmode=missing; \n";
   $characters_block .= " matrix \n";

   my $nchars = undef;
   for my $gobj (@$genotype_objects) {
      my ($gen, $id, $parents, $gtstring) = (
                                             $gobj->get_generation(),
                                             $gobj->get_id(),
                                             $gobj->get_parents(),
                                             $gobj->genotype_012string()
                                            );
      if (defined $nchars) {
         my $l = length $gtstring;
         die "sequence length inconsistency: $nchars  $l \n" if($l != $nchars);
      } else {
         $nchars = length $gtstring;
      }
      $taxa_block .= sprintf("%-10s \n", $id);
      #      $characters_block .= sprintf("%-10s  %s\n", $id, $gobj->genotype_string());
      $characters_block .= sprintf("%-10s  %s\n", $id, $gtstring);
   }
   $characters_block =~ s/NCHARS/$nchars/;
   $taxa_block .= "; \n" . "end;\n";
   $characters_block .= "; \n" . "end;\n";
   return "#NEXUS \n" . $taxa_block . $characters_block;
}


