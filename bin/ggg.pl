#!/usr/bin/perl
use strict;
use warnings;
use List::Util qw(min max sum);
# use Graphics::GnuplotIF qw(GnuplotIF);
use Getopt::Long;


use File::Basename 'dirname';
use Cwd 'abs_path';
my ( $bindir, $libdir );
BEGIN {     # this has to go in Begin block so happens at compile time
   $bindir =
     dirname( abs_path(__FILE__) ) ; # the directory containing this script
   $libdir = $bindir . '/../lib';
   $libdir = abs_path($libdir); # collapses the bin/../lib to just lib
}
use lib $libdir;


use Math::GSL::RNG  qw( :all );
use GenotypeSimulation qw ( :all );
use SimulatedGenotype;

# ggg.pl : "generate generations of genotypes"

# the idea here is to start with a set of unrelated genotypes
# and let them pair up randomly and produce offspring,
# for some specified number of generations.

{                               # main


   # defaults:
   my $n_snps = 10;             # number of snps in each sample.
   my $n_generations = 3;
   my $pop = 50;                # population size
   my $mafspec = 'delta,0.25';  # minor allele frequency 
   my $max_pedigree_depth = 2;
   my $rng_type = $gsl_rng_mt19937;
   my $rng_seed = 135791;
   my $n_gens_out = undef;
   my $sample_size = undef;
   my $fasta_character_set = '012'; # anything else gives aA, aa, AA.
   my $fasta_filename = undef;
   my $keep_prob = undef;

   GetOptions(
              'nsnps=i' => \$n_snps, # 
              'ngens|ngenerations=i' => \$n_generations,
              'population=i' => \$pop,
              'maf|minor_allele_frequence=s' => \$mafspec, # 0 < $maf <= 0.5
              'depth_limit=i' => \$max_pedigree_depth,
              'rng_type=s' => \$rng_type, # mt -> 'Mersenne Twister'
              'seed=i' => \$rng_seed,
              'nout=i' => \$n_gens_out,
              'sample_size=i' => \$sample_size,
              'character_set=s' => \$fasta_character_set,
              'output_fasta_filename=s' => \$fasta_filename,
              'keep_prob=f' => \$keep_prob,
             );
 if(!defined $n_gens_out){
     $n_gens_out = $n_generations;
   }
   if (!defined $sample_size) {
      $sample_size = $n_gens_out * $pop;
   }
   if (!defined $keep_prob) {
      $keep_prob = $sample_size / ($n_gens_out * $pop);
    }
  
   print STDERR "sample size: $sample_size  keep prob: $keep_prob \n";


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

   if (!defined $fasta_filename) {
      my $mfspc = $mafspec;
      $mfspc =~ s/,/-/g;
      my $ssbit = "_samplesize" . ((defined $sample_size)? $sample_size : $n_gens_out * $pop);
      $fasta_filename =
        "pop" . $pop . "_nsnps" . $n_snps .
          "_ngens" . $n_generations . "-" . $n_gens_out .
            "_maf-" . $mfspc . $ssbit;
   }
   print STDERR "output file name: ", $fasta_filename, "\n";


   # ##################################################################################################
   # ########### get the initial generation of SimulatedGenotype objects with the given set of mafs ############

   my $mafs = GenotypeSimulation::draw_mafs($the_rng, $mafspec, $n_snps);

   # generate a data set of genotypes, by choosing an initial population,
   # and generating several generations ...
   my $this_generation = []; # a set of SimulatedGenotype objects from 1 generation.

   my ($generation, $id) = (0, 0);
   open my $fh_out, ">", $fasta_filename . '.fasta';
   for my $k (1..$pop) {      # generate the initial population of genotypes
      my $gobj = SimulatedGenotype->new_from_mafs($the_rng, $mafs, $generation, 'A' . $id);
      print STDERR "gen 0.  i: $k \n" if($k % 100 == 0);
      #   print STDERR "$keep_prob  $n_gens_out  $n_generations ", gsl_rng_uniform($the_rng->raw() ), "\n";
      if ($n_gens_out == $n_generations) { # if should output the initial generation
         if (
             ($keep_prob >= 1)
             or
             (gsl_rng_uniform($the_rng->raw() ) < $keep_prob)
            ) {
            #    print "ABCDEF \n";
            print $fh_out genotypes_as_fasta([$gobj], $fasta_character_set);
         
         }
         #  print "don't output this one $id \n";
      }
      #   print "X gen, id:  ", $gobj->get_generation(), "  ", $gobj->get_id(), "\n";
      push @$this_generation, $gobj;
      $id++;
   }
   # exit;
   # ###################################################################################################


   # ###################################################################################################
   # ########### get specified number of generations, each one produced by random mating of ######
   # ########### individuals in the preceding generation.  #########
   # in each gen, choose $pop pairs u.a.r., one offspring from each

   #  my @this_generation = @initial_generation;
   for my $i_gen (1..$n_generations-1) {
      my $next_generation = [];
      for my $k (1..$pop) {

         my ($i, $j) = ( gsl_rng_uniform_int($the_rng->raw(), scalar @$this_generation),
                         gsl_rng_uniform_int($the_rng->raw(), scalar @$this_generation)
                       );

         my $gobj = SimulatedGenotype->new_offspring($the_rng, $this_generation->[$i], $this_generation->[$j], $i_gen, 'A' . $id ); # 1st offspring
         print STDERR "gen $i_gen.  k: $k \n" if($k % 100 == 0);
         if ($i_gen >= $n_generations - $n_gens_out) { # if this is one of the generations we want to keep.
            if (
                ($keep_prob >= 1)
                or
                (gsl_rng_uniform($the_rng->raw() ) < $keep_prob)
               ) {
               print $fh_out genotypes_as_fasta([$gobj], $fasta_character_set);
            }
         }
         push @$next_generation, $gobj;
         $id++;
      }                # end loop over individuals in this generation.

      $this_generation = $next_generation;
   }                            # end loop over generations
   close $fh_out;
   # ############### done generating genotype generations ################################################


   #    my @kept_generations = @stored_generations[-$n_gens_out..-1]; # discard all but last $n_gens_out generations
   #    my @dataset_genotype_objects = ();
   #    for my $gobjs (@kept_generations) {
   #       push @dataset_genotype_objects, @$gobjs;
   #    }
   #    @dataset_genotype_objects = $the_rng->shuffle(@dataset_genotype_objects);
   #    # my @dataset_genotype_objects;
   #    if (defined $sample_size  and $sample_size < scalar @dataset_genotype_objects) {
   #       print STDERR "Outputting  ", $sample_size, " samples out of  ", scalar @dataset_genotype_objects, "\n";
   #       @dataset_genotype_objects = @dataset_genotype_objects[0..$sample_size-1];
   #    } else {
   #       @dataset_genotype_objects = @dataset_genotype_objects;
   #       print STDERR "Outputting  ", scalar @dataset_genotype_objects, " samples. \n";
   #    }

   #    @dataset_genotype_objects = sort { $a->get_id() <=> $b->get_id() } @dataset_genotype_objects;

   #    # ###################################################################################################
   #    # ###############  output  ######
   #    # fasta
   #    open my $fh_out, ">", $fasta_filename . '.fasta';
   #    #  for my $a_gen (@dataset_genotype_objects) {
   #    print $fh_out genotypes_as_fasta(\@dataset_genotype_objects, $fasta_character_set);
   #    #   }
   #    close $fh_out;

}


#############################################################################################################

sub genotypes_as_fasta{
   my $genotype_objects = shift; # array ref
   my $character_set = shift // '012'; 
   my $separator = shift // '';
   my $fasta_string = '';

   for my $gobj (@$genotype_objects) {
      my ($gen, $id, $pedigree) = ($gobj->get_generation(), 
                                   $gobj->get_id(),
                                   #  $gobj->get_parents(), 
                                   $gobj->get_pedigree()
                                  );
      my $gstring =  ($character_set eq '012')?
#	$gobj->genotype_012string($separator) . "\n" .
	$gobj->get_gtstring() :
          $gobj->genotype_Aa_string($separator);
      $fasta_string .= ">$id   $gen   $pedigree \n";
      $fasta_string .= "$gstring \n";
   }
   return $fasta_string;
}

