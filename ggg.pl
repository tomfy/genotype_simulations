#!/usr/bin/perl
use strict;
use warnings;
use List::Util qw(min max sum);
# use Graphics::GnuplotIF qw(GnuplotIF);
use Getopt::Long;
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
   my $rng_seed = undef;
   my $n_gens_out = $n_generations;
   my $sample_size = undef;
   my $fasta_character_set = '012'; # anything else gives aA, aa, AA.
   my $fasta_filename = undef;

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

   if(!defined $fasta_filename){
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

   # my $g0 = SimulatedGenotype->new_from_mafs($the_rng, $mafs, 0, 666);
   # my $str012 = $g0->genotype_012string();
   # print "$str012 \n", $g0->genotype_string(), "\n";
   # my $g1 = SimulatedGenotype->new_from_012string($the_rng, $str012, 0, 777);
   # my $str012_1 = $g1->genotype_012string();
   # print "$str012_1 \n", $g1->genotype_string(), "\n";
   # exit;


   # generate a data set of genotypes, by choosing an initial population,
   # and generation several generations ... 
   my @initial_generation = (); # a set of SimulatedGenotype objects from 1 generation.
   my @stored_generations = ();
   my %id_genotypeobj = ();
   my ($generation, $id) = (0, 0);
   for (1..$pop) {      # generate the initial population of genotypes
      my $gobj = SimulatedGenotype->new_from_mafs($the_rng, $mafs, $generation, \$id);
      #   print "X gen, id:  ", $gobj->get_generation(), "  ", $gobj->get_id(), "\n";
      push @initial_generation, $gobj;
      $id_genotypeobj{$id} = $gobj;
      $id++;
   }
 
   push @stored_generations, \@initial_generation; # if($n_gens_out >= $n_generations);

   # ###################################################################################################


   # ###################################################################################################
   # ########### get specified number of generations, each one produced by random mating of ######
   # ########### individuals in the preceding generation.  #########
   # in each gen, choose $pop pairs u.a.r., one offspring from each

   my @this_generation = @initial_generation;
   for my $i_gen (1..$n_generations-1) {
      my @next_generation = ();
      for (1..$pop) {

         my ($i, $j) = ( gsl_rng_uniform_int($the_rng->raw(), scalar @this_generation),
                         gsl_rng_uniform_int($the_rng->raw(), scalar @this_generation)
                       );

         my $gobj = SimulatedGenotype->new_offspring($the_rng, $this_generation[$i], $this_generation[$j], $i_gen, \$id ); # 1st offspring
         $id_genotypeobj{$id} = $gobj;
         push @next_generation, $gobj;
         $id++;
      }
      push @stored_generations, \@next_generation; 
      @this_generation = @next_generation;
   }

   # ############### done generating genotype generations ################################################


   my @kept_generations = @stored_generations[-$n_gens_out..-1]; # discard all but last $n_gens_out generations
   my @dataset_genotype_objects = ();
   for my $gobjs (@kept_generations) {
      push @dataset_genotype_objects, @$gobjs;
   }
   @dataset_genotype_objects = $the_rng->shuffle(@dataset_genotype_objects);
  # my @dataset_genotype_objects;
   if (defined $sample_size  and $sample_size < scalar @dataset_genotype_objects) {
      print STDERR "Outputting  ", $sample_size, " samples out of  ", scalar @dataset_genotype_objects, "\n";
      @dataset_genotype_objects = @dataset_genotype_objects[0..$sample_size-1];
   } else {
      @dataset_genotype_objects = @dataset_genotype_objects;
      print STDERR "Outputting  ", scalar @dataset_genotype_objects, " samples. \n";
   }

   @dataset_genotype_objects = sort { $a->get_id() <=> $b->get_id() } @dataset_genotype_objects;

  # ###################################################################################################
  # ###############  output  ######
  # fasta
   open my $fh_out, ">", $fasta_filename . '.fasta';
 #  for my $a_gen (@dataset_genotype_objects) {
      print $fh_out genotypes_as_fasta(\@dataset_genotype_objects, $fasta_character_set);
#   }
   close $fh_out;

}

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
        $gobj->genotype_012string($separator) :
          $gobj->genotype_aAstring($separator);
      $fasta_string .= ">$id   $gen   $pedigree \n";
      $fasta_string .= "$gstring \n";
   }
   return $fasta_string;
}

