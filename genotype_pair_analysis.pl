#!/usr/bin/perl -w
use strict;
use warnings;
use List::Util qw(min max sum);
# use Graphics::GnuplotIF qw(GnuplotIF);
use Getopt::Long;
use Math::GSL::RNG  qw( :all );
use GenotypeSimulation qw ( :all );
use SimulatedGenotype;

# genotype_pair_analysis.pl : read genotypes in fasta format, 
# analyze them pairwise (agmr, hgmr, etc.)

{                               # main

   # defaults:
   my $sample_size = undef;
   my $rng_type = $gsl_rng_mt19937;
   my $rng_seed = undef;
   my $known_pairs_only = 0;

   GetOptions(
              'sample_size|size=i' => \$sample_size, # 
              'rng_type=s' => \$rng_type,
              'seed=s' => \$rng_seed,
              'known_pairs_only!' => \$known_pairs_only,
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

   # read in genotypes from fasta file,
   my @genotype_objects = ();
   while (<>) {
      if ( /^> (\S+) \s+ (\S+) \s+ (\S+) /x ) {
         my ($id, $gen, $pedigree) = ($1, $2, $3);
      #   print STDERR "$id $gen $pedigree \n";
         my $string = <>;
         my $gobj = SimulatedGenotype->new_from_012string($the_rng, $string, $gen, \$id, $pedigree);
     #    print STDERR $gobj->get_pedigree(), "\n";
         push @genotype_objects, $gobj;
      }
   }

   # ###########   get subset of size $sample_size pairs  ################################
   if (defined $sample_size  and  (2*$sample_size < scalar @genotype_objects)) {
   #   @genotype_objects = $the_rng->shuffle(@genotype_objects);
      @genotype_objects = @genotype_objects[0..2*$sample_size-1];
   }

   # ########################    compare pairs   ####################################

#   my $output_filename;
#   open my $fhout, ">", $output_filename;

   if($known_pairs_only){ # just compare 0 with 1, 2 with 3, etc.
      for (my $i=0; $i < scalar @genotype_objects; $i += 2) {
         my ($g1, $g2) = @genotype_objects[$i,$i+1];
         my ($an, $ad, $hn, $hd, $gpc_count) = GenotypeSimulation::agmr_hgmr($g1, $g2);
         printf  "%3i %3i %4i %4i %4i %4i  %5.4f %5.4f   ", 
           ($g1->get_id(), $g2->get_id(), $an, $ad, $hn, $hd, $an/$ad, ($hd > 0)? $hn/$hd : -100);
 
         
         printf  "%s   ", GenotypeSimulation::paircode_count_string($gpc_count);
         printf  "%s  %s \n", $g1->get_pedigree(), $g2->get_pedigree;

      }
   } else {                     # all vs all
      while (my ($i, $g1) = each @genotype_objects) {
         print STDERR "$i out of ", scalar @genotype_objects, "\n";
     
         for (my $j = $i; $j < scalar @genotype_objects; $j++) {
            #   print STDERR "$i $j \n";
            my $g2 = $genotype_objects[$j];
            #    GenotypeSimulation::print_genotypes_in_columns($g1, $g2);
            my ($an, $ad, $hn, $hd, $gpc_count) = GenotypeSimulation::agmr_hgmr($g1, $g2);
            printf  "%3i %3i %4i %4i %4i %4i  %5.4f %5.4f   ", 
              ($g1->get_id(), $g2->get_id(), $an, $ad, $hn, $hd, $an/$ad, ($hd > 0)? $hn/$hd : -100);
 
         
            printf  "%s   ", GenotypeSimulation::paircode_count_string($gpc_count);
            printf  "%s  %s \n", $g1->get_pedigree(), $g2->get_pedigree;
         }
      }
   }

} ############################     end of main   ####################################
