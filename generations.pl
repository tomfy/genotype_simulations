#!/usr/bin/perl
use strict;
use warnings;
use List::Util qw(min max sum);
# use Graphics::GnuplotIF qw(GnuplotIF);
use Getopt::Long;
use Math::GSL::RNG  qw( :all );
use GenotypeSim qw ( :all );
use Genotype;

# the idea here is to start with a set of unrelated genotypes
# and let them pair up and produce offspring,
# 

{
   my $the_rng = Math::GSL::RNG->new();

   # defaults:
   my $n_snps = 10;             # number of snps in each sample.
   #   my $relationship = 'PO'; # PO: parent-offspring. Or FS (full sib), HS (half sib), etc. ...
   #   my $n_reps = 100; # number or replicates (i.e. simulated sample pairs, each with $n_snps in both samples
   my $n_generations = 100;
   my $pop = 20;         # population size
   my $turnover = 5;     # this many births and deaths each generation
   my $mafspec = 'delta,0.25';  # minor allele frequency
   my $type = 3;
   my $n_gens_out = 4;
my $sample_size = undef;
 
   GetOptions(
              'nsnps=i' => \$n_snps, # 
              'ngenerations=i' => \$n_generations,
              'population=i' => \$pop,
              'turnover=i' => \$turnover, 
              'maf|minor_allele_frequence=s' => \$mafspec, # 0 < $maf <= 0.5
             'type=i' => \$type,
              'nout=i' => \$n_gens_out,
              'sample_size=i' => \$sample_size,
             );

   my $mafs = GenotypeSim::draw_mafs($the_rng, $mafspec, $n_snps);

   # generate a data set of genotypes, by choosing an initial population,
   # and generation several generations ... 
   my @genotypes = ();
   my @stored_genotypes = ();
   my %id_genotypeobj = ();
   my $id = undef;
   for $id = (1..$pop) {      # generate the initial population of genotypes
        my $gobj = Genotype->new_from_mafs($the_rng, $mafs, $_);
         push @genotypes, $gobj;
      $id_genotypeobj{$id} = 
   }

   if ($type == 1) { # one offspring, replacing one of old genotypes, per generation
      for (1..$n_generations) {
         my ($i, $j) = ( gsl_rng_uniform_int($the_rng->raw(), $pop),  gsl_rng_uniform_int($the_rng->raw(), $pop) );
         my $offspring = GenotypeSim::draw_offspring_genotype($the_rng, $genotypes[$i], $genotypes[$j] );
         my $i_die =  gsl_rng_uniform_int($the_rng->raw(), $n_snps); # this one 'dies', i.e. gets replaced by the new offspring genotype.
         $genotypes[$i_die] = $offspring;
         my $i_data = gsl_rng_uniform_int($the_rng->raw(), $pop); # this one gets 'observed' i.e. stored in the data set which will then be analyzed
         push @stored_genotypes, $genotypes[$i_data];
      }

   } elsif($type == 2) {
      for (1..$n_generations) {
         my ($i, $j, $k) = ( gsl_rng_uniform_int($the_rng->raw(), scalar @genotypes),
                             gsl_rng_uniform_int($the_rng->raw(), scalar @genotypes),
                             gsl_rng_uniform_int($the_rng->raw(), scalar @genotypes),
                           );
      
         push @genotypes, GenotypeSim::draw_offspring_genotype($the_rng, $genotypes[$i], $genotypes[$j] ); # 1st offspring
         if (1 or gsl_rng_uniform($the_rng->raw()) < 0.75) {
            push @genotypes, GenotypeSim::draw_offspring_genotype($the_rng, $genotypes[$i], $genotypes[$j] ); # if do this, it is FS of 1st offspring
         }
         push @genotypes, GenotypeSim::draw_offspring_genotype($the_rng, $genotypes[$i], $genotypes[$k] );
         #   my $i_die =  gsl_rng_uniform_int($the_rng->raw(), $n_snps); # this one 'dies', i.e. gets replaced by the new offspring genotype.
         #   $genotypes[$i_die] = $offspring;
         #   my $i_data = gsl_rng_uniform_int($the_rng->raw(), $pop); # this one gets 'observed' i.e. stored in the data set which will then be analyzed
         #      push @dataset_genotypes, $genotypes[$i_data];
      }
      @stored_genotypes = @genotypes;
   } elsif ($type == 3) { # in each gen, choose $pop pairs u.a.r., one offspring from each
      for my $i_gen (1..$n_generations) {
         my @new_genotypes = ();
         for my $i_mating (1..$pop) {
            my ($i, $j) = ( gsl_rng_uniform_int($the_rng->raw(), scalar @genotypes),
                            gsl_rng_uniform_int($the_rng->raw(), scalar @genotypes)
                          );

            push @new_genotypes, GenotypeSim::draw_offspring_genotype($the_rng, $genotypes[$i], $genotypes[$j] ); # 1st offspring
         }
         push @stored_genotypes, @new_genotypes if($i_gen > $n_generations - $n_gens_out);
         @genotypes = @new_genotypes;
      }
   }

   print STDERR "n samples: ", scalar @stored_genotypes, "\n";

   my @dataset_genotypes = $the_rng->shuffle(@stored_genotypes);
   
my @sample_genotypes =  (defined $sample_size)? @dataset_genotypes[0..$sample_size-1] : @dataset_genotypes;
   while (my ($i, $g1) = each @sample_genotypes) {
      for (my $j = $i+1; $j < scalar @sample_genotypes; $j++) {
         #   print STDERR "$i $j \n";
         my $g2 = $sample_genotypes[$j];
         #    GenotypeSim::print_genotypes_in_columns($g1, $g2);
         my ($an, $ad, $hn, $hd, $gpc_count) = GenotypeSim::agmr_hgmr($g1, $g2);
         print "$i$j   $an  $ad   $hn  $hd  ", $an/$ad, "  ", ($hd > 0)? $hn/$hd : '---',  "  ", GenotypeSim::paircode_count_string($gpc_count), "\n";
      }
   }
}
