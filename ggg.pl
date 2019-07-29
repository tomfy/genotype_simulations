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

{ # main
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
              'ngens|ngenerations=i' => \$n_generations,
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
   for $id (1..$pop) { # generate the initial population of genotypes
      my $gobj = Genotype->new_from_population($the_rng, $mafs, $id);
      push @genotypes, $gobj;
      $id_genotypeobj{$id} = $gobj;
   }
push @stored_genotypes, @genotypes if($n_gens_out >= $n_generations);

   # in each gen, choose $pop pairs u.a.r., one offspring from each
   for my $i_gen (1..$n_generations-1) {
      my @new_genotypes = ();
      for my $i_mating (1..$pop) {
         $id++;
         my ($i, $j) = ( gsl_rng_uniform_int($the_rng->raw(), scalar @genotypes),
                         gsl_rng_uniform_int($the_rng->raw(), scalar @genotypes)
                       );

# print $i, " ", ref $genotypes[$i], "\n";
         my ($an, $ad, $hn, $hd, $gpc_count) = GenotypeSim::agmr_hgmr($genotypes[$i], $genotypes[$j]);
       #  print "$i$j   $an  $ad   $hn  $hd  ", $an/$ad, "  ", ($hd > 0)? $hn/$hd : '---',  "  ", GenotypeSim::paircode_count_string($gpc_count), "\n";

         my $gobj = Genotype->new_offspring($the_rng, $genotypes[$i], $genotypes[$j], $id ); # 1st offspring
         $id_genotypeobj{$id} = $gobj;
         push @new_genotypes, $gobj;
      }
      push @stored_genotypes, @new_genotypes if($i_gen >= $n_generations - $n_gens_out);
      @genotypes = @new_genotypes;
   }

   @stored_genotypes = $the_rng->shuffle(@stored_genotypes);

   my @sample_genotypes;
   if (defined $sample_size  and $sample_size < scalar @stored_genotypes) {
     
      print STDERR "Showing results for ", $sample_size, " samples out of  ", scalar @stored_genotypes, "\n";
      @sample_genotypes = @stored_genotypes[0..$sample_size-1];
   } else {
      @sample_genotypes = @stored_genotypes;
      print STDERR "Showing results for all  ", scalar @stored_genotypes, " samples. \n";
   }
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
