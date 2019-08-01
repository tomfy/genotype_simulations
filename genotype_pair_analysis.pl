#!/usr/bin/perl -w
use strict;
use warnings;
use List::Util qw(min max sum);
# use Graphics::GnuplotIF qw(GnuplotIF);
use Getopt::Long;
use Math::GSL::RNG  qw( :all );
use GenotypeSim qw ( :all );
use Genotype;

# ggg.pl : "generate generations of genotypes"

# the idea here is to start with a set of unrelated genotypes
# and let them pair up randomly and produce offspring,
# for some specified number of generations.

{                               # main


   # defaults:
   my $sample_size = undef;

 GetOptions(
              'sample_size|size=i' => \$sample_size, # 
             );

# read in genotypes from fasta file,


# choose a sample of sample_size u.a.r. if sample_size defined.

# compare pairs

my $output_filename;
   open my $fhout, ">", $output_filename;
   while (my ($i, $g1) = each @sample_genotype_objects) {
      print STDERR "$i out of ", scalar @sample_genotype_objects, "\n";
      for (my $j = $i; $j < scalar @sample_genotype_objects; $j++) {
         #   print STDERR "$i $j \n";
         my $g2 = $sample_genotype_objects[$j];
         #    GenotypeSim::print_genotypes_in_columns($g1, $g2);
         my ($an, $ad, $hn, $hd, $gpc_count) = GenotypeSim::agmr_hgmr($g1, $g2);
         printf $fhout "%3i %3i %4i %4i %4i %4i  %5.4f %5.4f   ", ($g1->get_id(), $g2->get_id(), $an, $ad, $hn, $hd, $an/$ad, ($hd > 0)? $hn/$hd : -100);
         printf $fhout "%s \n", GenotypeSim::paircode_count_string($gpc_count);
      }
   }
}
