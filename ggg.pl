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

{                               # main
   my $the_rng = Math::GSL::RNG->new();

   # defaults:
   my $n_snps = 10;             # number of snps in each sample.
   #   my $relationship = 'PO'; # PO: parent-offspring. Or FS (full sib), HS (half sib), etc. ...
   #   my $n_reps = 100; # number or replicates (i.e. simulated sample pairs, each with $n_snps in both samples
   my $n_generations = 100;
   my $pop = 20;                # population size
   #   my $turnover = 5;     # this many births and deaths each generation
   my $mafspec = 'delta,0.25';  # minor allele frequency
   my $type = 3;
   my $n_gens_out = 4;
   my $sample_size = undef;

   GetOptions(
              'nsnps=i' => \$n_snps, # 
              'ngens|ngenerations=i' => \$n_generations,
              'population=i' => \$pop,
              #           'turnover=i' => \$turnover, 
              'maf|minor_allele_frequence=s' => \$mafspec, # 0 < $maf <= 0.5
              'type=i' => \$type,
              'nout=i' => \$n_gens_out,
              'sample_size=i' => \$sample_size,
             );

   my $mfspc = $mafspec;
   $mfspc =~ s/,/-/g;
   my $ssbit = "_samplesize" . ((defined $sample_size)? $sample_size : $n_gens_out * $pop);
   my $output_filename =
     "pop" . $pop . "_nsnps" . $n_snps .
       "_ngens" . $n_generations . "-" . $n_gens_out .
        "_maf-" . $mfspc . $ssbit;
   print STDERR "output file name: ", $output_filename, "\n";
# exit;

   my $mafs = GenotypeSim::draw_mafs($the_rng, $mafspec, $n_snps);

   # generate a data set of genotypes, by choosing an initial population,
   # and generation several generations ... 
   my @this_generation = (); # a set of Genotype objects from 1 generation.
   my @stored_generations = ();
   my %id_genotypeobj = ();
   my ($generation, $id) = (0, undef);
   for $id (1..$pop) {  # generate the initial population of genotypes
      my $gobj = Genotype->new_from_population($the_rng, $mafs, $generation, $id);
      push @this_generation, $gobj;
      $id_genotypeobj{$id} = $gobj;
   }
   push @stored_generations, \@this_generation; # if($n_gens_out >= $n_generations);

   # in each gen, choose $pop pairs u.a.r., one offspring from each
   for my $i_gen (1..$n_generations-1) {
      my @next_generation = ();
      for my $i_mating (1..$pop) {
         $id++;
         my ($i, $j) = ( gsl_rng_uniform_int($the_rng->raw(), scalar @this_generation),
                         gsl_rng_uniform_int($the_rng->raw(), scalar @this_generation)
                       );

         # print $i, " ", ref $this_generation[$i], "\n";
     #    my ($an, $ad, $hn, $hd, $gpc_count) = GenotypeSim::agmr_hgmr($this_generation[$i], $this_generation[$j]);
         #  print "$i$j   $an  $ad   $hn  $hd  ", $an/$ad, "  ", ($hd > 0)? $hn/$hd : '---',  "  ", GenotypeSim::paircode_count_string($gpc_count), "\n";

         my $gobj = Genotype->new_offspring($the_rng, $this_generation[$i], $this_generation[$j], $i_gen, $id ); # 1st offspring
         $id_genotypeobj{$id} = $gobj;
         push @next_generation, $gobj;
      }
      push @stored_generations, \@next_generation; # if($i_gen >= $n_generations - $n_gens_out);
      @this_generation = @next_generation;
   }

   for my $a_gen (@stored_generations){
      print genotypes_as_fasta($a_gen);
      # for my $gobj (@$a_gen){
      #    print $gobj->
      # }
   }

   my @kept_generations = @stored_generations[-$n_gens_out..-1]; # discard all but last $n_gens_out generations
   my @dataset_genotype_objects = ();
   for my $gobjs (@kept_generations){
      push @dataset_genotype_objects, @$gobjs;
   }
    @dataset_genotype_objects = $the_rng->shuffle(@dataset_genotype_objects);
   my @sample_genotype_objects;
   if (defined $sample_size  and $sample_size < scalar @dataset_genotype_objects) {
      print STDERR "Showing results for ", $sample_size, " samples out of  ", scalar @dataset_genotype_objects, "\n";
      @sample_genotype_objects = @dataset_genotype_objects[0..$sample_size-1];
   } else {
      @sample_genotype_objects = @dataset_genotype_objects;
      print STDERR "Showing results for all  ", scalar @dataset_genotype_objects, " samples. \n";
   }
   open my $fhout, ">", $output_filename;
   while (my ($i, $g1) = each @sample_genotype_objects) {
      print STDERR "$i out of ", scalar @sample_genotype_objects, "\n";
      for (my $j = $i+1; $j < scalar @sample_genotype_objects; $j++) {
         #   print STDERR "$i $j \n";
         my $g2 = $sample_genotype_objects[$j];
         #    GenotypeSim::print_genotypes_in_columns($g1, $g2);
         my ($an, $ad, $hn, $hd, $gpc_count) = GenotypeSim::agmr_hgmr($g1, $g2);
         printf $fhout "%3i %3i %4i %4i %4i %4i  %5.4f %5.4f   ", ($i, $j, $an, $ad, $hn, $hd, $an/$ad, ($hd > 0)? $hn/$hd : -100);
         printf $fhout "%s \n", GenotypeSim::paircode_count_string($gpc_count);
      }
   }
}

sub genotypes_as_fasta{
   my $set_of_genotypes = shift;
   my $fasta_string = '';

   for my $gobj (@$set_of_genotypes) {
      my ($id, $gtstring, $gen, $parents) = ( $gobj->get_id(), $gobj->genotype_string(), $gobj->get_generation(), $gobj->get_parents() );
      $fasta_string .= ">$id   $gen $parents \n";
      $fasta_string .= "$gtstring \n";
   }
   return $fasta_string;
}
