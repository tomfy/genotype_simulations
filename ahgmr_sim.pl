#!/usr/bin/perl
use strict;
use warnings;
use List::Util qw(min max sum);
# use Graphics::GnuplotIF qw(GnuplotIF);
use Getopt::Long;
use GenotypeSim qw ( :all );


{
   my $the_rng = Math::GSL::RNG->new();

   # defaults:
   my $n_snps = 10;             # number of snps in each sample.
   my $relationship = 'PO'; # PO: parent-offspring. Or FS (full sib), HS (half sib), etc. ...
   my $n_reps = 100; # number or replicates (i.e. simulated sample pairs, each with $n_snps in both samples
   my $mafspec = 'delta,0.25';      # minor allele frequency
 
   GetOptions(
              'nsnps=i' => \$n_snps, # 
              'relationship=s' => \$relationship,
              'nreps=i' => \$n_reps,
              'maf|minor_allele_frequence=s' => \$mafspec, # 0 < $maf <= 0.5
             );

   my $mafs = GenotypeSim::draw_mafs($the_rng, $mafspec, $n_snps);
 #  print STDERR "AAA: ", join(", ", @$mafs), "\n";

   for (1..$n_reps) {

      if ($relationship eq 'DPO') { # both parents are same individual
         my $p1_genotype = GenotypeSim::draw_genotype_from_mafs($the_rng, $mafs);
         my $o_genotype = GenotypeSim::draw_offspring_genotype($the_rng, $p1_genotype, $p1_genotype);
         my ($an, $ad, $hn, $hd, $gpc_count) = GenotypeSim::agmr_hgmr($p1_genotype, $o_genotype);
         print "DPO   $an  $ad   $hn  $hd  ", $an/$ad, "  ", ($hd > 0)? $hn/$hd : '---',  "  ", GenotypeSim::paircode_count_string($gpc_count), "\n";
      } elsif ($relationship eq 'PO') { # parent and offspring
         #print STDERR "ZZZZZZZZZZ \n";
         my $p1_genotype = GenotypeSim::draw_genotype_from_mafs($the_rng, $mafs);
         #print STDERR "XXXXXXXXXXXXXX \n";
         my $p2_genotype = GenotypeSim::draw_genotype_from_mafs($the_rng, $mafs);
         my $o_genotype = GenotypeSim::draw_offspring_genotype($the_rng, $p1_genotype, $p2_genotype);

      #   print GenotypeSim::genotype_string($p1_genotype), "\n";
        #    GenotypeSim::print_genotypes_in_columns($p1_genotype, $p2_genotype, $o_genotype);
         my ($an, $ad, $hn, $hd, $gpc_count) = GenotypeSim::agmr_hgmr($p1_genotype, $o_genotype);
         print "PO   $an  $ad   $hn  $hd  ", $an/$ad, "  ", ($hd > 0)? $hn/$hd : '---', "  ", GenotypeSim::paircode_count_string($gpc_count), "\n";
      } elsif ($relationship eq 'FS') { # full siblings
         my $p1_genotype = GenotypeSim::draw_genotype_from_mafs($the_rng, $mafs);
         my $p2_genotype = GenotypeSim::draw_genotype_from_mafs($the_rng, $mafs);
         my $o1_genotype = GenotypeSim::draw_offspring_genotype($the_rng, $p1_genotype, $p2_genotype);
         my $o2_genotype = GenotypeSim::draw_offspring_genotype($the_rng, $p1_genotype, $p2_genotype);
         my ($an, $ad, $hn, $hd, $gpc_count) = GenotypeSim::agmr_hgmr($o1_genotype, $o2_genotype);
         print "FS   $an  $ad   $hn  $hd  ", $an/$ad, "  ", ($hd > 0)? $hn/$hd : '---',  "  ", GenotypeSim::paircode_count_string($gpc_count), "\n";
      }elsif($relationship eq 'DFS') { # 'double full siblings ' i.e. 
         my $p1_genotype = GenotypeSim::draw_genotype_from_mafs($the_rng, $mafs);
         my $o1_genotype = GenotypeSim::draw_offspring_genotype($the_rng, $p1_genotype, $p1_genotype);
         my $o2_genotype = GenotypeSim::draw_offspring_genotype($the_rng, $p1_genotype, $p1_genotype);
         my ($an, $ad, $hn, $hd, $gpc_count) = GenotypeSim::agmr_hgmr($o1_genotype, $o2_genotype);
         print "DFS   $an  $ad   $hn  $hd  ", $an/$ad, "  ", ($hd > 0)? $hn/$hd : '---',  "  ", GenotypeSim::paircode_count_string($gpc_count), "\n";
      } elsif ($relationship eq 'FC') { # first cousins
         my $gp1_genotype = GenotypeSim::draw_genotype_from_mafs($the_rng, $mafs);
         my $gp2_genotype = GenotypeSim::draw_genotype_from_mafs($the_rng, $mafs);

         my $p1_genotype = GenotypeSim::draw_offspring_genotype($the_rng, $gp1_genotype, $gp2_genotype); # p1, p2 are full sibling
         my $p2_genotype = GenotypeSim::draw_offspring_genotype($the_rng, $gp1_genotype, $gp2_genotype);

         my $p3_genotype = GenotypeSim::draw_genotype_from_mafs($the_rng, $mafs); # p3, p4 are unrelated
         my $p4_genotype = GenotypeSim::draw_genotype_from_mafs($the_rng, $mafs);

         my $c1_genotype = GenotypeSim::draw_offspring_genotype($the_rng, $p1_genotype, $p3_genotype);
         my $c2_genotype = GenotypeSim::draw_offspring_genotype($the_rng, $p2_genotype, $p4_genotype);
         my ($an, $ad, $hn, $hd, $gpc_count) = GenotypeSim::agmr_hgmr($c1_genotype, $c2_genotype);
         print "FC   $an  $ad   $hn  $hd  ", $an/$ad, "  ", ($hd > 0)? $hn/$hd : '---',  "  ", GenotypeSim::paircode_count_string($gpc_count), "\n";
      } elsif ($relationship eq 'AUNN') { # aunt/uncle - niece/nephew
         my $gp1_genotype = GenotypeSim::draw_genotype_from_mafs($the_rng, $mafs);
         my $gp2_genotype = GenotypeSim::draw_genotype_from_mafs($the_rng, $mafs);

         my $p1_genotype = GenotypeSim::draw_offspring_genotype($the_rng, $gp1_genotype, $gp2_genotype); # p1, p2 are full siblings
         my $p2_genotype = GenotypeSim::draw_offspring_genotype($the_rng, $gp1_genotype, $gp2_genotype);

         my $p3_genotype = GenotypeSim::draw_genotype_from_mafs($the_rng, $mafs); # p3, p4 are unrelated
         #   my $p4_genotype = GenotypeSim::draw_genotype_from_mafs($the_rng, $mafs);

         my $nn_genotype = GenotypeSim::draw_offspring_genotype($the_rng, $p1_genotype, $p3_genotype); # offspring of p1 and p3, niece/nephew of p2
         #    my $c2_genotype = GenotypeSim::draw_offspring_genotype($the_rng, $p2_genotype, $p4_genotype);
         my ($an, $ad, $hn, $hd, $gpc_count) = GenotypeSim::agmr_hgmr($p2_genotype, $nn_genotype);
         print "AUNN   $an  $ad   $hn  $hd  ", $an/$ad, "  ", ($hd > 0)? $hn/$hd : '---',  "  ", GenotypeSim::paircode_count_string($gpc_count), "\n";
      } elsif ($relationship eq 'GPGC') { # grandparent - grandchild
         my $gp1_genotype = GenotypeSim::draw_genotype_from_mafs($the_rng, $mafs);
         my $gp2_genotype = GenotypeSim::draw_genotype_from_mafs($the_rng, $mafs);

         my $p1_genotype = GenotypeSim::draw_offspring_genotype($the_rng, $gp1_genotype, $gp2_genotype); # offspring of gp1, gp2
         my $p3_genotype = GenotypeSim::draw_genotype_from_mafs($the_rng, $mafs); # other parent - unrelated
         #   my $p4_genotype = GenotypeSim::draw_genotype_from_mafs($the_rng, $mafs);

         my $gc_genotype = GenotypeSim::draw_offspring_genotype($the_rng, $p1_genotype, $p3_genotype); # offspring of p1 and p3, gc of gp1 & gp2
         #    my $c2_genotype = GenotypeSim::draw_offspring_genotype($the_rng, $p2_genotype, $p4_genotype);
         my ($an, $ad, $hn, $hd, $gpc_count) = GenotypeSim::agmr_hgmr($gp2_genotype, $gc_genotype);
         print "GPGC   $an  $ad   $hn  $hd  ", $an/$ad, "  ", ($hd > 0)? $hn/$hd : '---',  "  ", GenotypeSim::paircode_count_string($gpc_count), "\n";
      } elsif ($relationship eq 'HS') { # half-siblings
         my $p1_genotype = GenotypeSim::draw_genotype_from_mafs($the_rng, $mafs);
         my $p2_genotype = GenotypeSim::draw_genotype_from_mafs($the_rng, $mafs);
         my $p3_genotype = GenotypeSim::draw_genotype_from_mafs($the_rng, $mafs);
         my $o12_genotype = GenotypeSim::draw_offspring_genotype($the_rng, $p1_genotype, $p2_genotype);
         my $o23_genotype = GenotypeSim::draw_offspring_genotype($the_rng, $p2_genotype, $p3_genotype);
         my ($an, $ad, $hn, $hd, $gpc_count) = GenotypeSim::agmr_hgmr($o12_genotype, $o23_genotype);
         print "HS   $an  $ad   $hn  $hd  ", $an/$ad, "  ", ($hd > 0)? $hn/$hd : '---',  "  ", GenotypeSim::paircode_count_string($gpc_count), "\n";
      } elsif ($relationship eq 'HAUNN') { # half-aunt/uncle - half-niece/nephew
         my $p1_genotype = GenotypeSim::draw_genotype_from_mafs($the_rng, $mafs);
         my $p2_genotype = GenotypeSim::draw_genotype_from_mafs($the_rng, $mafs);
         my $p3_genotype = GenotypeSim::draw_genotype_from_mafs($the_rng, $mafs);

         my $o12_genotype = GenotypeSim::draw_offspring_genotype($the_rng, $p1_genotype, $p2_genotype);
         my $o23_genotype = GenotypeSim::draw_offspring_genotype($the_rng, $p2_genotype, $p3_genotype);

         my $p4_genotype = GenotypeSim::draw_genotype_from_mafs($the_rng, $mafs);
         my $hnn_genotype = GenotypeSim::draw_offspring_genotype($the_rng, $o12_genotype, $p4_genotype);
         my ($an, $ad, $hn, $hd, $gpc_count) = GenotypeSim::agmr_hgmr($hnn_genotype, $o23_genotype);
         print "HAUNN   $an  $ad   $hn  $hd  ", $an/$ad, "  ", ($hd > 0)? $hn/$hd : '---',  "  ", GenotypeSim::paircode_count_string($gpc_count), "\n";
      } elsif ($relationship eq 'UN') { # unrelated
         my $genotype1 = GenotypeSim::draw_genotype_from_mafs($the_rng, $mafs);
         my $genotype2 = GenotypeSim::draw_genotype_from_mafs($the_rng, $mafs);
         my ($an, $ad, $hn, $hd, $gpc_count) = GenotypeSim::agmr_hgmr($genotype1, $genotype2);
         print "UN   $an  $ad   $hn  $hd  ", $an/$ad, "  ", ($hd > 0)? $hn/$hd : '---',  "  ", GenotypeSim::paircode_count_string($gpc_count), "\n";
      }
   }
}

