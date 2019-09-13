#!/usr/bin/perl
use strict;
use warnings;
use List::Util qw(min max sum);
# use Graphics::GnuplotIF qw(GnuplotIF);
use Getopt::Long;
use SimulatedGenotype qw ( :all );
use GenotypeSimulation qw ( :all );


{
   my $the_rng = Math::GSL::RNG->new();

   # defaults:
   my $n_snps = 10;             # number of snps in each sample.
   #  my $relationship = 'PO'; # PO: parent-offspring. Or FS (full sib), HS (half sib), etc. ...
   my $n_pairs = 100; # number of simulated sample pairs or each relationship specified by $rel_string;
   my $mafspec = 'delta,0.25';  # minor allele frequency
   my $rel_string = 'PO,FS,HS,FC'; # 
   my $fasta_character_set = '012';
 
   GetOptions(
              'nsnps=i' => \$n_snps, # 
              #  'relationship=s' => \$relationship,
              'npairs|nreps=i' => \$n_pairs,
              'maf|minor_allele_frequence=s' => \$mafspec, # 0 < $maf <= 0.5
              'rels|relationships=s' => \$rel_string, 
             );

   my $mafs = GenotypeSimulation::draw_mafs($the_rng, $mafspec, $n_snps);
   #  print STDERR "AAA: ", join(", ", @$mafs), "\n";
   my @rels = split(",", $rel_string);
   my @genotype_objects = ();
my ($gen, $id) = (0, 0);
   for my $relationship (@rels) {

      for (1..$n_pairs) { # simulate $n_pairs genotype pairs with $relationship

         if ($relationship eq 'DPO') { # both parents are same individual
            my $p1_genotype = Genotype->new_from_mafs($the_rng, $mafs, $gen, \$id);
            my $o_genotype = Genotype->new_offspring($the_rng, $p1_genotype, $p1_genotype, $gen, \$id);
            push @genotype_objects, ($p1_genotype, $o_genotype);
            #      my ($an, $ad, $hn, $hd, $gpc_count) = GenotypeSimulation::agmr_hgmr($p1_genotype, $o_genotype, $gen, \$id);
            #      print "DPO   $an  $ad   $hn  $hd  ", $an/$ad, "  ", ($hd > 0)? $hn/$hd : '---',  "  ", GenotypeSimulation::paircode_count_string($gpc_count), "\n";
         } elsif ($relationship eq 'PO') { # parent and offspring
            #print STDERR "ZZZZZZZZZZ \n";
            my $p1_genotype = SimulatedGenotype->new_from_mafs($the_rng, $mafs, $gen, \$id);
            #print STDERR "XXXXXXXXXXXXXX \n";
            my $p2_genotype = SimulatedGenotype->new_from_mafs($the_rng, $mafs, $gen, \$id);
            my $o_genotype = SimulatedGenotype->new_offspring($the_rng, $p1_genotype, $p2_genotype, $gen, \$id);
         #   print STDERR "refs of gobjs: ", ref $p1_genotype, "  ", ref $o_genotype, "\n";
            push @genotype_objects, ($p1_genotype, $o_genotype);
            #   print GenotypeSimulation::genotype_string($p1_genotype), "\n";
            #    GenotypeSimulation::print_genotypes_in_columns($p1_genotype, $p2_genotype, $o_genotype, $gen, \$id);
            #       my ($an, $ad, $hn, $hd, $gpc_count) = GenotypeSimulation::agmr_hgmr($p1_genotype, $o_genotype, $gen, \$id);
            #       print "PO   $an  $ad   $hn  $hd  ", $an/$ad, "  ", ($hd > 0)? $hn/$hd : '---', "  ", GenotypeSimulation::paircode_count_string($gpc_count), "\n";
         } elsif ($relationship eq 'FS') { # full siblings
            my $p1_genotype = SimulatedGenotype->new_from_mafs($the_rng, $mafs, $gen, \$id);
            my $p2_genotype = SimulatedGenotype->new_from_mafs($the_rng, $mafs, $gen, \$id);
            my $o1_genotype = SimulatedGenotype->new_offspring($the_rng, $p1_genotype, $p2_genotype, $gen, \$id);
            my $o2_genotype = SimulatedGenotype->new_offspring($the_rng, $p1_genotype, $p2_genotype, $gen, \$id);
            push @genotype_objects, ($o1_genotype, $o2_genotype);
            #        my ($an, $ad, $hn, $hd, $gpc_count) = GenotypeSimulation::agmr_hgmr($o1_genotype, $o2_genotype, $gen, \$id);
            #        print "FS   $an  $ad   $hn  $hd  ", $an/$ad, "  ", ($hd > 0)? $hn/$hd : '---',  "  ", GenotypeSimulation::paircode_count_string($gpc_count), "\n";
         } elsif ($relationship eq 'DFS') { # 'double full siblings ' i.e. 
            my $p1_genotype = SimulatedGenotype->new_from_mafs($the_rng, $mafs, $gen, \$id);
            my $o1_genotype = SimulatedGenotype->new_offspring($the_rng, $p1_genotype, $p1_genotype, $gen, \$id);
            my $o2_genotype = SimulatedGenotype->new_offspring($the_rng, $p1_genotype, $p1_genotype, $gen, \$id);
            push @genotype_objects, ($o1_genotype, $o2_genotype);
            #       my ($an, $ad, $hn, $hd, $gpc_count) = GenotypeSimulation::agmr_hgmr($o1_genotype, $o2_genotype, $gen, \$id);
            #       print "DFS   $an  $ad   $hn  $hd  ", $an/$ad, "  ", ($hd > 0)? $hn/$hd : '---',  "  ", GenotypeSimulation::paircode_count_string($gpc_count), "\n";
         } elsif ($relationship eq 'FC') { # first cousins
            my $gp1_genotype = SimulatedGenotype->new_from_mafs($the_rng, $mafs, $gen, \$id);
            my $gp2_genotype = SimulatedGenotype->new_from_mafs($the_rng, $mafs, $gen, \$id);

            my $p1_genotype = SimulatedGenotype->new_offspring($the_rng, $gp1_genotype, $gp2_genotype, $gen, \$id); # p1, p2 are full sibling
            my $p2_genotype = SimulatedGenotype->new_offspring($the_rng, $gp1_genotype, $gp2_genotype, $gen, \$id);

            my $p3_genotype = SimulatedGenotype->new_from_mafs($the_rng, $mafs, $gen, \$id); # p3, p4 are unrelated
            my $p4_genotype = SimulatedGenotype->new_from_mafs($the_rng, $mafs, $gen, \$id);

            my $c1_genotype = SimulatedGenotype->new_offspring($the_rng, $p1_genotype, $p3_genotype, $gen, \$id);
            my $c2_genotype = SimulatedGenotype->new_offspring($the_rng, $p2_genotype, $p4_genotype, $gen, \$id);
            push @genotype_objects, ($c1_genotype, $c2_genotype);
            #     my ($an, $ad, $hn, $hd, $gpc_count) = GenotypeSimulation::agmr_hgmr($c1_genotype, $c2_genotype, $gen, \$id);
            #      print "FC   $an  $ad   $hn  $hd  ", $an/$ad, "  ", ($hd > 0)? $hn/$hd : '---',  "  ", GenotypeSimulation::paircode_count_string($gpc_count), "\n";
         } elsif ($relationship eq 'AUNN') { # aunt/uncle - niece/nephew
            my $gp1_genotype = SimulatedGenotype->new_from_mafs($the_rng, $mafs, $gen, \$id);
            my $gp2_genotype = SimulatedGenotype->new_from_mafs($the_rng, $mafs, $gen, \$id);

            my $p1_genotype = SimulatedGenotype->new_offspring($the_rng, $gp1_genotype, $gp2_genotype, $gen, \$id); # p1, p2 are full siblings
            my $p2_genotype = SimulatedGenotype->new_offspring($the_rng, $gp1_genotype, $gp2_genotype, $gen, \$id);

            my $p3_genotype = SimulatedGenotype->new_from_mafs($the_rng, $mafs, $gen, \$id); # p3, p4 are unrelated
            #   my $p4_genotype = SimulatedGenotype->new_from_mafs($the_rng, $mafs, $gen, \$id);

            my $nn_genotype = SimulatedGenotype->new_offspring($the_rng, $p1_genotype, $p3_genotype, $gen, \$id); # offspring of p1 and p3, niece/nephew of p2
            push @genotype_objects, ($p2_genotype, $nn_genotype);
            #    my $c2_genotype = SimulatedGenotype->new_offspring($the_rng, $p2_genotype, $p4_genotype, $gen, \$id);
            #      my ($an, $ad, $hn, $hd, $gpc_count) = GenotypeSimulation::agmr_hgmr($p2_genotype, $nn_genotype, $gen, \$id);
            #      print "AUNN   $an  $ad   $hn  $hd  ", $an/$ad, "  ", ($hd > 0)? $hn/$hd : '---',  "  ", GenotypeSimulation::paircode_count_string($gpc_count), "\n";
         } elsif ($relationship eq 'GPGC') { # grandparent - grandchild
            my $gp1_genotype = SimulatedGenotype->new_from_mafs($the_rng, $mafs, $gen, \$id);
            my $gp2_genotype = SimulatedGenotype->new_from_mafs($the_rng, $mafs, $gen, \$id);

            my $p1_genotype = SimulatedGenotype->new_offspring($the_rng, $gp1_genotype, $gp2_genotype, $gen, \$id); # offspring of gp1, gp2
            my $p3_genotype = SimulatedGenotype->new_from_mafs($the_rng, $mafs, $gen, \$id); # other parent - unrelated
            #   my $p4_genotype = SimulatedGenotype->new_from_mafs($the_rng, $mafs, $gen, \$id);

            my $gc_genotype = SimulatedGenotype->new_offspring($the_rng, $p1_genotype, $p3_genotype, $gen, \$id); # offspring of p1 and p3, gc of gp1 & gp2
            push @genotype_objects, ($gp2_genotype, $gc_genotype);
            #    my $c2_genotype = SimulatedGenotype->new_offspring($the_rng, $p2_genotype, $p4_genotype, $gen, \$id);
            #     my ($an, $ad, $hn, $hd, $gpc_count) = GenotypeSimulation::agmr_hgmr($gp2_genotype, $gc_genotype, $gen, \$id);
            #     print "GPGC   $an  $ad   $hn  $hd  ", $an/$ad, "  ", ($hd > 0)? $hn/$hd : '---',  "  ", GenotypeSimulation::paircode_count_string($gpc_count), "\n";
         } elsif ($relationship eq 'HS') { # half-siblings
            my $p1_genotype = SimulatedGenotype->new_from_mafs($the_rng, $mafs, $gen, \$id);
            my $p2_genotype = SimulatedGenotype->new_from_mafs($the_rng, $mafs, $gen, \$id);
            my $p3_genotype = SimulatedGenotype->new_from_mafs($the_rng, $mafs, $gen, \$id);
            my $o12_genotype = SimulatedGenotype->new_offspring($the_rng, $p1_genotype, $p2_genotype, $gen, \$id);
            my $o23_genotype = SimulatedGenotype->new_offspring($the_rng, $p2_genotype, $p3_genotype, $gen, \$id);
        #    my ($an, $ad, $hn, $hd, $gpc_count) = GenotypeSimulation::agmr_hgmr($o12_genotype, $o23_genotype, $gen, \$id);
                  push @genotype_objects, ($o12_genotype, $o23_genotype);
            #      print "HS   $an  $ad   $hn  $hd  ", $an/$ad, "  ", ($hd > 0)? $hn/$hd : '---',  "  ", GenotypeSimulation::paircode_count_string($gpc_count), "\n";
         } elsif ($relationship eq 'HAUNN') { # half-aunt/uncle - half-niece/nephew
            my $p1_genotype = SimulatedGenotype->new_from_mafs($the_rng, $mafs, $gen, \$id);
            my $p2_genotype = SimulatedGenotype->new_from_mafs($the_rng, $mafs, $gen, \$id);
            my $p3_genotype = SimulatedGenotype->new_from_mafs($the_rng, $mafs, $gen, \$id);

            my $o12_genotype = SimulatedGenotype->new_offspring($the_rng, $p1_genotype, $p2_genotype, $gen, \$id);
            my $o23_genotype = SimulatedGenotype->new_offspring($the_rng, $p2_genotype, $p3_genotype, $gen, \$id);

            my $p4_genotype = SimulatedGenotype->new_from_mafs($the_rng, $mafs, $gen, \$id);
            my $hnn_genotype = SimulatedGenotype->new_offspring($the_rng, $o12_genotype, $p4_genotype, $gen, \$id);
            push @genotype_objects, ($hnn_genotype, $o23_genotype);
            #      my ($an, $ad, $hn, $hd, $gpc_count) = GenotypeSimulation::agmr_hgmr($hnn_genotype, $o23_genotype, $gen, \$id);
            #      print "HAUNN   $an  $ad   $hn  $hd  ", $an/$ad, "  ", ($hd > 0)? $hn/$hd : '---',  "  ", GenotypeSimulation::paircode_count_string($gpc_count), "\n";
         } elsif ($relationship eq 'UN') { # unrelated
            my $genotype1 = SimulatedGenotype->new_from_mafs($the_rng, $mafs, $gen, \$id);
            my $genotype2 = SimulatedGenotype->new_from_mafs($the_rng, $mafs, $gen, \$id);
            push @genotype_objects, ($genotype1, $genotype2);
            #      my ($an, $ad, $hn, $hd, $gpc_count) = GenotypeSimulation::agmr_hgmr($genotype1, $genotype2);
            #      print "UN   $an  $ad   $hn  $hd  ", $an/$ad, "  ", ($hd > 0)? $hn/$hd : '---',  "  ", GenotypeSimulation::paircode_count_string($gpc_count), "\n";
         } elsif ($relationship eq 'GP4ID') { # all 4 grandparents are the same individual, parents distinct.
             my $gp_genotype = SimulatedGenotype->new_from_mafs($the_rng, $mafs, $gen, \$id);
             my $p1_genotype = SimulatedGenotype->new_offspring($the_rng, $gp_genotype, $gp_genotype, $gen, \$id);
             my $p2_genotype = SimulatedGenotype->new_offspring($the_rng, $gp_genotype, $gp_genotype, $gen, \$id);
             my $o_genotype = SimulatedGenotype->new_offspring($the_rng, $p1_genotype, $p2_genotype, $gen, \$id);
             push @genotype_objects, ($gp_genotype, $o_genotype);
          }elsif ($relationship eq 'GP1122') { # 
             my $gp1_genotype = SimulatedGenotype->new_from_mafs($the_rng, $mafs, $gen, \$id);
             my $gp2_genotype = SimulatedGenotype->new_from_mafs($the_rng, $mafs, $gen, \$id);
             my $p1_genotype = SimulatedGenotype->new_offspring($the_rng, $gp1_genotype, $gp1_genotype, $gen, \$id);
             my $p2_genotype = SimulatedGenotype->new_offspring($the_rng, $gp2_genotype, $gp2_genotype, $gen, \$id);
             my $o_genotype = SimulatedGenotype->new_offspring($the_rng, $p1_genotype, $p2_genotype, $gen, \$id);
             push @genotype_objects, ($gp1_genotype, $o_genotype);
          }elsif ($relationship eq 'GP1212') { # parents are full siblings
             my $gp1_genotype = SimulatedGenotype->new_from_mafs($the_rng, $mafs, $gen, \$id);
             my $gp2_genotype = SimulatedGenotype->new_from_mafs($the_rng, $mafs, $gen, \$id);
             my $p1_genotype = SimulatedGenotype->new_offspring($the_rng, $gp1_genotype, $gp2_genotype, $gen, \$id);
             my $p2_genotype = SimulatedGenotype->new_offspring($the_rng, $gp1_genotype, $gp2_genotype, $gen, \$id);
             my $o_genotype = SimulatedGenotype->new_offspring($the_rng, $p1_genotype, $p2_genotype, $gen, \$id);
             push @genotype_objects, ($gp1_genotype, $o_genotype);
          }elsif ($relationship eq 'GP1112') { # three gp's are same individual - relation between gp 1 and gc.
             my $gp1_genotype = SimulatedGenotype->new_from_mafs($the_rng, $mafs, $gen, \$id);
             my $gp2_genotype = SimulatedGenotype->new_from_mafs($the_rng, $mafs, $gen, \$id);
             my $p1_genotype = SimulatedGenotype->new_offspring($the_rng, $gp1_genotype, $gp1_genotype, $gen, \$id);
             my $p2_genotype = SimulatedGenotype->new_offspring($the_rng, $gp1_genotype, $gp2_genotype, $gen, \$id);
             my $o_genotype = SimulatedGenotype->new_offspring($the_rng, $p1_genotype, $p2_genotype, $gen, \$id);
             push @genotype_objects, ($gp1_genotype, $o_genotype);
          }elsif ($relationship eq 'GP1111P22') { # all 4 grandparents are the one individual, and both parents are one indiv.
             my $gp_genotype = SimulatedGenotype->new_from_mafs($the_rng, $mafs, $gen, \$id);
             my $p_genotype = SimulatedGenotype->new_offspring($the_rng, $gp_genotype, $gp_genotype, $gen, \$id);
             my $o_genotype = SimulatedGenotype->new_offspring($the_rng, $p_genotype, $p_genotype, $gen, \$id);
             push @genotype_objects, ($gp_genotype, $o_genotype);
          }elsif ($relationship eq 'X00') { # rel. between two distinct indivs whose 4 parents are all distinct
             my $pa1_genotype = SimulatedGenotype->new_from_mafs($the_rng, $mafs, $gen, \$id);
             my $pa2_genotype = SimulatedGenotype->new_from_mafs($the_rng, $mafs, $gen, \$id);
             my $pb1_genotype = SimulatedGenotype->new_from_mafs($the_rng, $mafs, $gen, \$id);
             my $pb2_genotype = SimulatedGenotype->new_from_mafs($the_rng, $mafs, $gen, \$id);

             my $a_genotype = SimulatedGenotype->new_offspring($the_rng, $pa1_genotype, $pa2_genotype, $gen, \$id);
             my $b_genotype = SimulatedGenotype->new_offspring($the_rng, $pb1_genotype, $pb2_genotype, $gen, \$id);
             push @genotype_objects, ($a_genotype, $b_genotype);
          }elsif ($relationship eq 'X01') { # rel. between two distinct indivs. 3 distinct parents, a1=b1 half-sibling
             my $pa1_genotype = SimulatedGenotype->new_from_mafs($the_rng, $mafs, $gen, \$id);
             my $pa2_genotype = SimulatedGenotype->new_from_mafs($the_rng, $mafs, $gen, \$id);
             my $pb1_genotype = $pa1_genotype; # SimulatedGenotype->new_from_mafs($the_rng, $mafs, $gen, \$id);
             my $pb2_genotype = SimulatedGenotype->new_from_mafs($the_rng, $mafs, $gen, \$id);

             my $a_genotype = SimulatedGenotype->new_offspring($the_rng, $pa1_genotype, $pa2_genotype, $gen, \$id);
             my $b_genotype = SimulatedGenotype->new_offspring($the_rng, $pb1_genotype, $pb2_genotype, $gen, \$id);
             push @genotype_objects, ($a_genotype, $b_genotype);
          }elsif ($relationship eq 'X02') { # rel. between two distinct indivs. 2 distinct parents, a1=b1, a2=b2 full sibling
             my $pa1_genotype = SimulatedGenotype->new_from_mafs($the_rng, $mafs, $gen, \$id);
             my $pa2_genotype = SimulatedGenotype->new_from_mafs($the_rng, $mafs, $gen, \$id);
             my $pb1_genotype = $pa1_genotype; # SimulatedGenotype->new_from_mafs($the_rng, $mafs, $gen, \$id);
             my $pb2_genotype = $pa2_genotype; # SimulatedGenotype->new_from_mafs($the_rng, $mafs, $gen, \$id);

             my $a_genotype = SimulatedGenotype->new_offspring($the_rng, $pa1_genotype, $pa2_genotype, $gen, \$id);
             my $b_genotype = SimulatedGenotype->new_offspring($the_rng, $pb1_genotype, $pb2_genotype, $gen, \$id);
             push @genotype_objects, ($a_genotype, $b_genotype)
          }elsif ($relationship eq 'X10') { # rel. between two distinct indivs. 3 distinct parents. a1, a2 are identical
             my $pa1_genotype = SimulatedGenotype->new_from_mafs($the_rng, $mafs, $gen, \$id);
             my $pa2_genotype = $pa1_genotype; # SimulatedGenotype->new_from_mafs($the_rng, $mafs, $gen, \$id);
             my $pb1_genotype = SimulatedGenotype->new_from_mafs($the_rng, $mafs, $gen, \$id);
             my $pb2_genotype = SimulatedGenotype->new_from_mafs($the_rng, $mafs, $gen, \$id);

             my $a_genotype = SimulatedGenotype->new_offspring($the_rng, $pa1_genotype, $pa2_genotype, $gen, \$id);
             my $b_genotype = SimulatedGenotype->new_offspring($the_rng, $pb1_genotype, $pb2_genotype, $gen, \$id);
             push @genotype_objects, ($a_genotype, $b_genotype);
          }elsif ($relationship eq 'X20') { # rel. between two distinct indivs. 2 distinct parent, a1=a2, b1=b2
             my $pa1_genotype = SimulatedGenotype->new_from_mafs($the_rng, $mafs, $gen, \$id);
             my $pa2_genotype = $pa1_genotype; # SimulatedGenotype->new_from_mafs($the_rng, $mafs, $gen, \$id);
             my $pb1_genotype = SimulatedGenotype->new_from_mafs($the_rng, $mafs, $gen, \$id);
             my $pb2_genotype = $pb1_genotype; # SimulatedGenotype->new_from_mafs($the_rng, $mafs, $gen, \$id);

             my $a_genotype = SimulatedGenotype->new_offspring($the_rng, $pa1_genotype, $pa2_genotype, $gen, \$id);
             my $b_genotype = SimulatedGenotype->new_offspring($the_rng, $pb1_genotype, $pb2_genotype, $gen, \$id);
             push @genotype_objects, ($a_genotype, $b_genotype);
          }elsif ($relationship eq 'X12') { # rel. between two distinct indivs. 3 identical parents, a1=a2=b1
             my $pa1_genotype = SimulatedGenotype->new_from_mafs($the_rng, $mafs, $gen, \$id);
             my $pa2_genotype = $pa1_genotype; # SimulatedGenotype->new_from_mafs($the_rng, $mafs, $gen, \$id);
             my $pb1_genotype = $pa1_genotype; # SimulatedGenotype->new_from_mafs($the_rng, $mafs, $gen, \$id);
             my $pb2_genotype = SimulatedGenotype->new_from_mafs($the_rng, $mafs, $gen, \$id);

             my $a_genotype = SimulatedGenotype->new_offspring($the_rng, $pa1_genotype, $pa2_genotype, $gen, \$id);
             my $b_genotype = SimulatedGenotype->new_offspring($the_rng, $pb1_genotype, $pb2_genotype, $gen, \$id);
             push @genotype_objects, ($a_genotype, $b_genotype);
          }elsif ($relationship eq 'X24') { # rel. between two distinct indivs. 2 distinct parent, a1=a2, b1=b2
             my $pa1_genotype = SimulatedGenotype->new_from_mafs($the_rng, $mafs, $gen, \$id);
             my $pa2_genotype = $pa1_genotype; # SimulatedGenotype->new_from_mafs($the_rng, $mafs, $gen, \$id);
             my $pb1_genotype = $pa1_genotype;
             my $pb2_genotype = $pa1_genotype;

             my $a_genotype = SimulatedGenotype->new_offspring($the_rng, $pa1_genotype, $pa2_genotype, $gen, \$id);
             my $b_genotype = SimulatedGenotype->new_offspring($the_rng, $pb1_genotype, $pb2_genotype, $gen, \$id);
             push @genotype_objects, ($a_genotype, $b_genotype);
          }
      }    # end loop over simulated pairs of one type (relationship).
      #  print "n gobjs: ", scalar @genotype_objects, "\n";
      # ############## output fasta ##########################
      print "# relationship: $relationship   mafspec $mafspec   nsnps: $n_snps   npairs: $n_pairs \n";
      for my $gobj (@genotype_objects) {
         #    print STDERR "ref gobj: ", ref $gobj, "\n";
         my ($idline, $seqline) = $gobj->fasta($fasta_character_set);
         print ">$idline\n", "$seqline\n";
      }

   }                            # end loop over relationships

} # ############################# end of main #########################################


 #   while (my ($i, $g1) = each @genotypes) {
#       for (my $j = $i+1; $j < scalar @genotypes; $j++) {
#          #   print STDERR "$i $j \n";
#          my $g2 = $genotypes[$j];
#          #    GenotypeSimulation::print_genotypes_in_columns($g1, $g2);
#          my ($an, $ad, $hn, $hd, $gpc_count) = GenotypeSimulation::agmr_hgmr($g1, $g2);
#          my $pair_id = $i . '_' . $j;
# print "$pair_id   $an  $ad   $hn  $hd  ", $an/$ad, "  ", ($hd > 0)? $hn/$hd : '---',  "  ", GenotypeSimulation::paircode_count_string($gpc_count), "\n";
#       }
#    }


# }

