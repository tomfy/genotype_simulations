#!/usr/bin/perl -w
use strict;
use lib '/usr/local/share/perl/5.14.2';
use List::Util qw(min max sum);
use Graphics::GnuplotIF qw(GnuplotIF);
use Math::GSL::SF  qw( :all );
use Math::GSL::RNG  qw( :all );
use Math::GSL::Randist  qw( :all );
use Math::GSL::CDF  qw( :all );
use Getopt::Long;

# globals
my $the_rng = Math::GSL::RNG->new();

{
   my $n_snps = 10;             # number of snps in each sample.
   my $relationship = 'PO'; # PO: parent-offspring. Or FS (full sib), HS (half sib), etc. ...
   my $n_reps = 100; # number or replicates (i.e. simulated sample pairs, each with $n_snps in both samples
   my $maf = 'delta,0.25';      # minor allele frequency
 
   GetOptions(
              'nsnps=i' => \$n_snps, # 
              'relationship=s' => \$relationship,
              'nreps=i' => \$n_reps,
              'maf|minor_allele_frequence=s' => \$maf, # 0 < $maf <= 0.5
             );
   for (1..$n_reps) {

      if($relationship eq 'DPO'){ # both parents are same individual
my $p1_genotype = draw_genotype_from_population($maf, $n_snps);
  my $o_genotype = draw_offspring_genotype($p1_genotype, $p1_genotype);
 my ($an, $ad, $hn, $hd, $gpc_count) = agmr_hgmr($p1_genotype, $o_genotype);
         print "DPO   $an  $ad   $hn  $hd  ", $an/$ad, "  ", ($hd > 0)? $hn/$hd : '---',  "  ", paircode_count_string($gpc_count), "\n";
}elsif($relationship eq 'PO') { # parent and offspring
         my $p1_genotype = draw_genotype_from_population($maf, $n_snps);
         my $p2_genotype = draw_genotype_from_population($maf, $n_snps);
         my $o_genotype = draw_offspring_genotype($p1_genotype, $p2_genotype);

         #  print_genotype($p1_genotype, $p2_genotype, $o_genotype);
         my ($an, $ad, $hn, $hd, $gpc_count) = agmr_hgmr($p1_genotype, $o_genotype);
         print "PO   $an  $ad   $hn  $hd  ", $an/$ad, "  ", ($hd > 0)? $hn/$hd : '---', "  ", paircode_count_string($gpc_count), "\n";
      } elsif ($relationship eq 'FS') { # full siblings
         my $p1_genotype = draw_genotype_from_population($maf, $n_snps);
         my $p2_genotype = draw_genotype_from_population($maf, $n_snps);
         my $o1_genotype = draw_offspring_genotype($p1_genotype, $p2_genotype);
         my $o2_genotype = draw_offspring_genotype($p1_genotype, $p2_genotype);
         my ($an, $ad, $hn, $hd, $gpc_count) = agmr_hgmr($o1_genotype, $o2_genotype);
         print "FS   $an  $ad   $hn  $hd  ", $an/$ad, "  ", ($hd > 0)? $hn/$hd : '---',  "  ", paircode_count_string($gpc_count), "\n";
      } elsif ($relationship eq 'FC') { # first cousins
         my $gp1_genotype = draw_genotype_from_population($maf, $n_snps);
         my $gp2_genotype = draw_genotype_from_population($maf, $n_snps);

         my $p1_genotype = draw_offspring_genotype($gp1_genotype, $gp2_genotype); # p1, p2 are full sibling
         my $p2_genotype = draw_offspring_genotype($gp1_genotype, $gp2_genotype);

         my $p3_genotype = draw_genotype_from_population($maf, $n_snps); # p3, p4 are unrelated
         my $p4_genotype = draw_genotype_from_population($maf, $n_snps);

         my $c1_genotype = draw_offspring_genotype($p1_genotype, $p3_genotype);
         my $c2_genotype = draw_offspring_genotype($p2_genotype, $p4_genotype);
         my ($an, $ad, $hn, $hd, $gpc_count) = agmr_hgmr($c1_genotype, $c2_genotype);
         print "FC   $an  $ad   $hn  $hd  ", $an/$ad, "  ", ($hd > 0)? $hn/$hd : '---',  "  ", paircode_count_string($gpc_count), "\n";
      } elsif ($relationship eq 'AUNN') { # aunt/uncle - niece/nephew
           my $gp1_genotype = draw_genotype_from_population($maf, $n_snps);
         my $gp2_genotype = draw_genotype_from_population($maf, $n_snps);

         my $p1_genotype = draw_offspring_genotype($gp1_genotype, $gp2_genotype); # p1, p2 are full siblings
         my $p2_genotype = draw_offspring_genotype($gp1_genotype, $gp2_genotype);

         my $p3_genotype = draw_genotype_from_population($maf, $n_snps); # p3, p4 are unrelated
         #   my $p4_genotype = draw_genotype_from_population($maf, $n_snps);

         my $nn_genotype = draw_offspring_genotype($p1_genotype, $p3_genotype); # offspring of p1 and p3, niece/nephew of p2
         #    my $c2_genotype = draw_offspring_genotype($p2_genotype, $p4_genotype);
         my ($an, $ad, $hn, $hd, $gpc_count) = agmr_hgmr($p2_genotype, $nn_genotype);
         print "AUNN   $an  $ad   $hn  $hd  ", $an/$ad, "  ", ($hd > 0)? $hn/$hd : '---',  "  ", paircode_count_string($gpc_count), "\n";
        }elsif($relationship eq 'GPGC'){ # grandparent - grandchild
           my $gp1_genotype = draw_genotype_from_population($maf, $n_snps);
           my $gp2_genotype = draw_genotype_from_population($maf, $n_snps);

         my $p1_genotype = draw_offspring_genotype($gp1_genotype, $gp2_genotype); # offspring of gp1, gp2
         my $p3_genotype = draw_genotype_from_population($maf, $n_snps); # other parent - unrelated
         #   my $p4_genotype = draw_genotype_from_population($maf, $n_snps);

         my $gc_genotype = draw_offspring_genotype($p1_genotype, $p3_genotype); # offspring of p1 and p3, gc of gp1 & gp2
         #    my $c2_genotype = draw_offspring_genotype($p2_genotype, $p4_genotype);
         my ($an, $ad, $hn, $hd, $gpc_count) = agmr_hgmr($gp2_genotype, $gc_genotype);
         print "GPGC   $an  $ad   $hn  $hd  ", $an/$ad, "  ", ($hd > 0)? $hn/$hd : '---',  "  ", paircode_count_string($gpc_count), "\n";
        }elsif ($relationship eq 'HS') { # half-siblings
         my $p1_genotype = draw_genotype_from_population($maf, $n_snps);
         my $p2_genotype = draw_genotype_from_population($maf, $n_snps);
         my $p3_genotype = draw_genotype_from_population($maf, $n_snps);
         my $o12_genotype = draw_offspring_genotype($p1_genotype, $p2_genotype);
         my $o23_genotype = draw_offspring_genotype($p2_genotype, $p3_genotype);
         my ($an, $ad, $hn, $hd, $gpc_count) = agmr_hgmr($o12_genotype, $o23_genotype);
         print "HS   $an  $ad   $hn  $hd  ", $an/$ad, "  ", ($hd > 0)? $hn/$hd : '---',  "  ", paircode_count_string($gpc_count), "\n";
      } elsif ($relationship eq 'HAUNN') { # half-aunt/uncle - half-niece/nephew
         my $p1_genotype = draw_genotype_from_population($maf, $n_snps);
         my $p2_genotype = draw_genotype_from_population($maf, $n_snps);
         my $p3_genotype = draw_genotype_from_population($maf, $n_snps);

         my $o12_genotype = draw_offspring_genotype($p1_genotype, $p2_genotype);
         my $o23_genotype = draw_offspring_genotype($p2_genotype, $p3_genotype);

         my $p4_genotype = draw_genotype_from_population($maf, $n_snps);
         my $hnn_genotype = draw_offspring_genotype($o12_genotype, $p4_genotype);
         my ($an, $ad, $hn, $hd, $gpc_count) = agmr_hgmr($hnn_genotype, $o23_genotype);
         print "HAUNN   $an  $ad   $hn  $hd  ", $an/$ad, "  ", ($hd > 0)? $hn/$hd : '---',  "  ", paircode_count_string($gpc_count), "\n";
      }elsif($relationship eq 'UN'){ # unrelated
         my $genotype1 = draw_genotype_from_population($maf, $n_snps);
         my $genotype2 = draw_genotype_from_population($maf, $n_snps);
         my ($an, $ad, $hn, $hd, $gpc_count) = agmr_hgmr($genotype1, $genotype2);
         print "UN   $an  $ad   $hn  $hd  ", $an/$ad, "  ", ($hd > 0)? $hn/$hd : '---',  "  ", paircode_count_string($gpc_count), "\n";
      }
   }
}

sub draw_alleles{
   my $maf_distribution = shift;
   my $n_alleles_to_draw = shift;
   my ($distrib, $param1, $param2) = split(',', $maf_distribution);
   my @alleles = ();
   my $maf = undef;
   if ($distrib eq 'delta') { # delta function at param1 (param2 ignored if present)
      for (1..$n_alleles_to_draw) {
         $maf = $param1;
         push @alleles, (gsl_rng_uniform($the_rng->raw() ) < $maf)? 'a' : 'A';
         push @alleles, (gsl_rng_uniform($the_rng->raw() ) < $maf)? 'a' : 'A';
      }
   } elsif ($distrib eq 'flat') {
      my ($min, $max) = ($param1 // 0.0, $param2 // 0.5);
      for (1..$n_alleles_to_draw) {
         $maf = $min + gsl_rng_uniform($the_rng->raw()) * ($max - $min);
         push @alleles, (gsl_rng_uniform($the_rng->raw() ) < $maf)? 'a' : 'A';
         push @alleles, (gsl_rng_uniform($the_rng->raw() ) < $maf)? 'a' : 'A';
         #   print STDERR "maf: $maf \n";
      }
   } elsif ($distrib eq 'inv') { # density propto 1/$maf (0 outside of [$min, $max]
      my ($min, $max) = ($param1 // 0.0, $param2 // 0.5);
      for (1..$n_alleles_to_draw) {
         $maf =  $min * exp(gsl_rng_uniform($the_rng->raw()) * log($max/$min));
         # $min + gsl_rng_uniform($the_rng->raw()) * ($max - $min);
         push @alleles, (gsl_rng_uniform($the_rng->raw() ) < $maf)? 'a' : 'A';
         push @alleles, (gsl_rng_uniform($the_rng->raw() ) < $maf)? 'a' : 'A';
         #   print STDERR "maf: $maf \n";
      }
   } else {
      die "Unknown minor allele frequence distribution: $distrib \n";
   }
   # print STDERR "number of alleles drawn: ", scalar @alleles, "\n";
   return \@alleles;
}

sub draw_genotype_from_population{
   # draw ('a','a'), ('a','A'), etc.
   my $maf = shift;             # minor allele frequency
   my $n_snps = shift;
   my @genotypes = ();
   my $alleles = draw_alleles($maf, $n_snps);
   for my $i_snp (0..$n_snps-1) {
      my $gt = [$alleles->[2*$i_snp], $alleles->[2*$i_snp + 1]];
      # for (1..2) {
      #    push  @$gt, 
      #      draw_allele($maf); # 
      #    # ( gsl_rng_uniform($the_rng->raw() ) < 0.3)? 'a' : 'A';    
      # }
      push @genotypes, $gt;
   }
   return \@genotypes;
}

sub draw_offspring_genotype{
   my $pg1 = shift;             # parent genotype 1
   my $pg2 = shift;             # parent genotype 2
   my @genotypes = ();
   die "number of snps is different in the 2 parents. bye. \n" if(scalar @$pg1 != scalar @$pg2);
   # print STDERR "n snps: ", scalar @$pg1, "\n";
   while (my ($i, $g1) = each @$pg1) {
      #   print "i: $i \n";
      my $g2 = $pg2->[$i];
      my $og = [];              # offspring genotype at this snp
      push @$og, ( gsl_rng_uniform($the_rng->raw() ) < 0.5)? $g1->[0] : $g1->[1];
      push @$og, ( gsl_rng_uniform($the_rng->raw() ) < 0.5)? $g2->[0] : $g2->[1];
      #   print "offspring snp:  $i  ", join(",", @$og), "\n";
      push @genotypes, $og;
   }
   return \@genotypes;
}


sub agmr_hgmr{
   my $g1 = shift;
   my $g2 = shift;
   die "number of snps is different in the 2 genotypes. bye. \n" if(scalar @$g1 != scalar @$g2);
   my ($an, $ad, $hn, $hd) = (0, 0, 0, 0);
# encode genotype: 0: AA, 1: aA, 2: aa  i.e. a count of the number of minor alleles present.
# encode the pair of genotypes as 00, 01, 10, 02, 20, etc.
my %paircode_count = ();
   while (my($i, $s1) = each $g1) {
      my $s2 = $g2->[$i]; # corresponding snp from sample 2
      my $ss1 = join("", @$s1);
      my $ss2 = join("", @$s2);
      if ($ss1 eq 'aa') {
         if ($ss2 eq 'aa') {
            $ad++; $hd++;
            $paircode_count{'22'}++; # aa-aa -> 22
         } elsif ($ss2 eq 'AA') {
            $ad++; $an++; $hd++; $hn++;
            $paircode_count{'20'}++; # aa-AA -> 20
         } else {               # s2 is heterozygous
            $ad++; $an++;
            $paircode_count{'21'}++; # aa-Aa -> 21
         }
      } elsif ($ss1 eq 'AA') {
         if ($ss2 eq 'aa') {
            $ad++; $an++; $hd++; $hn++;
            $paircode_count{'02'}++;
         } elsif ($ss2 eq 'AA') {
            $ad++; $hd++;
            $paircode_count{'00'}++;
         } else {               # s2 is heterozygous
            $ad++; $an++;
            $paircode_count{'01'}++;
         }
      } else {                  # s1 is heterozygous
         if ($ss2 eq 'aa') {
            $ad++; $an++;
            $paircode_count{'12'}++;
         } elsif ($ss2 eq 'AA') {
            $ad++; $an++;
            $paircode_count{'10'}++;
         } else {               # s2 is heterozygous
            $ad++;
            $paircode_count{'11'}++;
         }
      }
   }
   return ($an, $ad, $hn, $hd, \%paircode_count); 
}


sub print_genotype{
   my @gts = @_;
   for my $i ( 0.. scalar @{$gts[0]}-1 ) {
      for my $agt (@gts) {
         my $snp = $agt->[$i];
         print join("", @$snp), "   ";
      }
      print "\n";
   }
}
  sub paircode_count_string{
     my $pc_c = shift;
     my $str = ''; # "[ ";
     for my $g1 ('0','1','2') {
        for my $g2 ('0','1','2') {
           my $g1g2 = $g1 . $g2;
           $str .= sprintf("%4i ", $pc_c->{$g1g2} // 0);
        }
        $str .= " "; # "] [";
     }
     #$str .= "]";
     return $str;
  }
