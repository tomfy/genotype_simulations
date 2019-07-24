package GenotypeSim;
use strict;
use 5.012;
use List::Util qw (min max sum);
use Math::GSL::SF  qw( :all );
use Math::GSL::RNG  qw( :all );
use Math::GSL::Randist  qw( :all );
use Math::GSL::CDF  qw( :all );


sub draw_alleles{
   my $the_rng = shift;
   my $maf_distribution = shift; # describes distribution of minor allele frequencies. e.g. 'delta,0.25' (all mafs are 0.25), or 'inv,0.15,0.5' min maf is 0.15, max is 0.5, between prob. is propto 1/maf.
   my $n_SNPs = shift;           # for each SNP, draw 2 alleles

   # print STDERR "maf distrib:  $maf_distribution \n";
   my ($distrib, $param1, $param2) = split(',', $maf_distribution);
   my @alleles = ();
   my $maf = undef;
   if ($distrib eq 'delta') { # delta function at param1 (param2 ignored if present)
      for (1..$n_SNPs) {
         $maf = $param1;
         push @alleles, (gsl_rng_uniform($the_rng->raw() ) < $maf)? 'a' : 'A';
         push @alleles, (gsl_rng_uniform($the_rng->raw() ) < $maf)? 'a' : 'A';
      }
   } elsif ($distrib eq 'flat') {
      my ($min, $max) = ($param1 // 0.0, $param2 // 0.5);
      for (1..$n_SNPs) {
         $maf = $min + gsl_rng_uniform($the_rng->raw()) * ($max - $min);
         push @alleles, (gsl_rng_uniform($the_rng->raw() ) < $maf)? 'a' : 'A';
         push @alleles, (gsl_rng_uniform($the_rng->raw() ) < $maf)? 'a' : 'A';
         #   print STDERR "maf: $maf \n";
      }
   } elsif ($distrib eq 'inv') { # density propto 1/$maf (0 outside of [$min, $max]
      my ($min, $max) = ($param1 // 0.0, $param2 // 0.5);
      for (1..$n_SNPs) {
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
   my $the_rng = shift;
   my $maf = shift;             # minor allele frequency
   my $n_snps = shift;
   # print STDERR "$the_rng ; $maf ;  $n_snps \n";
   my @genotypes = ();
   my $alleles = draw_alleles($the_rng, $maf, $n_snps); # a array ref holding 2*$n_snps alleles.
   for my $i_snp (0..$n_snps-1) {
      my $gt = [$alleles->[2*$i_snp], $alleles->[2*$i_snp + 1]];
      push @genotypes, $gt;
   }
   return \@genotypes;
}

sub draw_offspring_genotype{
   my $the_rng = shift; 
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
      my $s2 = $g2->[$i];       # corresponding snp from sample 2
      my $ss1 = join("", @$s1);
      my $ss2 = join("", @$s2);
      if ($ss1 eq 'aa') {
         if ($ss2 eq 'aa') {
            $ad++; $hd++;
            $paircode_count{'22'}++; # aa-aa -> 22
         } elsif ($ss2 eq 'AA') {
            $ad++; $an++; $hd++; $hn++;
            $paircode_count{'20'}++; # aa-AA -> 20
         } else {                    # s2 is heterozygous
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
   for my $i ( 0 .. scalar @{$gts[0]}-1 ) {
      for my $agt (@gts) {
         my $snp = $agt->[$i];
         print join("", @$snp), "   ";
      }
      print "\n";
   }
}


  sub paircode_count_string{
     my $pc_c = shift;
     my $str = '';              # "[ ";
     for my $g1 ('0','1','2') {
        for my $g2 ('0','1','2') {
           my $g1g2 = $g1 . $g2;
           $str .= sprintf("%4i ", $pc_c->{$g1g2} // 0);
        }
        $str .= " ";            # "] [";
     }
     #$str .= "]";
     return $str;
  }

# end of package
1; 
