package Genotype;
use Math::GSL::RNG  qw( :all );


sub new_from_population{
   my $class = shift;
   #   my $args = shift;
   my $self  = bless {}, $class;
   my $the_rng = shift;
   my $mafs = shift;            # array ref of maf values
   my $n_snps = scalar @$mafs;
   my $generation = shift;
   my $id = shift;
   # print STDERR "$the_rng ; $maf ;  $n_snps \n";
   
   #  print STDERR "ZZZ: ", join(", ", @$mafs), "\n";
   #  $self->set_mafs($mafs);
   $self->set_rng($the_rng);
   $self->set_generation($generation);
   $self->set_id($id);
   $self->set_parents('X', 'X');

   my $alleles = $self->draw_alleles($mafs); # an array ref holding 2*$n_snps alleles.
   my @genotypes = ();
   for my $i_snp (0..$n_snps-1) {
      my $gt =  $alleles->[2*$i_snp] . $alleles->[2*$i_snp + 1]; # [$alleles->[2*$i_snp], $alleles->[2*$i_snp + 1]];
      push @genotypes, $gt;
   }
   #   print "in new_from_pop. [", ref \@genotypes, "] \n";
   $self->set_genotypes(\@genotypes);
   #   print "in new_from_pop. [", ref $self->get_genotypes(), "] \n";
   return $self;
}

sub new_offspring{
   my $class = shift;
   my $self = bless {}, $class;
   my $the_rng = shift;
   my $pgobj1 = shift;          # parent genotype obj. 1
   my $pgobj2 = shift;          # parent genotype obj. 2
   my $generation = shift;
   my $id = shift;
   my @genotypes = ();
   my $pg1 = $pgobj1->get_genotypes();
   my $pg2 = $pgobj2->get_genotypes();
   die "number of snps is different in the 2 parents. bye. \n" if(scalar @$pg1 != scalar @$pg2);
   # print STDERR "n snps: ", scalar @$pg1, "\n";
   $self->set_rng($the_rng);
   $self->set_generation($generation);
   $self->set_id($id);
   $self->set_parents($pgobj1->get_id(), $pgobj2->get_id());

   while (my ($i, $g1) = each @$pg1) {
      #   print "i: $i \n";
      my $g2 = $pg2->[$i];
      # my $og = [];              # offspring genotype at this snp
      # print STDERR "AAA:  [", ref $g1, "]\n"; 
      # push @$og, ( gsl_rng_uniform($the_rng->raw() ) < 0.5)? $g1->[0] : $g1->[1];
      # push @$og, ( gsl_rng_uniform($the_rng->raw() ) < 0.5)? $g2->[0] : $g2->[1];
      # #   print "offspring snp:  $i  ", join(",", @$og), "\n";

      my $og = '';
      my ($a,$b) = split('', $g1);
      $og .= ( gsl_rng_uniform($the_rng->raw() ) < 0.5)? $a : $b;
      ($a,$b) = split('', $g2);
      $og .= ( gsl_rng_uniform($the_rng->raw() ) < 0.5)? $a : $b;
      push @genotypes, $og;
   }
   $self->set_genotypes(\@genotypes);
   return $self;
}

sub draw_alleles{
   my $self = shift;
   my $mafs = shift;            # array ref, $n_SNPa.
   #   my $n_SNPs = scalar @$mafs;  # for each SNP, draw 2 alleles
   my $the_rng = $self->get_rng();

   my @alleles = ();
   for my $amaf (@$mafs) {
      #   print STDERR "amaf: $amaf \n";
      push @alleles, ((gsl_rng_uniform($the_rng->raw() ) < $amaf)? 'a' : 'A');
      push @alleles, ((gsl_rng_uniform($the_rng->raw() ) < $amaf)? 'a' : 'A');
   }
   # print STDERR "number of alleles drawn: ", scalar @alleles, "\n";
   #  print STDERR "alleles: ", join(", ", @alleles), "\n";
   return \@alleles;
}

sub genotype_string{
   my $self = shift;
   my $gt = $self->get_genotypes();
   # print "ref gt: ", ref $gt, "  ", scalar @$gt, "\n";
   my $gstring = '';
   while (my($i, $snp) = each @$gt) {
      # my $snp = $s; # join('', @$s);
      $snp = 'aA' if($snp eq 'Aa');
      $gstring .= $snp;
   }
   return $gstring;
}

sub set_rng{
   my $self = shift;
   my $rng = shift;
   $self->{rng} = $rng;
}

sub get_rng{
   my $self = shift;
   return $self->{rng};
}

sub set_genotypes{
   my $self = shift;
   my $genotypes = shift;
   $self->{genotypes} = $genotypes;
}

sub get_genotypes{
   my $self = shift;
   return $self->{genotypes};
}

sub set_generation{
   my $self = shift;
   my $generation = shift;
   $self->{generation} = $generation;
}

sub get_generation{
   my $self = shift;
   return $self->{generation};
}

sub set_id{
   my $self = shift;
   my $id = shift;
   $self->{id} = $id;
}

sub get_id{
   my $self = shift;
   return $self->{id};
}

sub set_parents{
   my $self = shift;
   my @parent_ids = @_;
   $self->{parents} = '[' . $parent_ids[0] . ']' . '[' . $parent_ids[1] . ']';
}

sub get_parents{
   my $self = shift;
   return $self->{parents};
}


1;

