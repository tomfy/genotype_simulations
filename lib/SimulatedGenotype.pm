package SimulatedGenotype;
use Math::GSL::RNG  qw( :all );


sub new_from_mafs{
   my $class = shift;
   #   my $args = shift;
   my $self  = bless {}, $class;
   my $the_rng = shift;
   my $mafs = shift;            # array ref of maf values
   my $n_snps = scalar @$mafs;
   my $generation = shift;
   my $id = shift; # ref to scalar 
my $max_pedigree_depth = shift // 2; # 2 -> keep just back to grandparents, etc.
   # print STDERR "$the_rng ; $maf ;  $n_snps \n";

   #  print STDERR "ZZZ: ", join(", ", @$mafs), "\n";
   #  $self->set_mafs($mafs);
   $self->set_rng($the_rng);
   $self->set_generation($generation);
   $self->set_id($$id);
   $self->set_parents('X', 'X');
   my $pedigree = '()' . $$id;
   $self->set_pedigree($pedigree);
   $self->set_max_pedigree_depth($max_pedigree_depth);

   my $alleles = $self->draw_alleles($mafs); # an array ref holding 2*$n_snps alleles.
   my @genotypes = ();
   for my $i_snp (0..$n_snps-1) {
      my $gt =  $alleles->[2*$i_snp] . $alleles->[2*$i_snp + 1]; # [$alleles->[2*$i_snp], $alleles->[2*$i_snp + 1]];
      push @genotypes, $gt;
   }
   #   print "in new_from_pop. [", ref \@genotypes, "] \n";
   $self->set_genotypes(\@genotypes); # e.g. ['AA','Aa','AA', 'aa']
   #   print "in new_from_pop. [", ref $self->get_genotypes(), "] \n";
   $$id++; # increment the id number.
   return $self;
}

sub new_offspring{
   my $class = shift;
   my $self = bless {}, $class;
   my $the_rng = shift;
   my $pgobj1 = shift;          # parent genotype obj. 1
   my $pgobj2 = shift;          # parent genotype obj. 2
   my $generation = shift;
   my $id = shift; # ref to scalar
   my $max_pedigree_depth = shift // 2; # 2 -> keep just back to grandparents, etc.
   my @genotypes = ();
   my $pg1 = $pgobj1->get_genotypes();
   my $pg2 = $pgobj2->get_genotypes();
   die "number of snps is different in the 2 parents. bye. \n" if(scalar @$pg1 != scalar @$pg2);
   # print STDERR "n snps: ", scalar @$pg1, "\n";
   $self->set_rng($the_rng);
   $self->set_generation($generation);
   $self->set_id($$id);
   $self->set_parents($pgobj1->get_id(), $pgobj2->get_id());
   $self->set_max_pedigree_depth($max_pedigree_depth);

   my $pedigree = '(' . $pgobj1->get_pedigree() . ',' . $pgobj2->get_pedigree() . ')' . $$id;
   $self->reduce_and_set_pedigree($pedigree);
 
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
$$id++; # increment the id number;
   return $self;
}

sub new_from_012string{
   my $class = shift;
   my $self = bless {}, $class;
   my $the_rng = shift;
   my $string = shift;          # string of 0's, 1's, and 2's
   my $generation = shift;
   my $id = shift;
   my $pedigree = shift // "()$$id";
   my $max_pedigree_depth = shift // 2;
   $self->set_rng($the_rng);
   $self->set_generation($generation);
   $self->set_id($$id);
   $self->set_parents('X', 'X');

   $self->set_pedigree($pedigree);
   $self->set_max_pedigree_depth($max_pedigree_depth);
 #  print STDERR "012string: $string \n";
   $string =~ s/\s+//g;
   my @zots = split('', $string);
   my @genotypes = ();
   for my $azot (@zots) {
  #    print STDERR "[$azot]\n";
      if ($azot == 0) {
         push @genotypes, 'AA';
      } elsif ($azot == 1) {
         push @genotypes, (gsl_rng_uniform($the_rng->raw()) < 0.5)? 'Aa' : 'aA';
      } elsif ($azot == 2) {
         push @genotypes, 'aa'
      } else {
         die "Unexpected char in 012string: $azot \n";
      }
   }
 #  print STDERR '[', join(", ", @genotypes), "]\n";
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

sub genotype_aAstring{
   my $self = shift;
   my $separator = shift // '';
   my $gt = $self->get_genotypes();
   # print "ref gt: ", ref $gt, "  ", scalar @$gt, "\n";
   my $gstring = '';
   while (my($i, $snp) = each @$gt) {
      # my $snp = $s; # join('', @$s);
      $snp = 'aA' if($snp eq 'Aa');
      $gstring .= $snp . $separator;
   }
   return $gstring;
}

sub genotype_012string{         # AA->0, aA, Aa -> 1, aa->2
   my $self = shift;
   my $separator = shift // '';
   my $gt = $self->get_genotypes();
   # print "ref gt: ", ref $gt, "  ", scalar @$gt, "\n";
   my $gstring = '';
   while (my($i, $snp) = each @$gt) {
      # my $snp = $s; # join('', @$s);
      if ($snp eq 'AA') {
         $gstring .= '0' . $separator;
      } elsif ($snp eq 'aa') {
         $gstring .= '2' . $separator;
      } else {                  # aA or Aa
         $gstring .= '1' . $separator;
      }
   }
   return $gstring;
}

sub fasta{ # 
   my $self = shift;
   my $character_set = shift // '012';
   my $separator = shift // '';

      my ($gen, $id, $pedigree) = ($self->get_generation(), 
                                   $self->get_id(),
                                   $self->get_pedigree()
                                  );
      my $gstring =  ($character_set eq '012')?
        $self->genotype_012string($separator) :
          $self->genotype_aAstring($separator);

   return ("$id   $gen   $pedigree", $gstring);
}

sub reduce_and_set_pedigree{ # if depth of pedigree is greater than specified max, reduce it
   # and set the pedigree of the SimulatedGenotype object to this new, reduced value.
  my $self = shift;
my $pedigree = shift // $self->get_pedigree();
my $max_pedigree_depth = $self->get_max_pedigree_depth();
  $pedigree =~ /^ ( [(]+ )/x; # get initial left parens
 #  print STDERR "pedigree: $pedigree \n";
   my $pedigree_depth = (length $1) - 1; # depth == 2 <=> we have info on grandparents, but not great-grandparent.
 #  print STDERR "  $1   depth $pedigree_depth \n";

   while($pedigree_depth > $max_pedigree_depth) {
      $pedigree =~ s/\(\)\d+,\(\)\d+//g;
      $pedigree =~ /^ ( [(]+ )/x; # get initial left parens
      $pedigree_depth = (length $1) - 1;
   }
   $self->set_pedigree($pedigree);
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

sub set_pedigree{
   my $self = shift;
   my $pedigree = shift;
   $self->{pedigree} = $pedigree;
}

sub get_pedigree{
   my $self = shift;
   return $self->{pedigree};
}

sub set_max_pedigree_depth{
   my $self = shift;
   my $max_pedigree_depth = shift;
   $self->{max_pedigree_depth} = $max_pedigree_depth;
}

sub get_max_pedigree_depth{
   my $self = shift;
   return $self->{max_pedigree_depth};
}


1;

