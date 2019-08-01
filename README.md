Simulating genotypes, and analyzing pairs of genotypes.

Genotype.pm  Genotype class
      Can construct either from 
          * a set of minor allele frequencies ( new_from_mafs ), or
          * two parent Genotype objects ( new_offspring ), or 
          * a string of 0s, 1s, and 2s, 0 = AA, 1 = heterozygou, 2 = aa ( new_from012string ).

ggg.pl  generate generations of genotypes
        e.g.
        ggg.pl  -pop 50 -nsnp 500 -ngen 3 -nout 3 -maf 'delta,0.25' -sample 100

        This will make 3 generations of a population of size 50, with each genotype having 500 snps. A random 100 of the 150 genotypes are output in fasta format. 