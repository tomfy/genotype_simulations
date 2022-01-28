#!/usr/bin/perl -w
use strict;



print "MARKER \n";
while(<>){
    die if(! /^>(\S+)/);
    my $id = $1;
    
    my $next_line = <>;
    $next_line =~ /^\s*(\S+)/;
    my $genotypes = $1;
    my $nmarkers = length $genotypes;
    print "$id  $genotypes\n";
}
	
