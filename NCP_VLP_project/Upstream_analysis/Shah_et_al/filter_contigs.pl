#!/usr/bin/perl
use strict;
use warnings;

# copied from https://www.biostars.org/p/79202/#79467
# Author: Vivek Krishnakumar

my $minlen = shift or die "Error: `minlen` parameter not provided\n";
{
    local $/=">";
    while(<>) {
    	chomp;
    	next unless /\w/;
    	s/>$//gs;
    	my @chunk = split /\n/;
    	my $header = shift @chunk;
    	my $seqlen = length join "", @chunk;
    	print ">$_" if($seqlen >= $minlen);
    }
    local $/="\n";
}

