#!/usr/bin/perl -w
# programmer: Yu Wang
# usage:
# input:
# output:
use strict;
use warnings;
my $genename = "";
my @region;
my $nRegions = 0;
my $n =0;
my $start = 14362;
my $stop = 0;
my $exon = 1;
open(IN,$ARGV[0]) or die("Fail to open IN.\n"); #read target stdin
while (<IN>)
{
    chomp;
    my @box = split;
    ($genename = $box[3]) =~ s/.*\:.*\://;
    $region[$nRegions]->{'genename'}=$genename;
    $region[$nRegions]->{'line'}= $_;
    $region[$nRegions]->{'index'}=$box[3];
    $region[$nRegions]->{'chr'}=$box[0];
    $region[$nRegions]->{'start'}=$box[1];
    $region[$nRegions]->{'stop'}=$box[2];
    $nRegions++;
}
close IN;
open(OUT,">$ARGV[1]") or die("Fail to open OUT.\n"); #read target stdout
for(@region){
    if ($n < $nRegions-1){

    if ($region[$n]{'genename'} ne $region[$n+1]{'genename'}){
        $stop = $region[$n]{'stop'};
        print OUT "$region[$n]{'genename'}\t$region[$n]{'chr'}\t$start\t$stop\t$exon\n";
        $start= $region[$n+1]{'start'};
        $exon = 1;
    }else{
        $exon++;
    }
    }elsif($n == $nRegions-1){
         $stop = $region[$n]{'stop'};
        print OUT "$region[$n]{'genename'}\t$region[$n]{'chr'}\t$start\t$stop\t$exon\n";
    }
    $n++;
    }
close OUT;


