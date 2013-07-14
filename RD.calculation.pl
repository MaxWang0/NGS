#!/usr/bin/perl -w                                                                                                       
# programmer: Yu Wang                                                                                                    
# usage:                                                                                                                 
# input:                                                                                                                 
# output:                                                                                                                
use strict;                                                                                                              
use warnings;                                                                                                            
use Scalar::Util qw(looks_like_number);                                                                                  
my @region;                                                                                                              
my @read;                                                                                                                
my $nReads = 0;                                                                                                          
my $nRegions = 0;                                                                                                        
my $n = 0;                                                                                                               
my $genename = "";                                                                                                       
my $read_depth = 0;                                                                                                      
my $i=0;                                                                                                                 
my $j = 0;                                                                                                               
my $chrnumber = "";                                                                                                      
open(IN,$ARGV[0]) or die("Fail to open IN.\n"); # read geneprobe stdin                                                   
while (<IN>)                                                                                                             
{                                                                                                                        
    chomp;                                                                                                               
    my @box = split/\t/;                                                                                                 
    $region[$nRegions]->{'chrnumber'}=$box[0];                                                                           
    $region[$nRegions]->{'genename'}=$box[1];                                                                            
    $region[$nRegions]->{'chr'}=$box[2];                                                                                 
    $region[$nRegions]->{'start'}=$box[3];                                                                               
    $region[$nRegions]->{'stop'}=$box[4];                                                                                
    $region[$nRegions]->{'average'}=($box[3]+$box[4])/2;                                                                 
    $nRegions++;                                                                                                         
}                                                                                                                        
close IN;                                                                                                                
open (OUT,">$ARGV[2]") or die("Fail to open OUT.\n"); #read stdout                                                       
open (IN2,$ARGV[1]) or die("Fail to open IN2.\n"); #read rpkm stdin                                                      
while(<IN2>)                                                                                                             
{                                                                                                                        
    if ($nReads>0)                                                                                                       
    {                                                                                                                    
    chomp;                                                                                                               
    my @box2 =split;                                                                                                     
    ($chrnumber = $box2[1]) =~ s/chr//;                                                                                  
    $read[$nReads]->{'chrnumber'} = $chrnumber;                                                                          
    $read[$nReads]->{'start'} = $box2[2];                                                                                
    $read[$nReads]->{'stop'} = $box2[3];                                                                                 
    $read[$nReads]->{'mappedreads'} = $box2[4];                                                                          
    }                                                                                                                    
    $nReads++;                                                                                                           
}                                                                                                                        
                                                                                                                         
print "nRegions:$nRegions\tnReads:$nReads\n";                                                                            
while($i < $nReads and $j < $nRegions)                                                                                   
{                                                                                                                        
   if (looks_like_number($read[$i]{'chrnumber'}) and looks_like_number($region[$j]{'chrnumber'}) and ($read[$i]{'chrnumber'} < $region[$j]{'chrnumber'}))                                                                                       
     {                                                                                                                    
        $i++;                                                                                                            
    }elsif(looks_like_number($read[$i]{'chrnumber'}) and ($region[$j]{'chrnumber'} eq "X" or $region[$j]{'chrnumber'} eq "Y"))                                                                                                                   
    {                                                                                                                    
        $i++;                                                                                                            
    }elsif($read[$i]{'chrnumber'} eq "X" and $region[$j]{'chrnumber'} eq "Y")                                            
    {                                                                                                                    
        $i++;                                                                                                            
    }elsif( looks_like_number($read[$i]{'chrnumber'}) and looks_like_number($region[$j]{'chrnumber'}) and ($read[$i]{'chr    number'} > $region[$j]{'chrnumber'}))                                                                                    
    {                                                                                                                    
        $j++;                                                                                                            
    }elsif(looks_like_number($region[$j]{'chrnumber'}) and ($read[$i]{'chrnumber'} eq "X" or $read[$i]{'chrnumber'} eq "Y    "))                                                                                                                      
    {                                                                                                                    
        $j++;                                                                                                            
    }elsif($read[$i]{'chrnumber'} eq "Y" and $region[$j]{'chrnumber'} eq "X")                                            
    {                                                                                                                    
        $j++;                                                                                                            
    }elsif($read[$i]{'chrnumber'}eq $region[$j]{'chrnumber'}){                                                           
        if($read[$i]{'average'}<$region[$j]{'start'}){                                                                   
            $i++;                                                                                                        
        }elsif($read[$i]{'average'}>=$region[$j]{'stop'}){                                                               
            $j++;                                                                                                        
        }elsif($read[$i]{'average'}>=$region[$j]{'start'} and $read[$i]{'average'}<$region[$j]{'stop'}){                 
            $read_depth = $read_depth + $read[$i]{'mappedreads'};                                                        
            if (($read[$i+1]{'average'}>=$region[$j]{'stop'}) or ($read[$i+1]{'chrnumber'} ne $region[$j]{'chrnumber'})or     ($i == $#read)){                                                                                                        
                print OUT "$region[$j]{'genename'}\t$region[$j]{'chr'}\t$region[$j]{'start'}\t$region[$j]{'stop'}\t$read_    depth\n";                                                                                                                
                $read_depth = 0;                                                                                         
           }                                                                                                            
            $i++;                                                                                                        
        }                                                                                                                
    }                                                                                                                    
}                                                                                                                        
close IN2;                                                                                                               
close OUT; 