#!/usr/bin/perl 

use strict ;
use warnings;

use Data::Dumper;
use Getopt::Long;

my $venn;

GetOptions ("venn=s" => \$venn);

if( ! defined($venn)) {
  print STDERR "USAGE: $0 -v <venn.txt>\n";
  exit();
}


open (INPUT, "<$venn") or die "$!";
my $venn_str = "";
foreach my $line (<INPUT>) 
  {
    $venn_str .= $line;
  }
close(INPUT);

my @venn = split("//", $venn_str);
  print "set\tlength\tstrand\ttype\n";

foreach my $v (@venn) {
  my @lines = split(/\n/, $v) ;
  if(scalar(@lines) <= 0){
    next;
  }
  while($lines[0] =~ /^\s*$/)
    {
      shift(@lines);
    }
  if(scalar(@lines) <= 0){
    next;
  }
  my ($subset, $lixo) = split(/\t/,$lines[0]);
  $subset =~ s#\|#.#g;
  chomp($subset);
  for(my $i = 1; $i < scalar(@lines); $i++) 
    {
      my $l = $lines[$i];
      $l =~ s/\t//g;
      $l =~ /(.+)?:(\d+)-(\d+),(\+|\-),(\w+)/;
      print $subset."\t".($3 - $2 + 1)."\t".$4."\t".$5."\n";
    }
}
