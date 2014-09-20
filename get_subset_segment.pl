#!/usr/bin/perl

use strict ;
use warnings;

use Data::Dumper;
use Getopt::Long;

my $venn;
my $subset;
GetOptions ("venn=s" => \$venn, "subset=s" => \$subset);

if( ! defined($venn) || ! defined($subset) ) {
  print STDERR "USAGE: $0 -v <venn.txt> -s <subset name>\n";
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
  my ($subset2, $lixo) = split(/\t/,$lines[0]);
  $subset2 =~ s#\|#.#g;
  chomp($subset2);
  if($subset2 =~ /$subset/) {
    for(my $i = 1; $i < scalar(@lines); $i++)
      {
        my $l = $lines[$i];
        $l =~ s/\t//g;
        $l =~ /(.+)?:(\d+)-(\d+),(\+|\-)(,(\w+))?/;
        print $1." ".($2)." ".($3)."\n";
      }
  }
}
