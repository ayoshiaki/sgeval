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
  print "set\tsource\tname\tcount\ttranscript_id\n";
  my $id = 0;

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
      $id++;
      foreach my $entry (split(/;/, $l)) 
	{
	  if( $entry =~ /(.+)?:(.+)?,(\d+)/) {
	    print $subset."\t".$1."\t".$2."\t".$3."\t".$id."\n";
	  }
	}


    }
}
