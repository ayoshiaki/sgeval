#!/usr/bin/perl

use strict ;
use warnings;

use Data::Dumper;
use Getopt::Long;

my @directory;
my $output;

GetOptions ("d=s{,}" => \@directory,
           "o=s" => \$output);

if( ($#directory <0) || !defined $output) {
  print STDERR "USAGE: $0 -d <dir1> <dir2> ...\n";
  exit();
}

my @files = ("gene_exact_accuracy.txt",
             "gene_overlaped_accuracy.txt",
             "exon_exact_accuracy.txt",
             "intron_exact_accuracy.txt",
             "exon_overlaped_accuracy.txt",
             "intron_exact_accuracy.txt",
             "nucleotide_exon_accuracy.txt",
             "nucleotide_intron_accuracy.txt",
             "start_accuracy.txt",
             "stop_accuracy.txt",
             "acceptor_accuracy.txt",
             "donor_accuracy.txt");
my %results;
my %total;
print "dir\tfile\tpart\ttotal\ttp\tfp\tfn\tppv\tsn\tf\n";

foreach my $file (@files)
  {
    $total{$file}{FP} = 0;
    $total{$file}{FN} = 0;
    $total{$file}{TP} = 0;
    
    $total{$file}{PPV} = 0;
    $total{$file}{Sensitivity} = 0;
    $total{$file}{F} = 0;
    
    
    
    foreach my $dir (@directory)
      {
	for (my $k  = 0; $k < 5; $k++) {
	  open (IN, "<$dir/compara$k/$file") or die "$!";
	  my $txt;
	  foreach my $line (<IN>)
	    {
	      $txt .= $line;
	    }
	  close(IN);
	  
	  my @entries = split("//", $txt);
	  foreach my $r (@entries)
	    {
	      if ($r =~ /^\s*$/)
		{
		  last;
		}
	      if($r =~ m/(.+)\t(\d+)\n\tTP\t(\d+)\n\tFP\t(\d+)\n\tFN\t(\d+)\n\tPPV\t(\d+\.\d+)\n\tSensitivity\t(\d+\.\d+)\n\tF\t(\d+\.\d+).*/){
		print "$dir\t$file\t$k\t$2\t$3\$4\t$5\t$6\t$7\t$8\n";
	      }
	    }
	}
      }
  }
