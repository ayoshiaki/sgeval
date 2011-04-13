#!/usr/bin/perl

use strict ;
use warnings;

use Data::Dumper;
use Getopt::Long;

my @directory;

GetOptions ("d=s{,}" => \@directory);

if( ($#directory <0)) {
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
print "pred\tdir\tfile\tpart\ttotal\ttp\tfp\tfn\tppv\tsn\tf\n";

foreach my $dir (@directory)
  {
    foreach my $file (@files)
      {
        $dir =~ /cross_validation_(\d+)/;
        my $training_set_size = $1*100;
        my $type = $file;
        $type =~ s/_accuracy.txt//g;
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
		print "$1\t$training_set_size\t$type\t$k\t$2\t$3\t$4\t$5\t$6\t$7\t$8\n";
	      }
	    }
	}
      }
  }
