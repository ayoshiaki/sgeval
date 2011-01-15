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
	     "exon_overlaped_accuracy.txt",
	     "nucleotide_accuracy.txt",
	     "start_accuracy.txt",
	     "stop_accuracy.txt",
	     "acceptor_accuracy.txt",
	     "donor_accuracy.txt");
my %results;
foreach my $file (@files) 
  {
    foreach my $dir (@directory) 
      {
	open (IN, "<$dir/$file") or die "$!";
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
	    $r =~ m/(.+)\t(\d+)\n\tTP\t(\d+)\n\tFP\t(\d+)\n\tFN\t(\d+)\n.*/;
	    $results{$file}{$1}{total} += $2;
	    $results{$file}{$1}{TP} += $3;
	    $results{$file}{$1}{FP} += $4;
	    $results{$file}{$1}{FN} += $5;
	  }
      }
  }
mkdir $output;
foreach my $file (keys %results) 
  {
    open (OUT, ">$output/$file") or die "$!";
    foreach my $entry (keys %{$results{$file}}) 
      {
	my $tp = $results{$file}{$entry}{TP};
	my $fp = $results{$file}{$entry}{FP};
	my $fn = $results{$file}{$entry}{FN};
	print OUT "$entry\t$results{$file}{$entry}{total}\n";
	print OUT "\tTP\t$tp\n";
	print OUT "\tFP\t$fp\n";
	print OUT "\tFN\t$fn\n";
	printf OUT ("\tSpecificity\t%.2f\n\tSensitivity\t%.2f\n", (100.0*($tp/($tp + $fp))),(100.0*($tp/($tp + $fn))));
	print OUT "//\n";
      }
    close(OUT);
  }
