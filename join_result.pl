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
             "nucleotide_accuracy.txt",
             "start_accuracy.txt",
             "stop_accuracy.txt",
             "acceptor_accuracy.txt",
             "donor_accuracy.txt");
my %results;
my %total;
foreach my $file (@files)
  {
    $total{$file}{FP} = 0;
    $total{$file}{FN} = 0;
    $total{$file}{TP} = 0;

    $total{$file}{Specificity} = 0;
    $total{$file}{Sensitivity} = 0;
    $total{$file}{F} = 0;



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
            $r =~ m/(.+)\t(\d+)\n\tTP\t(\d+)\n\tFP\t(\d+)\n\tFN\t(\d+)\n\tSpecificity\t(\d+)\n\tSensitivity\t(\d+)\n\tF\t(\d+).*/;
            push @{$results{$file}{$1}{total}}, $2;
            push @{$results{$file}{$1}{TP}}, $3;
            push @{$results{$file}{$1}{FP}}, $4;
            push @{$results{$file}{$1}{FN}}, $5;
            push @{$results{$file}{$1}{Specificity}}, $6;
            push @{$results{$file}{$1}{Sensitivity}}, $7;
            push @{$results{$file}{$1}{F}}, $8;
          }
      }
  }

mkdir $output;
foreach my $file (keys %results)
  {
    open (OUT, ">$output/$file") or die "$!";
    foreach my $entry (keys %{$results{$file}})
      {

        my $tp = mean($results{$file}{$entry}{TP});
        my $var_tp = var($results{$file}{$entry}{TP});

        my $fp = mean($results{$file}{$entry}{FP});
        my $var_fp = var($results{$file}{$entry}{FP});

        my $fn = mean($results{$file}{$entry}{FN});
        my $var_fn = var($results{$file}{$entry}{FN});

        my $sp = mean ($results{$file}{$entry}{Specificity});
        my $var_sp = var($results{$file}{$entry}{Specificity});

        my $sn = mean($results{$file}{$entry}{Sensitivity});
        my $var_sn = var($results{$file}{$entry}{Sensitivity});

        my $total = mean ($results{$file}{$entry}{total});
        my $var_total = var (@{$results{$file}{$entry}{total}});

        my $F = mean ($results{$file}{$entry}{F});
        my $var_F = var ($results{$file}{$entry}{F});

        print OUT "$entry\t$total\t$var_total\n";
        print OUT "\tTP\t$tp\t$var_tp\n";
        print OUT "\tFP\t$fp\t$var_fp\n";
        print OUT "\tFN\t$fn\t$var_fn\n";
        printf OUT ("\tSpecificity\t%.2f\t%.2f\n\tSensitivity\t%.2f\t%.2f\n\tF\t%.2f\t%.2f", $sp,$var_sp, $sn, $var_sn, $F, $var_F);
        print OUT "//\n";
      }
    close(OUT);
  }

sub var {
  my $ref = shift;
  my @values = @{$ref};
  my $n    = 0;
  my $sum1 = 0;
  my $sum2 = 0;

  foreach my $x (@values) {
    $n    = $n + 1;
    $sum1 = $sum1 + $x;
  }
  my $mean = $sum1/$n;

  foreach my $x (@values) {
    $sum2 = $sum2 + ($x - $mean)*($x - $mean);
  }
  my $variance = $sum2/($n - 1);
  return $variance;
}

sub mean {
  my $ref = shift;
  my @values = @{$ref};
  my $total = 0;
  foreach my $v (@values){
    $total += $v;
  }
  return $total/($#values + 1);
}
