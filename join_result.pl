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
             "exon_exact_accuracy.txt",
             "intron_exact_accuracy.txt",
             "intron_exact_accuracy.txt",
             "nucleotide_exon_accuracy.txt",
             "nucleotide_exon_partial_accuracy.txt",
             "nucleotide_intron_accuracy.txt",
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

    $total{$file}{PPV} = 0;
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
            if($r =~ m/(.+)\t(\d+)\n\tTP\t(\d+)\n\tFP\t(\d+)\n\tFN\t(\d+)\n\tPPV\t(\d+\.\d+)\n\tSensitivity\t(\d+\.\d+)\n\tF\t(\d+\.\d+).*/){
              push @{$results{$file}{$1}{total}}, $2;
              push @{$results{$file}{$1}{TP}}, $3;
              push @{$results{$file}{$1}{FP}}, $4;
              push @{$results{$file}{$1}{FN}}, $5;
              push @{$results{$file}{$1}{PPV}}, $6;
              push @{$results{$file}{$1}{Sensitivity}}, $7;
              push @{$results{$file}{$1}{F}}, $8;
            }
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
        my $min_tp = min($results{$file}{$entry}{TP});
        my $max_tp = max($results{$file}{$entry}{TP});

        my $fp = mean($results{$file}{$entry}{FP});
        my $var_fp = var($results{$file}{$entry}{FP});
        my $min_fp = min($results{$file}{$entry}{FP});
        my $max_fp = max($results{$file}{$entry}{FP});

        my $fn = mean($results{$file}{$entry}{FN});
        my $var_fn = var($results{$file}{$entry}{FN});
        my $min_fn = min($results{$file}{$entry}{FN});
        my $max_fn = max($results{$file}{$entry}{FN});

        my $sp = mean ($results{$file}{$entry}{PPV});
        my $var_sp = var($results{$file}{$entry}{PPV});
        my $min_sp = min($results{$file}{$entry}{PPV});
        my $max_sp = max($results{$file}{$entry}{PPV});

        my $sn = mean($results{$file}{$entry}{Sensitivity});
        my $var_sn = var($results{$file}{$entry}{Sensitivity});
        my $min_sn = min($results{$file}{$entry}{Sensitivity});
        my $max_sn = max($results{$file}{$entry}{Sensitivity});

        my $total = mean ($results{$file}{$entry}{total});
        my $var_total = var ($results{$file}{$entry}{total});
        my $min_total = min ($results{$file}{$entry}{total});
        my $max_total = max ($results{$file}{$entry}{total});


        my $F = mean ($results{$file}{$entry}{F});
        my $var_F = var ($results{$file}{$entry}{F});
        my $min_F = min ($results{$file}{$entry}{F});
        my $max_F = max ($results{$file}{$entry}{F});

        print OUT "$entry\n";
#        printf OUT ("\t  \tmean:%.2f\tvar:%.2f\tmin:%.2f\tmax:%.2f\n", $total, $var_total, $min_total, $max_total);
 #       printf OUT ("\tTP\tmean:%.2f\tvar:%.2f\tmin:%.2f\tmax:%.2f\n", $tp,$var_tp, $min_tp, $max_tp);
  #      printf OUT ("\tFP\tmean:%.2f\tvar:%.2f\tmin:%.2f\tmax:%.2f\n", $fp, $var_fp, $min_fp, $max_fp);
   #     printf OUT ("\tFN\tmean:%.2f\tvar:%.2f\tmin:%.2f\tmax:%.2f\n", $fn, $var_fn, $min_fn, $max_fn);
        printf OUT ("\tSP\tmean:%.2f\tvar:%.2f\tmin:%.2f\tmax:%.2f\n", $sp, $var_sp, $min_sp, $max_sp);
        printf OUT ("\tSN\tmean:%.2f\tvar:%.2f\tmin:%.2f\tmax:%.2f\n", $sn, $var_sn, $min_sn, $max_sn);
        printf OUT ("\tF \tmean:%.2f\tvar:%.2f\tmin:%.2f\tmax:%.2f\n", $F, $var_F, $min_F, $max_F);
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

sub min {
  my $ref = shift;
  my @values = @{$ref};
  my $total = 0;
  my $min = 1e1000;
  foreach my $v (@values){
    if($min > $v){
      $min = $v;
    }
  }
  return $min;
}

sub max {
  my $ref = shift;
  my @values = @{$ref};
  my $total = 0;
  my $max = -1e1000;
  foreach my $v (@values){
    if($max < $v){
      $max = $v;
    }
  }
  return $max;
}
