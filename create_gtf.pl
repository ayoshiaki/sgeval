#!/usr/bin/perl 

use strict ;
use warnings;

use GTF;
use Data::Dumper;
use Getopt::Long;

my @gtf_files;
my $venn;
my $output_dir;

GetOptions ("gtf=s{,}" => \@gtf_files, 
	    "venn=s" => \$venn,
	    "out=s" => \$output_dir);

if($#gtf_files < 0 || !defined ($output_dir)) {
  print STDERR "USAGE: $0 -o <output_directory> -v <venn.txt> -g <reference.gtf> <prediction1.gtf> <prediction2.gtf> ...\n";
  exit();
}

my %gtf_entry_by_geneid;

foreach my $gtf_file (@gtf_files) 
  {
    my $gtf = GTF::new({gtf_filename => $gtf_file});
    my $source = $gtf_file;
    $source =~ s/\.gtf//g;
    $source =~ s%.+/(.+)$%$1%g;
    foreach my $gene (@{$gtf->genes()}) 
      {
	$gtf_entry_by_geneid{$source}{$gene->id()} = $gene;
      }
  }

mkdir $output_dir;
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



  my ($subset, $lixo) = split(/\t/,$lines[0]);
  $subset =~ s#\|#.#g;
  chomp($subset);
  mkdir "$output_dir/$subset";

  my @sources = split(/\./, $subset);
  my %fhs;
  foreach my $source (@sources) {
    my $fh;
    open ($fh, ">$output_dir/$subset/$source.gtf")  or die "$!: $output_dir/$subset/$source.gtf";
    $fhs{$source} = $fh;

  }
  for(my $i = 1; $i < scalar(@lines); $i++) 
    {
      my $l = $lines[$i];
      $l =~ s/\t//g;
      my @transcripts = split(/;/, $l) ;
      foreach my $transcript (@transcripts) {
	my ($source, $geneid) = split(/:/, $transcript);
	$geneid =~ /(.+)?,/;
	$geneid=$1;
	$gtf_entry_by_geneid{$source}{$geneid}->output_gtf($fhs{$source});
      }
    }
  foreach my $source (@sources) {
    close($fhs{$source} );
  }
  
}
