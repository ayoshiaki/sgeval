#!/usr/bin/perl -w 

use strict;
use warnings;

use GTF;
use Data::Dumper;
use Getopt::Long;

my @gtf_files;
my $output_dir;

GetOptions ("gtf=s{,}" => \@gtf_files, 
	    "out=s" => \$output_dir);

if($#gtf_files < 0 || !defined ($output_dir)) {
  print STDERR "USAGE: $0 -o <output_directory> -g <reference.gtf> <prediction1.gtf> <prediction2.gtf> ...\n";
  exit();
}


# reading all GTFs
my %sites;
my %component;
my $ref_source;
my $first_source = 1;
foreach my $gtf_file (@gtf_files) 
  {
    my $gtf = GTF::new({gtf_filename => $gtf_file});
    my $source = $gtf_file;
    $source =~ s/\.gtf//g;
    $source =~ s%.+/(.+)$%$1%g;
    if($first_source == 1) {
      $ref_source = $source;
      $first_source = 0;
    }
    foreach my $gene (@{$gtf->genes()}) 
      {
	
	if($gene->strand() eq "+") 
	  {
	    process_forward($gene, \%sites, $source);
	  } else {
	    process_reverse($gene, \%sites, $source);
	  }
      }
  }
build_components(\%component, \%sites);


my %component_by_seqname;
foreach my $entry  (sort {$a <=> $b} keys %component)
  {
    ${$component{$entry}}[0] =~ m/(.+)?:.*/;
    push @{$component_by_seqname{$1}},$entry;	
  }



gene_venn();



sub gene_venn { 
  my %subsets;
  foreach my $seqname (keys %component_by_seqname) 
    {
      foreach my $component (@{$component_by_seqname{$seqname}}) 
	{
	  my %recticulate;
	  foreach my $node (@{$component{$component}})
	    {
	      my $label = "";
	      foreach my $source (sort {$a cmp $b} (@{$sites{$node}->{Source}})) 
		{
		  $label .= "<$source>";
		}
	      $recticulate{$label} = {};
	    }
	  foreach my $from (keys %recticulate) {
	    foreach my $to (keys %recticulate) {
	      if($from eq $to){
		next;
	      }
	      if($to =~ /$from/) 
		{
		  push @{$recticulate{$from}->{Next}}, $to;
		  push @{$recticulate{$to}->{From}}, $from;
		}
	    }
	  }
	  print $seqname."\n";
	  foreach my $from (keys %recticulate) {
	    if(!defined $recticulate{$from}->{From}) {
	      $from =~ s/></;/g;
	      $from =~ s/<//g;
	      $from =~ s/>//g;
	      
	      print " ".$from."\n";
	    }
	  }
	}
    }
}





sub process_forward {
    my $gene = shift;
    my $sites_r = shift;
    my $gtf_filename = shift;
    my $geneid = $gene->id();
    foreach my $tx (@{$gene->transcripts()}) 
    {
	my @start_codons = @{$tx->start_codons()};
	my @stop_codons = @{$tx->stop_codons()};
	my $last_right_site;
	foreach my $cds (@{$tx->cds()} ) 
	{
	    my $left_site;
	    my $right_site;
	    my $is_start_codon = 0;
	    foreach my $start_codon (@start_codons) 
	    {
		if($cds->start() ==  $start_codon->start()) 
		{
		    $is_start_codon = 1;
		    last;
		}
	    }
	    my $is_stop_codon = 0;
	    foreach my $stop_codon (@stop_codons) 
	    {
		if(($cds->stop() + 1) ==  $stop_codon->start()) 
		{
		    $is_stop_codon = 1;
		    last;
		}
	    }
	    
	    if($is_start_codon) {
		$left_site = 
		{
		    SeqName => $gene->seqname(),
		    Position => $cds->start(),
		    Type=> "start_codon",
		    Frame => $cds->frame(),
		    Strand => "+"
		};
	    } else {
		$left_site = 
		{
		    SeqName =>$gene->seqname(), 
		    Position => $cds->start(),
		    Type => "acceptor",
		    Frame => $cds->frame(),
		    Strand => "+"
		}
	    }
	    if($is_stop_codon) {
		$right_site = 
		{
		    SeqName =>$gene->seqname(), 
		    Position => $cds->stop(),
		    Type => "stop_codon",
		    Frame => $cds->frame(),
		    Strand => "+"
		};
	    } else {
		$right_site = 
		{
		    SeqName =>$gene->seqname(), 
		    Position => $cds->stop(),
		    Type => "donor", 
		    Frame => $cds->frame(),
		    Strand => "+"
		};
	    }
	    my $source = $gtf_filename;
	    $left_site = add_site($left_site, $sites_r, $source, $gene->id());
	    $right_site = add_site($right_site, $sites_r, $source, $gene->id());
	    if(defined $last_right_site) 
	    {
		push @{$last_right_site->{Next}->{get_key_from_site($left_site)}}, $source.":". $gene->id();
		push @{$left_site->{From}->{get_key_from_site($last_right_site)}}, $source.":". $gene->id();
	    }
	    push @{$left_site->{Next}->{get_key_from_site($right_site)}}, $source.":". $gene->id();
	    push @{$right_site->{From}->{get_key_from_site($left_site)}}, $source.":". $gene->id();
	    $last_right_site = $right_site;
	}
    }
}


sub process_reverse {
    my $gene = shift;
    my $sites_r = shift;
    my $gtf_filename = shift;
    foreach my $tx (@{$gene->transcripts()}) 
    {
	my @start_codons = @{$tx->start_codons()};
	my @stop_codons = @{$tx->stop_codons()};
	my $last_right_site;
	foreach my $cds (@{$tx->cds()} ) 
	{
	    my $left_site;
	    my $right_site;
	    my $is_start_codon = 0;
	    foreach my $start_codon (@start_codons) 
	    {
		if($cds->stop() ==  $start_codon->stop()) 
		{
		    $is_start_codon = 1;
		    last;
		}
	    }
	    my $is_stop_codon = 0;
	    foreach my $stop_codon (@stop_codons) 
	    {
		if(($cds->start() - 1) ==  $stop_codon->stop()) 
		{
		    $is_stop_codon = 1;
		    last;
		}
	    }
	    
	    if($is_start_codon) {
		$right_site = 
		{
		    SeqName => $gene->seqname(),
		    Position => $cds->stop(),
		    Type=> "start_codon",
		    Frame => $cds->frame(),
		    Strand => "-"
		};
	    } else {
		$right_site = 
		{
		    SeqName =>$gene->seqname(), 
		    Position => $cds->stop(),
		    Type => "acceptor",
		    Frame => $cds->frame(),
		    Strand => "-"
		}
	    }
	    if($is_stop_codon) {
		$left_site = 
		{
		    SeqName =>$gene->seqname(), 
		    Position => $cds->start(),
		    Type => "stop_codon",
		    Frame => $cds->frame(),
		    Strand => "-"
		};
	    } else {
		$left_site = 
		{
		    SeqName =>$gene->seqname(), 
		    Position => $cds->start(),
		    Type => "donor", 
		    Frame => $cds->frame(),
		    Strand => "-"
		};
	    }
	    my $source = $gtf_filename;
	    $left_site = add_site($left_site, $sites_r, $source, $gene->id());
	    $right_site = add_site($right_site, $sites_r, $source, $gene->id());
	    if(defined $last_right_site) 
	    {
		push @{$last_right_site->{Next}->{get_key_from_site($left_site)}}, $source.":".$gene->id();
		push @{$left_site->{From}->{get_key_from_site($last_right_site)}}, $source.":". $gene->id();
	    }
	    push @{$left_site->{Next}->{get_key_from_site($right_site)}}, $source.":". $gene->id();
	    push @{$right_site->{From}->{get_key_from_site($left_site)}}, $source.":".$gene->id();
	    $last_right_site = $right_site;
	}
    }
}

sub add_site {
    my $site = shift;
    my $sites_r = shift;
    my $gene_source = shift;
    my $gene_id = shift;
    my $key = get_key_from_site($site);
    if(! defined $sites_r-> {$key}) 
    {
	$site = $sites_r->{$key} = $site;
    } else {
	$site = $sites_r->{$key};
    }
    push @{$site->{Source}}, "$gene_source:$gene_id";
    return $site;
}

sub get_key_from_site {
    my $site = shift;
    return $site->{SeqName}.":".$site->{Position}.",".$site->{Type};
}



sub build_components {
  my $component_ref = shift;
  my $sites_ref = shift;
  my $id = 1;
  my @nodes = keys (%{$sites_ref});
  my %marked;
  foreach my $n (@nodes) {
    if (defined $marked{$n} ) 
      {
	next;
      }
    
    my @queue; 
    push @queue, ($n);
    push @{$component_ref->{$id}},$n;
    $marked{$n} = 1;
    while (scalar (@queue) > 0) 
      {
	my $node = pop @queue;
	my @edges = (keys %{$sites_ref->{$node}->{Next}});
	if(defined  $sites_ref->{$node}->{From}){
	  foreach my $y (keys %{ $sites_ref->{$node}->{From} })
	    {
	      push @edges, $y;
	    }
	}
	foreach my $x (@edges) {
	  if(!defined ($marked{$x})){
	    $marked{$x} = 1;
	    push @queue, $x;
	    push @{$component_ref->{$id}},$x;
	  }
	}
      }
    $id ++;
  }
}
