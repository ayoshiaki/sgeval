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
mkdir $output_dir;

# reading all GTFs
my %sites;
my %component;
my @sources;
my $ref_source;
my $first_source = 1;
my %transcripts;
foreach my $gtf_file (@gtf_files)
  {
    my $gtf = GTF::new({gtf_filename => $gtf_file});
    my $source = $gtf_file;
    $source =~ s/\.gtf//g;
    $source =~ s%.+/(.+)$%$1%g;
    push @sources, $source;

    if($first_source == 1) {
      $ref_source = $source;
      $first_source = 0;
    }
    foreach my $gene (@{$gtf->genes()})
      {
        foreach my $tx (@{$gene->transcripts})
          {
            $transcripts{$tx->id} = $tx;
          }

        if($gene->strand() eq "+")
          {
            process_forward($gene, \%sites, $source);
          } else {
            process_reverse($gene, \%sites, $source);
          }
      }
  }
# build connected components
build_components(\%component, \%sites);

my %component_by_seqname;
foreach my $entry  (sort {$a <=> $b} keys %component)
  {
    ${$component{$entry}}[0] =~ m/(.+)?:.*/;
    push @{$component_by_seqname{$1}},$entry;
  }

my $number_of_transcripts = 0;
my %gvenn_overlaped = gene_overlap_venn();
my %gvenn_exact = gene_exact_venn();
my %exon_exact = exon_exact_venn();
my %intron_exact = intron_exact_venn();
my %start = start_codon_venn();
my %stop = stop_codon_venn();
my %acceptor = acceptor_venn();
my %donor = donor_venn();

my %exon_overlaped = exon_overlaped_venn();
my %nucleotide = nucleotide_venn();

generate_result("gene_overlaped", \%gvenn_overlaped);
generate_result("gene_exact", \%gvenn_exact);
generate_result("exon_exact", \%exon_exact);
generate_result("intron_exact", \%intron_exact);
generate_result("exon_overlaped", \%exon_overlaped);
generate_result("start", \%start);
generate_result("stop", \%stop);
generate_result("acceptor", \%acceptor);
generate_result("donor", \%donor);
generate_result("nucleotide", \%nucleotide);


sub generate_result {
  my $output_filename = shift;
  my $ref_venn = shift;
  my %venn = %{$ref_venn};
  open (OUTPUT, ">$output_dir/$output_filename"."_venn.txt") or die "$!";
  foreach my $key (keys %venn) {
    print OUTPUT $key."\t".scalar @{$venn{$key}}."\n";
    foreach my $el ( @{$venn{$key}} )
      {
        print OUTPUT "\t".$el."\n";
      }
    print OUTPUT "//\n";
  }
  close(OUTPUT);
  open (OUTPUT, ">$output_dir/$output_filename"."_accuracy.txt") or die "$!";
  foreach my $source (@sources) {
    my $tp = 0;
    my $fp = 0;
    my $fn = 0;
    if($source eq $ref_source) {
      next;
    }
    foreach my $subset (keys %venn)
      {
        my $count = scalar(@{$venn{$subset}});
        my $a = ($subset =~ /^$source$/) || ($subset =~ /^$source(\|)/)|| ($subset =~ /(\|)$source(\|)/) || ($subset =~ /(|)$source$/);
        my $b = ($subset =~ /^$ref_source$/) || ($subset =~ /^$ref_source(\|)/)|| ($subset =~ /(\|)$ref_source(\|)/) || ($subset =~ /(|)$ref_source$/);
        if($a && $b) {
          $tp += $count;
        } elsif(!($a) && ( $b)) {
          $fn += $count;
        } elsif(($a) && !( $b)) {
          $fp += $count;
        } else {
        }
      }
    my $sp = (100.0*($tp/($tp + $fp)));
    my $sn = (100.0*($tp/($tp + $fn)));
    my $f = 0;
    if(($sp + $sn) != 0) {
      $f = 2 * $sp * $sn / ($sp + $sn);
    }
    print OUTPUT $source."\t".($tp+$fp)."\n";
    print OUTPUT "\tTP\t$tp\n\tFP\t$fp\n\tFN\t$fn\n";
    printf OUTPUT ("\tSpecificity\t%.2f\n\tSensitivity\t%.2f\n", (100.0*($tp/($tp + $fp))),(100.0*($tp/($tp + $fn))));
    printf OUTPUT ("\tF\t%.2f\n", $f);
    print OUTPUT "//\n";
  }
  close(OUTPUT);

}



sub gene_exact_venn {
  my %subsets;
  my %gvenn;
  foreach my $seqname (keys %component_by_seqname)
    {
      foreach my $component (@{$component_by_seqname{$seqname}})
        {
          my %recticulate;
          foreach my $node (@{$component{$component}})
            {
              my $label = "";
              my $first = 1;
              foreach my $next_node  (keys %{$sites{$node}->{Next}})
                {
                  my @transcripts = @{$sites{$node}->{Next}->{$next_node}};
                  foreach my $source  (sort {$a cmp $b} (@transcripts))
                    {
                      if($first ) {
                        $label .= "$source";
                        $first = 0;
                      }
                      else {
                        $label .= ";$source";
                      }
                    }
                  $recticulate{$label} = {};
                }
            }


          %recticulate = %{build_recticulate(\%recticulate)};



          my @sorted= sort { my @aa = split(";", $a); my @bb = split(";", $b); @aa <=> @bb } keys %recticulate ;
          my $k = 0;
          while (($k < scalar(@sorted)) && (scalar(@sorted) > 0))
            {
              my $from = $sorted[$k];
              if(scalar (@{$recticulate{$from}->{From}}) <= 1) {
                my $subsets = build_subset_string($from);
                my $str = "";
                my $nexon = 0;
                foreach my $source (split(";", $from))
                  {
		      
		      if($source =~ m/(.+)?:(.+)/)
		      {
			  $nexon = count_exon($2);
			  $str .= "$1:$2,$nexon;";
		      } else {
			  print STDERR "Something wrong: $from\n";
		      }
                  }

                push @{$gvenn{$subsets}}, $str;
                %recticulate = () ;
                foreach my $el (@sorted)
                  {
                    my @remove_subset = split(";", $from);
		    

                    foreach my $xx (@remove_subset) {
                      $el =~ s/;$xx$//g;
                      $el =~ s/^$xx;//g;
                      $el =~ s/;$xx;/;/g;
                      $el =~ s/^$xx$//g;
                    }

                    if(!$el =~/^\s*$/){
                      $recticulate{$el} = {};
                    }

                  }
                %recticulate = %{build_recticulate(\%recticulate)};


                @sorted= sort { my @aa = split(";", $a); my @bb = split(";", $b); @aa <=> @bb } keys %recticulate ;
                $k = 0;
              } else {
                $k++;
              }
            }
        }

    }

  return %gvenn;
}






sub build_subset_string {
  my $str = shift;
  my @sets = split(";", $str);
  my $firsttime = 1;
  my $subsets = "";
  my %aux = ();
  foreach my $set (sort {$a cmp $b} (@sets))
    {
      if($set =~/(.+)?:(.+)/) {
        if(defined $aux{$1}) {
          next;
        }
        $aux{$1} = 1;
        if($firsttime) {
          $subsets .= $1;
          $firsttime = 0;
        } else {
          $subsets .= "|".$1;
        }
      }
    }
  return $subsets;
}
sub gene_overlap_venn {
  my %subsets;
  my %gvenn;
  foreach my $seqname (keys %component_by_seqname)
    {
      foreach my $c (@{$component_by_seqname{$seqname}})
        {
          my %aux;
          my $list ="";
          my $first = 1;

          foreach my $x (@{$component{$c}}) {
            foreach my $s (@{ $sites{$x}->{Source}}){
              $aux{$s} = 1;
            }
          }
          foreach my $x (keys %aux) {
            if($x =~ m/(.+)?:(.+)/){
              my $nexon = count_exon($2);
              if($first) {
                $list .= "$x,$nexon";
                $first = 0;
              } else{
                $list .= ";$x,$nexon";
              }
            }

            my $subset = build_subset_string($list);
            push @{$gvenn{$subset}}, $list;
          }
        }

    }
  return %gvenn;
}


sub count_exon {
  my $id = shift;
  if($id) {
    my $tx = $transcripts{$id};
    return scalar(@{$tx->cds()});
  }
    return 0;
}

sub build_recticulate {
  my $ref = shift;
  my %recticulate = %{$ref};
  foreach my $from (keys %recticulate) {
    foreach my $to (keys %recticulate) {
      my @firstset = split(";", $from);
      my @secondset = split(";", $to);
      my $is_subset = 1;
      foreach my $el (@firstset)  {
        my $found = 0;
        foreach my $el2 (@secondset) {
          if($el eq $el2) {
            $found = 1;
            last;
          }
        }
        if (!$found) {
          $is_subset = 0;
        }
      }
      if($is_subset)
        {
          push @{$recticulate{$from}->{Next}}, $to;
          push @{$recticulate{$to}->{From}}, $from;
        }
    }
  }
  return \%recticulate;
}


sub nucleotide_venn {
  my %nucleotide_venn;

  foreach my $seqname (keys %component_by_seqname)
    {
      my %intervals;

      foreach my $c (@{$component_by_seqname{$seqname}})
        {
          foreach my $node (@{$component{$c}})
            {
              if(
                 ($node =~ /start/ && $sites{$node}->{Strand} eq "+")
                 ||
                 ($node =~ /acceptor/ && $sites{$node}->{Strand} eq "+")
                 ||
                 ($node =~ /stop/ && $sites{$node}->{Strand} eq "-")
                 ||
                 ($node =~ /donor/ && $sites{$node}->{Strand} eq "-")
                )
                {
                  my $strand = $sites{$node}->{Strand};
                  foreach my $next_node (keys %{$sites{$node}->{Next}})
                    {
                      $node =~ m/(.+)?:(\d+),(.+)/;
                      my $start = $2;
                      $next_node =~ m/(.+)?:(\d+),(.+)/;
                      my $end = $2;
                      for (my $i = $start;  $i <= $end; $i++)
                        {
                          foreach my $source (@{$sites{$node}->{Next}->{$next_node}}) {
                            push @{$intervals{$strand}{$i}}, $source;
                          }
                        }
                    }
                }
            }
          foreach my $strand (keys %intervals) {
            foreach my $i (sort {$a <=> $b} (keys %{$intervals{$strand}}))
              {
                my %aux;
                foreach my $source (@{$intervals{$strand}{$i}})
                  {
                    $source =~ /(.+)?:(.+)/;
                    $aux{$1} = 1;
                  }
                my $subset;
                my $first = 1;
                foreach my $k (sort {$a cmp $b} (keys %aux)) {
                  if($first) {
                    $first = 0;
                    $subset .= $k;
                  } else {
                    $subset .= "|$k";
                  }
                }

                push @{$nucleotide_venn{$subset}}, "$seqname,$strand:$i";
              }
          }
        }
    }
  return %nucleotide_venn;
}


sub exon_overlaped_venn {
  my %exon_overlaped_venn;
  foreach my $seqname (keys %component_by_seqname)
    {
      my %intervals;

      foreach my $c (@{$component_by_seqname{$seqname}})
        {
          foreach my $node (@{$component{$c}})
            {
              if(
                 ($node =~ /start/ && $sites{$node}->{Strand} eq "+")
                 ||
                 ($node =~ /acceptor/ && $sites{$node}->{Strand} eq "+")
                 ||
                 ($node =~ /stop/ && $sites{$node}->{Strand} eq "-")
                 ||
                 ($node =~ /donor/ && $sites{$node}->{Strand} eq "-")
                )
                {
                  my $strand = $sites{$node}->{Strand};
                  foreach my $next_node (keys %{$sites{$node}->{Next}})
                    {
                      $node =~ m/(.+)?:(\d+),(.+)/;
                      my $start = $2;
                      $next_node =~ m/(.+)?:(\d+),(.+)/;
                      my $end = $2;
                      for (my $i = $start;  $i <= $end; $i++)
                        {
                          foreach my $source (@{$sites{$node}->{Next}->{$next_node}}) {
                            push @{$intervals{$strand}{$i}}, $source;
                          }
                        }
                    }
                }
            }



          foreach my $node (@{$component{$c}})
            {
              if(
                 ($node =~ /start/ && $sites{$node}->{Strand} eq "+")
                 ||
                 ($node =~ /acceptor/ && $sites{$node}->{Strand} eq "+")
                 ||
                 ($node =~ /stop/ && $sites{$node}->{Strand} eq "-")
                 ||
                 ($node =~ /donor/ && $sites{$node}->{Strand} eq "-")
                )
                {
                  my $strand = $sites{$node}->{Strand};
                  foreach my $next_node (keys %{$sites{$node}->{Next}})
                    {
                      $node =~ m/(.+)?:(\d+),(.+)/;
                      my $start = $2;
                      $next_node =~ m/(.+)?:(\d+),(.+)/;
                      my $end = $2;
                      my %aux;
                      for (my $i = $start;  $i <= $end; $i++)
                        {
                          foreach my $source (@ {$intervals{$strand}{$i}}){
                            $source =~ /(.+)?:(.+)/;
                            $aux{$1} = 1;
                          }
                        }
                      my $subset;
                      my $first = 1;
                      foreach my $k (sort {$a cmp $b} (keys %aux)) {
                        if($first) {
                          $first = 0;
                          $subset .= $k;
                        } else {
                          $subset .= "|$k";
                        }
                      }
                      my $type;
                      if($node =~ /start/ && $sites{$node}->{Strand} eq "+")
                        {
                          $type = "initial";
                          if ($next_node =~/stop/){
                            $type = "single";
                          }

                        }
                      elsif ($node =~ /acceptor/ && $sites{$node}->{Strand} eq "+")
                        {
                          $type = "internal";
                          if($next_node =~ /stop/) {
                            $type = "final";
                          }
                        }
                      elsif ($node =~ /stop/ && $sites{$node}->{Strand} eq "-")
                        {
                          $type = "final";
                          if($next_node =~ /start/) {
                            $type = "single";
                          }
                        }
                      elsif ($node =~ /donor/ && $sites{$node}->{Strand} eq "-")
                        {
                          $type = "internal";
                          if($next_node =~ /start/) {
                            $type = "initial";
                          }
                        }

                      push @{$exon_overlaped_venn{$subset}} , "$seqname:".$sites{$node}->{Position}."-".$sites{$next_node}->{Position}.",".$sites{$node}->{Strand}.",".$type;



                    }
                }
            }


        }
    }

  return %exon_overlaped_venn;
}



sub exon_exact_venn {
  my %exon_venn;
  foreach my $seqname (keys %component_by_seqname)
    {
      foreach my $c (@{$component_by_seqname{$seqname}})
        {
          foreach my $node (@{$component{$c}})
            {
              if(
                 ($node =~ /start/ && $sites{$node}->{Strand} eq "+")
                 ||
                 ($node =~ /acceptor/ && $sites{$node}->{Strand} eq "+")
                 ||
                 ($node =~ /stop/ && $sites{$node}->{Strand} eq "-")
                 ||
                 ($node =~ /donor/ && $sites{$node}->{Strand} eq "-")
                )
                {

                  foreach my $next_node (keys %{$sites{$node}->{Next}})
                    {
                      my %aux;
                      foreach my $source (@{$sites{$node}->{Next}->{$next_node}})
                        {
                          $source =~ /(.+)?:(.+)/;
                          $aux{$1} = 1;
                        }
                      my $subset;
                      my $first = 1;
                      foreach my $k (sort {$a cmp $b} (keys %aux)) {
                        if($first) {
                          $first = 0;
                          $subset .= $k;
                        } else {
                          $subset .= "|$k";
                        }
                      }
                      my $type;
                      if($node =~ /start/ && $sites{$node}->{Strand} eq "+")
                        {
                          $type = "initial";
                          if ($next_node =~/stop/){
                            $type = "single";
                          }

                        }
                      elsif ($node =~ /acceptor/ && $sites{$node}->{Strand} eq "+")
                        {
                          $type = "internal";
                          if($next_node =~ /stop/) {
                            $type = "final";
                          }
                        }
                      elsif ($node =~ /stop/ && $sites{$node}->{Strand} eq "-")
                        {
                          $type = "final";
                          if($next_node =~ /start/) {
                            $type = "single";
                          }
                        }
                      elsif ($node =~ /donor/ && $sites{$node}->{Strand} eq "-")
                        {
                          $type = "internal";
                          if($next_node =~ /start/) {
                            $type = "initial";
                          }
                        }

                      push @{$exon_venn{$subset}} , "$seqname:".$sites{$node}->{Position}."-".$sites{$next_node}->{Position}.",".$sites{$node}->{Strand}.",".$type;

                    }
                }
            }
        }
    }
  return %exon_venn;
}





sub intron_exact_venn {
  my %intron_venn;
  foreach my $seqname (keys %component_by_seqname)
    {
      foreach my $c (@{$component_by_seqname{$seqname}})
        {
          foreach my $node (@{$component{$c}})
            {
              if(
                 ($node =~ /donor/ && $sites{$node}->{Strand} eq "+")
                 ||
                 ($node =~ /acceptor/ && $sites{$node}->{Strand} eq "-")
                )
                {
                  foreach my $next_node (keys %{$sites{$node}->{Next}})
                    {
                      my %aux;
                      foreach my $source (@{$sites{$node}->{Next}->{$next_node}})
                        {
                          $source =~ /(.+)?:(.+)/;
                          $aux{$1} = 1;
                        }
                      my $subset;
                      my $first = 1;
                      foreach my $k (sort {$a cmp $b} (keys %aux)) {
                        if($first) {
                          $first = 0;
                          $subset .= $k;
                        } else {
                          $subset .= "|$k";
                        }
                      }

                      push @{$intron_venn{$subset}} , "$seqname:".$sites{$node}->{Position}."-".$sites{$next_node}->{Position}.",".$sites{$node}->{Strand};

                    }
                }
            }
        }
    }
  return %intron_venn;
}









sub donor_venn {
  my %donor_venn;
  foreach my $seqname (keys %component_by_seqname)
    {
      foreach my $c (@{$component_by_seqname{$seqname}})
        {
          foreach my $node (@{$component{$c}})
            {
              if($node =~ /donor/)
                {

                  my %aux;
                  foreach my $source (@{$sites{$node}->{Source}})
                    {
                      $source =~ /(.+)?:(.+)/;
                      $aux{$1} = 1;
                    }

                  my $subset;
                  my $first = 1;
                  foreach my $k (sort {$a cmp $b} (keys %aux)) {
                    if($first) {
                      $first = 0;
                      $subset .= $k;
                    } else {
                      $subset .= "|$k";
                    }
                  }

                  push @{$donor_venn{$subset}} , "$seqname:".$sites{$node}->{Position}.",".$sites{$node}->{Strand};

                }
            }
        }
    }
  return %donor_venn;

}




sub acceptor_venn {
  my %acceptor_venn;
  foreach my $seqname (keys %component_by_seqname)
    {
      foreach my $c (@{$component_by_seqname{$seqname}})
        {
          foreach my $node (@{$component{$c}})
            {
              if($node =~ /acceptor/)
                {

                  my %aux;
                  foreach my $source (@{$sites{$node}->{Source}})
                    {
                      $source =~ /(.+)?:(.+)/;
                      $aux{$1} = 1;
                    }

                  my $subset;
                  my $first = 1;
                  foreach my $k (sort {$a cmp $b} (keys %aux)) {
                    if($first) {
                      $first = 0;
                      $subset .= $k;
                    } else {
                      $subset .= "|$k";
                    }
                  }

                  push @{$acceptor_venn{$subset}} , "$seqname:".$sites{$node}->{Position}.",".$sites{$node}->{Strand};

                }
            }
        }
    }
  return %acceptor_venn;

}



sub stop_codon_venn {
  my %stop_codon_venn;
  foreach my $seqname (keys %component_by_seqname)
    {
      foreach my $c (@{$component_by_seqname{$seqname}})
        {
          foreach my $node (@{$component{$c}})
            {
              if($node =~ /stop/)
                {

                  my %aux;
                  foreach my $source (@{$sites{$node}->{Source}})
                    {
                      $source =~ /(.+)?:(.+)/;
                      $aux{$1} = 1;
                    }

                  my $subset;
                  my $first = 1;
                  foreach my $k (sort {$a cmp $b} (keys %aux)) {
                    if($first) {
                      $first = 0;
                      $subset .= $k;
                    } else {
                      $subset .= "|$k";
                    }
                  }

                  push @{$stop_codon_venn{$subset}} , "$seqname:".$sites{$node}->{Position}.",".$sites{$node}->{Strand};

                }
            }
        }
    }
  return %stop_codon_venn;

}


sub start_codon_venn {
  my %start_codon_venn;
  foreach my $seqname (keys %component_by_seqname)
    {
      foreach my $c (@{$component_by_seqname{$seqname}})
        {
          foreach my $node (@{$component{$c}})
            {
              if($node =~ /start/)
                {

                  my %aux;
                  foreach my $source (@{$sites{$node}->{Source}})
                    {
                      $source =~ /(.+)?:(.+)/;
                      $aux{$1} = 1;
                    }

                  my $subset;
                  my $first = 1;
                  foreach my $k (sort {$a cmp $b} (keys %aux)) {
                    if($first) {
                      $first = 0;
                      $subset .= $k;
                    } else {
                      $subset .= "|$k";
                    }
                  }

                  push @{$start_codon_venn{$subset}} , "$seqname:".$sites{$node}->{Position}.",".$sites{$node}->{Strand};

                }
            }
        }
    }
  return %start_codon_venn;

}




sub process_forward {
    my $gene = shift;
    my $sites_r = shift;
    my $gtf_filename = shift;
#    my $geneid = $gene->id();
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
            $left_site = add_site($left_site, $sites_r, $source, $tx->id);
            $right_site = add_site($right_site, $sites_r, $source, $tx->id);

            if(defined $last_right_site)
            {
		
		push @{$last_right_site->{Next}->{get_key_from_site($left_site)}}, $source.":". $tx->id;
                push @{$left_site->{From}->{get_key_from_site($last_right_site)}}, $source.":". $tx->id;
            }
            push @{$left_site->{Next}->{get_key_from_site($right_site)}}, $source.":". $tx->id;
            push @{$right_site->{From}->{get_key_from_site($left_site)}}, $source.":". $tx->id;
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
            $left_site = add_site($left_site, $sites_r, $source, $tx->id);
            $right_site = add_site($right_site, $sites_r, $source, $tx->id);
            if(defined $last_right_site)
            {
                push @{$last_right_site->{Next}->{get_key_from_site($left_site)}}, $source.":".$tx->id;
                push @{$left_site->{From}->{get_key_from_site($last_right_site)}}, $source.":". $tx->id;
            }
            push @{$left_site->{Next}->{get_key_from_site($right_site)}}, $source.":". $tx->id;
            push @{$right_site->{From}->{get_key_from_site($left_site)}}, $source.":".$tx->id;
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
