#!/usr/bin/perl -w

use strict;
use warnings;

use GTF;
use Data::Dumper;
use Getopt::Long;

my $gtf_file;

GetOptions ("gtf=s" => \$gtf_file);

if(!defined $gtf_file) {
  print STDERR "USAGE: $0 -g <sampling.gtf>\n";
  exit();
}

# reading all GTFs
my %sites;
my %component;
my %transcripts;
my $gtf = GTF::new({gtf_filename => $gtf_file});
my $source = $gtf_file;
$source =~ s/\.gtf//g;
$source =~ s%.+/(.+)$%$1%g;

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


# build connected components
build_components(\%component, \%sites);

my %component_by_seqname;
foreach my $entry  (sort {$a <=> $b} keys %component)
  {
    ${$component{$entry}}[0] =~ m/(.+)?:.*/;
    push @{$component_by_seqname{$1}},$entry;
  }



my %genes = make_gene_cluster();

foreach my $seqname (keys %component_by_seqname)
  {
    my %transcript_to_freq;
    foreach my $entry (keys %{$genes{$seqname}}) {
      foreach my $gene (@{$genes{$seqname}{$entry}})
        {
          my ($atranscript, @other_transcripts)  = split(/;/, $gene);
          my $n = 1 + scalar(@other_transcripts);
          $transcript_to_freq{$atranscript} = ($n/scalar(@{$genes{$seqname}{$entry}}))*100.0;
        }
    }
    my $count;
    foreach my $entry (sort {my $xx = $transcript_to_freq{$b}; my $yy =  $transcript_to_freq{$a}; $xx <=> $yy } (keys %transcript_to_freq))
      {
        $entry =~ s/^.+?://g;
        $entry =~ s/,\d+$//g;
        print $transcripts{$entry}->output_gtf()."\n";
        $count ++;
        if($count >= 3){
          last;
        }
      }
  }




sub make_gene_cluster {
  my %subsets;
  my %gvenn;
  my %cluster;

  foreach my $seqname (keys %component_by_seqname)
    {
      my %supported_by;
      my %transcript_name;

      foreach my $component (@{$component_by_seqname{$seqname}})
        {
          foreach my $node (@{$component{$component}})
            {
              my $label = "";
              my $first = 1;
              foreach my $next_node  (keys %{$sites{$node}->{Next}})
                {
                  my @transcripts = @{$sites{$node}->{Next}->{$next_node}};
                  foreach my $source  (sort {$a cmp $b} (@transcripts))
                    {
                      $transcript_name{$source} = 1;
                      if($first ) {
                        $label .= "$source";
                        $first = 0;
                      }
                      else {
                        $label .= ";$source";
                      }
                    }
                  $supported_by{$label} = 1;
                }
            }

        }


      foreach my $tx_name (keys %transcript_name)
        {
          my $label_str;
          foreach my $label ((sort {$a cmp $b} (keys %supported_by))) {
            if(($label =~ m/$tx_name;/) || ($label =~ m/$tx_name$/) )
              {
                $label_str .= "<".$label.">";
              }
          }
          push @{$cluster{$seqname}{$label_str}}, $tx_name;
        }

      foreach my $key ( keys %{$cluster{$seqname}})
        {
          my $first = 1;
          my $list = "";
          foreach my $x (@{$cluster{$seqname}{$key}}) {
            if($x =~ m/(.+)?:(.+)/){
              my $nexon = count_exon($2);
              if($first) {
                $list .= "$x,$nexon";
                $first = 0;
              } else{
                $list .= ";$x,$nexon";
              }
            }
          }
          my $subset = build_subset_string($list);
          push @{$gvenn{$seqname}{$subset}}, $list;
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
          last;
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
