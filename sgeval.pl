#!/usr/bin/perl -w

use strict;
use warnings;

use GTF;
use Digest::MD5 qw(md5_hex);
use Data::Dumper;
use Getopt::Long;

my $print_graph = 0;
my @gtf_files;
my $output_dir;
my $compare  = 0;
GetOptions ("gtf=s{,}" => \@gtf_files,
            "out=s" => \$output_dir,
            "print_graph" => \$print_graph,
            "compare" => \$compare);

if($#gtf_files < 0 || !defined ($output_dir)) {
  print STDERR "USAGE: $0 [-c] [-p] -o <output_directory> -g <reference.gtf> <prediction1.gtf> <prediction2.gtf> ...\n";
  print STDERR "\t-c: venn diagrams without accuracy metrics\n";
  print STDERR "\t-p: print splicing graphs (tgf format)   \n";
  exit();
}
mkdir $output_dir;


#  hashes that stores md5 code for the entries.
my %names_hash;

my $tx_counter = 1;
foreach my $gtf_file (@gtf_files) {
    my $gtf = GTF::new({gtf_filename => $gtf_file});
    my $source = $gtf_file;
    $source =~ s/\.gtf//g;
    $source =~ s%.+/(.+)$%$1%g;
    $names_hash{md5_hex($source)} = $source;
    foreach my $gene (@{$gtf->genes()})
      {
        $names_hash{md5_hex($gene->seqname)} = $gene->seqname;
        foreach my $tx (@{$gene->transcripts}){
            my $txid= md5_hex($tx->id);
            if(defined ($names_hash{$txid})) {
              my $new_id = ($tx->id).".".$tx_counter;
              $tx_counter++;
              $tx->set_id($new_id);
              $txid = md5_hex($new_id);
            }
            $names_hash{$txid} = $tx->id;
        }
      }

  }

my %sites;
my %component;
my @sources;
my $ref_source;
my $first_source = 1;
my %transcripts;
$tx_counter = 1;
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
            my $txid= md5_hex($tx->id);
            if(defined ($transcripts{$txid})) {
              my $new_id = ($tx->id).".".$tx_counter;
              $tx_counter++;
              $tx->set_id($new_id);
              $txid = md5_hex($new_id);
            }
            $transcripts{$txid} = $tx;
          }
        if($gene->strand() eq "+")
          {
            process_forward($gene, \%sites, md5_hex($source));
          } else {
            process_reverse($gene, \%sites, md5_hex($source));
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

my %gvenn_exact = gene_exact_venn();
my %exon_exact = exon_exact_venn();
my %intron_exact = intron_exact_venn();
my %start = start_codon_venn();
my %stop = stop_codon_venn();
my %acceptor = acceptor_venn();
my %donor = donor_venn();

my %nucleotide = nucleotide_venn();
my %nucleotide_intron = nucleotide_intron_venn();
my %nucleotide_exon_with_intron_partial = nucleotide_venn_with_intron_partial();
# my %exon_overlaped = exon_overlaped_venn(\%nucleotide);

generate_result("nucleotide_exon", \%nucleotide);
generate_result("nucleotide_exon_with_intron_partial", \%nucleotide_exon_with_intron_partial);
generate_result("nucleotide_intron", \%nucleotide_intron);
# generate_result("exon_overlaped", \%exon_overlaped);

generate_result("gene_exact", \%gvenn_exact);
generate_result("exon_exact", \%exon_exact);
generate_result("intron_exact", \%intron_exact);

generate_result("start", \%start);
generate_result("stop", \%stop);
generate_result("acceptor", \%acceptor);
generate_result("donor", \%donor);

if($print_graph) {
  print_graph();
}

sub print_graph {
  mkdir ("$output_dir/splicing_graphs");
  foreach my $seqname (keys %component_by_seqname)
    {

      my $name =     $names_hash{$seqname};
      open (OUT, ">$output_dir/splicing_graphs/$name.tgf") or die "cant open $output_dir/splicing_graphs/$name.tgf:$!\n";
      my $id = 0;
      my %node_to_id;
      foreach my $c ( @{$component_by_seqname{$seqname}})
        {
          foreach my $node (@{$component{$c}})
            {
              $node_to_id{$node} = $id;
              my $x = $node;
              $x =~ s/((\w|\d){32}):(\d+,\w+)/$names_hash{$1}:$3/g;
              print OUT $id." ".$x."\n";
              $id++;
            }
        }
      print OUT "#\n";
      my %printed;
      foreach my $c ( @{$component_by_seqname{$seqname}})
        {
          foreach my $from (@{$component{$c}})
            {
              foreach my $to  (keys %{$sites{$from}->{Next}})
                {
                  my @transcripts = @{$sites{$from}->{Next}->{$to}};
                  foreach my $source  (sort {$a cmp $b} (@transcripts)){
                    my $x = $source;
                    $x =~ s/((\w|\d){32}):((\w|\d){32})/$names_hash{$1}:$names_hash{$3}/g;

                    print OUT $node_to_id{$from}." ".$node_to_id{$to}." ".$x."\n" if !(defined $printed{$node_to_id{$from}." ".$node_to_id{$to}." ".$x."\n"});
                    $printed{$node_to_id{$from}." ".$node_to_id{$to}." ".$x."\n"} = 1;


                  }
                }
            }
        }
      close(OUT);
    }
}


sub subset_string {
  my $aux_ref = shift;

  my %aux = %{$aux_ref};
  my $subset = "";
  my $first = 1;
  foreach my $k (sort {$a cmp $b} (keys %aux)) {
    if ($first) {
      $first = 0;
      $subset .= $k;
    } else {
      $subset .= "|$k";
    }
  }
  chomp($subset);
  return $subset;
}


sub generate_result {
  my $output_filename = shift;
  my $ref_venn = shift;
  my %venn = %{$ref_venn};
  open (OUTPUT, ">$output_dir/$output_filename"."_venn.txt") or die "$!";
  foreach my $key (keys %venn) {
    my $keyhash = $key;
    my @sources = split(/\|/, $keyhash);
    $key = $names_hash{$sources[0]};
    shift @sources;
    foreach  my $x (@sources) {
      $key .= "|".($names_hash{$x});
    }
      print OUTPUT $key."\t".$venn{$keyhash}->{"count"}."\n";
    foreach my $el ( @{$venn{$keyhash}->{"elements"}} )
      {
        if($el =~ /((\d|\w){32}):\d+,(-|\+)/){
          $el =~ s/((\d|\w){32}):(\d+),(-|\+)/<$names_hash{$1}>:$3,$4/;
        }elsif($el =~ /((\d|\w){32}):\d+-\d+,(-|\+),(\w+)/){
          $el =~ s/((\d|\w){32}):(\d+-\d+),(-|\+),(\w+)/<$names_hash{$1}>:$3,$4,$5/;
        }elsif($el =~ /((((\d|\w){32}):((\d|\w){32}),\d+);?)+/) {
          my @entries = split (/;/, $el);
          $el = $entries[0];
          $el =~ s/((\d|\w){32}):((\d|\w){32}),(\d+)/<$names_hash{$1}:$names_hash{$3},$5>/;
          shift @entries;
          foreach my $entry (@entries) {
            my $x = $entry;
            $x =~ s/((\d|\w){32}):((\d|\w){32}),(\d+)/<$names_hash{$1}:$names_hash{$3},$5>/;
            $el .= ";$x";
          }
        }elsif($el =~ /((\d|\w){32}):\d+-\d+,(-|\+)/){
          $el =~ s/((\d|\w){32}):(\d+-\d+),(-|\+)/<$names_hash{$1}>:$3,$4/;
        }

        print OUTPUT "\t".$el."\n";
      }
    print OUTPUT "//\n";
  }
  close(OUTPUT);
  if($compare == 0) {
  open (OUTPUT, ">$output_dir/$output_filename"."_accuracy.txt") or die "$!";
  foreach my $source (@sources) {
    my $ref_source2 = $ref_source;
    my $source2 = $source;
    my $source = md5_hex($source);
    my $ref_source = md5_hex($ref_source);
    my $tp = 0;
    my $fp = 0;
    my $fn = 0;
    if($source eq $ref_source) {
      next;
    }
    foreach my $subset (keys %venn)
      {

        my $count = $venn{$subset}->{"count"};
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


    my $sp = 0;

    if($tp + $fp != 0)
      {$sp = (100.0*($tp/($tp + $fp)))};
    my $sn = 0;
    if($tp + $fn != 0)
      { $sn = (100.0*($tp/($tp + $fn)));}

    my $f = 0;
    if(($sp + $sn) != 0) {
      $f = 2 * $sp * $sn / ($sp + $sn);
    }
    print OUTPUT $source2."\t".($tp+$fp)."\n";
    print OUTPUT "\tTP\t$tp\n\tFP\t$fp\n\tFN\t$fn\n";
    printf OUTPUT ("\tPPV\t%.5f\n\tSensitivity\t%.5f\n", $sp,$sn);
    printf OUTPUT ("\tF\t%.5f\n", $f);
    print OUTPUT "//\n";
  }
  close(OUTPUT);
 }
}


sub gene_exact_venn {
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
          push @{$cluster{$label_str}}, $tx_name;
        }
    }

  foreach my $key ( keys %cluster)
    {
      my $first = 1;
      my $list = "";
      foreach my $x (@{$cluster{$key}}) {
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
      push @{$gvenn{$subset}->{"elements"}}, $list;
      $gvenn{$subset}->{"count"} += 1;
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
            push @{$gvenn{$subset}->{"elements"}}, $list;
            $gvenn{$subset}->{"count"} += 1;

          }
        }

    }
  return %gvenn;
}


sub count_exon {
  my $id = shift;
  if($id) {
    my $tx = $transcripts{$id};
    if(!defined $tx) {
      die "ERROR: ".$id."\n";
    }
    return scalar(@{$tx->cds()});
  }
    return 0;
}


sub nucleotide_intron_venn {
  my %nucleotide_venn;

  foreach my $seqname (keys %component_by_seqname)
    {
      foreach my $c (@{$component_by_seqname{$seqname}})
        {
          my %intervals;
          foreach my $node (@{$component{$c}})
            {
              if(
                 ($node =~ /donor/ && $sites{$node}->{Strand} eq "+")
                 ||
                 ($node =~ /acceptor/ && $sites{$node}->{Strand} eq "-")
                )
                {
                  my $strand = $sites{$node}->{Strand};
                  foreach my $next_node (keys %{$sites{$node}->{Next}})
                    {
                      $node =~ m/(.+)?:(\d+),(.+)/;
                      my $start = $2;
                      $next_node =~ m/(.+)?:(\d+),(.+)/;
                      my $end = $2;
                      foreach my $source (@{$sites{$node}->{Next}->{$next_node}})
                        {
                          my $x;
                          $x->{"source"} = $source;
                          $x->{"start"} = $start;
                          $x->{"end"} = $end;
                          push @{$intervals{$strand}{$start}}, $x;
                          push @{$intervals{$strand}{$end}}, $x;
                        }
                    }
                }
            }

          foreach my $strand (keys %intervals)
            {
              my @sorted = sort {$a <=> $b}( keys %{$intervals{$strand}});
              my @sites;
              my @source_by_position;

              foreach my $p (@sorted)
                {
                  push @sites, $p;
                  if($#sites-1 >= 0 && $sites[$#sites] == $sites[$#sites-1])
                    {
                      pop @sites;
                    }
                }


              for (my $k = 0; $k <= $#sites; $k++)
                {
                  if($k > 0) {
                    foreach my $key (keys %{$source_by_position[$k-1]}){
                      ${$source_by_position[$k]}{$key} =  ${$source_by_position[$k-1]}{$key}  ;
                    }
                  }
                 foreach my $i (@{$intervals{$strand}{$sites[$k]}})
                    {
                      if($i->{"end"} == $sites[$k])
                        {
                          ${$source_by_position[$k]}{$i->{"source"}} = 0;
                        }
                    }
                  foreach my $i (@{$intervals{$strand}{$sites[$k]}})
                    {
                      if($i->{"start"} == $sites[$k])
                        {
                          ${$source_by_position[$k]}{$i->{"source"}} = 1;
                        }
                    }
                }
              my $last_end  = 0;
              for (my $p = 1; $p <= $#source_by_position; $p+=1) {
                my %left_aux;
                foreach my $source (keys %{$source_by_position[$p-1]}) {
                  if(${$source_by_position[$p-1]}{$source} == 1) {
                    $source =~ /(.+)?:(.+)/;

                    $left_aux{$1} = 1;
                  }
                }

                my %right_aux;
                foreach my $source (keys %{$source_by_position[$p]}) {
                  if(${$source_by_position[$p]}{$source} == 1) {
                    $source =~ /(.+)?:(.+)/;
                    $right_aux{$1} = 1;
                  }
                }

                my $subset = subset_string(\%left_aux);
                my $subset2 = subset_string(\%right_aux);

                if(!$subset eq ""){
                  my $start = $sites[$p-1];
                  my $end = $sites[$p];

		   if(($p-2 >= 0) && ($start == $last_end))
		    {
		      $start += 1;
		    }
                  $last_end = $end;

                  push @{$nucleotide_venn{$subset}->{"elements"}},$seqname.":".$start."-".$end.",".$strand.",".($end - $start + 1);
                  push @{$nucleotide_venn{$subset}->{"interval"}->{$seqname}->{$strand}},$start."-".$end;
                  $nucleotide_venn{$subset}->{"count"}  += $end - $start + 1;
                }
              }





            }

        }
    }
  return %nucleotide_venn;
}




sub nucleotide_venn {
  my %nucleotide_venn;

  foreach my $seqname (keys %component_by_seqname)
    {
      foreach my $c (@{$component_by_seqname{$seqname}})
        {
          my %intervals;

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
                      foreach my $source (@{$sites{$node}->{Next}->{$next_node}})
                        {
                          $source =~ /(.+)?:(.+)/;
                          my $x;
                          $x->{"source"} = $source;
                          $x->{"start"} = $start;
                          $x->{"end"} = $end;
                          push @{$intervals{$strand}{$start}}, $x;
                          push @{$intervals{$strand}{$end}}, $x;
                        }
                    }
                }
            }
          foreach my $strand (keys %intervals)
            {
              my @sorted = sort {$a <=> $b}( keys %{$intervals{$strand}});
              my @sites;
              my @source_by_position;

              #  uniq
              foreach my $p (@sorted)
                {
                  push @sites, $p;
                  if($#sites-1 >= 0 && $sites[$#sites] == $sites[$#sites-1])
                    {
                      pop @sites;
                    }
                }
		for (my $k = 0; $k <= $#sites; $k++)
                {
                  if($k > 0) {
                    foreach my $key (keys %{$source_by_position[$k-1]}){
                      ${$source_by_position[$k]}{$key} =  ${$source_by_position[$k-1]}{$key}  ;
                    }
                  }
                  foreach my $i  (@{$intervals{$strand}{$sites[$k]}})
                    {
                      if($i->{"end"} == $sites[$k])
                        {
                          ${$source_by_position[$k]}{$i->{"source"}} = 0;
                        }
                    }

                  foreach my $i  (@{$intervals{$strand}{$sites[$k]}})
                    {
                      if($i->{"start"} == $sites[$k])
                        {
                          ${$source_by_position[$k]}{$i->{"source"}} = 1;
                        }
                    }
                }
		my %already_counted;
              for (my $p = 1; $p <= $#source_by_position; $p+=1) {
                my %left_aux;
                foreach my $source (keys %{$source_by_position[$p-1]}) {
                  if(${$source_by_position[$p-1]}{$source} == 1) {
                          $source =~ /(.+)?:(.+)/;
                          $left_aux{$1} = 1;
                  }
                }

                my %right_aux;
                foreach my $source (keys %{$source_by_position[$p]}) {
                  if(${$source_by_position[$p]}{$source} == 1) {
                    $source =~ /(.+)?:(.+)/;
                    $right_aux{$1} = 1;
                  }
                }

                my $subset = subset_string(\%left_aux);
                my $subset2 = subset_string(\%right_aux);
                if(!$subset eq ""){
                  my $start = $sites[$p-1];
                  my $end = $sites[$p];
		  if ($subset2 !~ /^\s*$/ && $subset =~ /$subset2/g && length($subset) < length($subset2)) {
		    $end--;
		  }
		  if (defined $already_counted{$subset.$start}) {
		    $start ++;
		  }
		  $already_counted{$subset.$start} = 1;
		  $already_counted{$subset2.$end} = 1;
                  push @{$nucleotide_venn{$subset}->{"elements"}},$seqname.":".$start."-".$end.",".$strand.",".($end - $start + 1);
                  push @{$nucleotide_venn{$subset}->{"interval"}->{$seqname}->{$strand}},$start."-".$end;
                  $nucleotide_venn{$subset}->{"count"}  += $end - $start + 1;
                }
              }
            }

        }
    }
  return %nucleotide_venn;
}



sub nucleotide_venn_with_intron_partial {
  my %nucleotide_venn;
  foreach my $seqname (keys %component_by_seqname)
    {
      foreach my $c (@{$component_by_seqname{$seqname}})
        {
          my %intervals;
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
                      foreach my $source (@{$sites{$node}->{Next}->{$next_node}})
                        {
                          my $x;
                          $x->{"source"} = $source;
                          $x->{"start"} = $start;
                          $x->{"end"} = $end;
                          push @{$intervals{$strand}{$start}}, $x;
                          push @{$intervals{$strand}{$end}}, $x;
                        }
                    }
                }
              if(
                 ($node =~ /donor/ && $sites{$node}->{Strand} eq "+")
                 ||
                 ($node =~ /acceptor/ && $sites{$node}->{Strand} eq "-")
                )
                {
                  my $strand = $sites{$node}->{Strand};
                  foreach my $next_node (keys %{$sites{$node}->{Next}})
                    {
                      $node =~ m/(.+)?:(\d+),(.+)/;
                      my $start = $2;
                      $next_node =~ m/(.+)?:(\d+),(.+)/;
                      my $end = $2;
                      foreach my $source (@{$sites{$node}->{Next}->{$next_node}})
                        {
			  my $x;
                          my $partA = $source;
                          my $partB;
                          if($source =~ /:/) {
                            ($partA, $partB) = split(/:/, $source);
                          }
                          if($partA eq md5_hex($ref_source)) {
                            my $clean_name = $names_hash{$partA};
                            my $clean_intron_name = $clean_name."_intron";
                            $x->{"source"} = $source;
                            $x->{"start"} = $start + 1;
                            $x->{"end"} = $end - 1;
                            $x->{"source"} = md5_hex($clean_intron_name);
			    $names_hash{$x->{"source"}} = $clean_intron_name; 
                            push @{$intervals{$strand}{$start + 1}}, $x;
                            push @{$intervals{$strand}{$end - 1}}, $x;
                          }
			  
                        }
                    }
                }
            }


          foreach my $strand (keys %intervals)
            {
              my @sorted = sort {$a <=> $b}( keys %{$intervals{$strand}});
              my @sites;
              my @source_by_position;

              foreach my $p (@sorted)
                {
                  push @sites, $p;
                  if($#sites-1 >= 0 && $sites[$#sites] == $sites[$#sites-1])
                    {
                      pop @sites;
                    }
                }


              for (my $k = 0; $k <= $#sites; $k++)
                {
                  if($k > 0) {
                    foreach my $key (keys %{$source_by_position[$k-1]}){
                      ${$source_by_position[$k]}{$key} =  ${$source_by_position[$k-1]}{$key}  ;
                    }
                  }
                  foreach my $i (@{$intervals{$strand}{$sites[$k]}})
                    {
                      if($i->{"start"} == $sites[$k])
                        {
                          ${$source_by_position[$k]}{$i->{"source"}} = 1;
                        }
                      if($i->{"end"} == $sites[$k])
                        {
                          ${$source_by_position[$k]}{$i->{"source"}} = 0;
                        }
                    }
                }

		my %already_counted;

              for (my $p = 1; $p <= $#source_by_position; $p+=1) {
                my %left_aux;

		foreach my $source (keys %{$source_by_position[$p-1]}) {
                  if(${$source_by_position[$p-1]}{$source} == 1) {
                    my $partA = $source;
                    my $partB;
                    if($source =~ /:/) {
                      ($partA, $partB) = split(/:/, $source);
                    }
                    my $n = $names_hash{$partA};
                    $names_hash {md5_hex($n)} = $n;
                    $left_aux{md5_hex($n)} = 1;
                  }
                }
                my %right_aux;
                foreach my $source (keys %{$source_by_position[$p]}) {
                  if(${$source_by_position[$p]}{$source} == 1) {
                    my $partA = $source;
                    my $partB;
                    if($source =~ /:/) {
                      ($partA, $partB) = split(/:/, $source);
                    }
                    my $n = $names_hash{$partA};
                    $names_hash{md5_hex($n)} = $n;
                    $right_aux{md5_hex($n)} = 1;
                  }
                }
		
                my $subset = subset_string(\%left_aux);
                my $subset2 = subset_string(\%right_aux);
                 
                if(!$subset eq ""){
                  my $start = $sites[$p-1];
                  my $end = $sites[$p];
		  # print "<$subset> <$subset2> <$start> <$end>\n";
		  if ($subset2 !~ /^\s*$/ && $subset =~ /$subset2/g && length($subset) < length($subset2)) {
		    $end--;
		  }
		  if (defined $already_counted{$subset.$start}) {
		    $start ++;
		  }
		  $already_counted{$subset.$start} = 1;
		  $already_counted{$subset2.$end} = 1;
                
                    push @{$nucleotide_venn{$subset}->{"elements"}},$seqname.":".$start."-".$end.",".$strand.",".($end - $start + 1);
                    push @{$nucleotide_venn{$subset}->{"interval"}->{$seqname}->{$strand}},$start."-".$end;
                    $nucleotide_venn{$subset}->{"count"}  += $end - $start + 1;
                }
              }
            }
        }
    }
  return %nucleotide_venn;
}




              
sub exon_overlaped_venn {
  my $nucleotide_venn_ref = shift;
  my %nucleotide_venn = %{$nucleotide_venn_ref};
  my %exon_overlaped_venn;
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
                  my $strand = $sites{$node}->{Strand};
                  foreach my $next_node (keys %{$sites{$node}->{Next}})
                    {
                      $node =~ m/(.+)?:(\d+),(.+)/;
                      my $start = $2;
                      $next_node =~ m/(.+)?:(\d+),(.+)/;
                      my $end = $2;
                      my $strand = $sites{$node}->{Strand};
                      my %aux;
                      foreach my $subset1 (keys %nucleotide_venn)
                        {
                          my $is_overlap = 0;
                          foreach my $interval (@{$nucleotide_venn{$subset1}->{"interval"}->{$seqname}->{$strand}})
                            {
                              my ($start2, $end2) = split(/-/, $interval);
#                              print "<".$subset1."> ".$start2." ".$end2." ".$start." ".$end."\n";
                              if((($start2 >= $start) && ($start2 <= $end))||
                                 (($start >= $start2) && ($start <= $end2)))
                                {
#                                  print " overlaped\n";
                                  $is_overlap = 1;
                                  last;
                                }
                            }
                          if($is_overlap) {
                            foreach my $set (split(/\|/, $subset1))
                              {
                                $aux{$set} = 1;
                              }
                          }
                        }

                      my $subset = subset_string(\%aux);

#                      print "<".$subset.">\n";
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

                      push @{$exon_overlaped_venn{$subset}->{"elements"}} , $seqname.":".$sites{$node}->{Position}."-".$sites{$next_node}->{Position}.",".$sites{$node}->{Strand}.",".$type;
                      $exon_overlaped_venn{$subset}->{"count"} +=1;
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
                      my $subset = subset_string(\%aux);
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

                      push @{$exon_venn{$subset}->{"elements"}} , $seqname.":".$sites{$node}->{Position}."-".$sites{$next_node}->{Position}.",".$sites{$node}->{Strand}.",".$type;
                      $exon_venn{$subset}->{"count"} += 1;


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

                      my $subset = subset_string(\%aux);


                      push @{$intron_venn{$subset}->{"elements"}} , $seqname.":".$sites{$node}->{Position}."-".$sites{$next_node}->{Position}.",".$sites{$node}->{Strand};
                      $intron_venn{$subset}->{"count"} += 1;

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

                  my $subset = subset_string(\%aux);
                  push @{$donor_venn{$subset}->{"elements"}} , $seqname.":".$sites{$node}->{Position}.",".$sites{$node}->{Strand};
                  $donor_venn{$subset}->{"count"} += 1;
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
                  my $subset = subset_string(\%aux);

                  push @{$acceptor_venn{$subset}->{"elements"}} , $seqname.":".$sites{$node}->{Position}.",".$sites{$node}->{Strand};
                  $acceptor_venn{$subset} ->{"count"} += 1;
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

                  my $subset = subset_string(\%aux);

                  push @{$stop_codon_venn{$subset}->{"elements"}} , $seqname.":".$sites{$node}->{Position}.",".$sites{$node}->{Strand};
                  $stop_codon_venn{$subset}->{"count"} += 1;
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

                  my $subset = subset_string(\%aux);

                  push @{$start_codon_venn{$subset}->{"elements"}} , $seqname.":".$sites{$node}->{Position}.",".$sites{$node}->{Strand};
                  $start_codon_venn{$subset}->{"count"} += 1;
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
                    SeqName => md5_hex($gene->seqname()),
                    Position => $cds->start(),
                    Type=> "start_codon",
                    Frame => $cds->frame(),
                    Strand => "+"
                };
            } else {
                $left_site =
                {
                    SeqName =>md5_hex($gene->seqname()),
                    Position => $cds->start(),
                    Type => "acceptor",
                    Frame => $cds->frame(),
                    Strand => "+"
                }
            }
            if($is_stop_codon) {
                $right_site =
                {
                    SeqName =>md5_hex($gene->seqname()),
                    Position => $cds->stop(),
                    Type => "stop_codon",
                    Frame => $cds->frame(),
                    Strand => "+"
                };
            } else {
                $right_site =
                {
                    SeqName =>md5_hex($gene->seqname()),
                    Position => $cds->stop(),
                    Type => "donor",
                    Frame => $cds->frame(),
                    Strand => "+"
                };
            }
            my $txid = md5_hex($tx->id);
            my $source = $gtf_filename;
            $left_site = add_site($left_site, $sites_r, $source, $txid);
            $right_site = add_site($right_site, $sites_r, $source, $txid);

            if(defined $last_right_site)
            {

                push @{$last_right_site->{Next}->{get_key_from_site($left_site)}}, $source.":". $txid;
                push @{$left_site->{From}->{get_key_from_site($last_right_site)}}, $source.":". $txid;
            }
            push @{$left_site->{Next}->{get_key_from_site($right_site)}}, $source.":". $txid;
            push @{$right_site->{From}->{get_key_from_site($left_site)}}, $source.":". $txid;
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
                    SeqName => md5_hex($gene->seqname()),
                    Position => $cds->stop(),
                    Type=> "start_codon",
                    Frame => $cds->frame(),
                    Strand => "-"
                };
            } else {
                $right_site =
                {
                    SeqName =>md5_hex($gene->seqname()),
                    Position => $cds->stop(),
                    Type => "acceptor",
                    Frame => $cds->frame(),
                    Strand => "-"
                }
            }
            if($is_stop_codon) {
                $left_site =
                {
                    SeqName =>md5_hex($gene->seqname()),
                    Position => $cds->start(),
                    Type => "stop_codon",
                    Frame => $cds->frame(),
                    Strand => "-"
                };
            } else {
                $left_site =
                {
                    SeqName =>md5_hex($gene->seqname()),
                    Position => $cds->start(),
                    Type => "donor",
                    Frame => $cds->frame(),
                    Strand => "-"
                };
            }
            my $source = $gtf_filename;
            my $txid = md5_hex($tx->id);
            $left_site = add_site($left_site, $sites_r, $source, $txid);
            $right_site = add_site($right_site, $sites_r, $source, $txid);
            if(defined $last_right_site)
            {
                push @{$last_right_site->{Next}->{get_key_from_site($left_site)}}, $source.":".($txid);
                push @{$left_site->{From}->{get_key_from_site($last_right_site)}}, $source.":". ($txid);
            }
            push @{$left_site->{Next}->{get_key_from_site($right_site)}}, $source.":".($txid);
            push @{$right_site->{From}->{get_key_from_site($left_site)}}, $source.":".($txid);
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

