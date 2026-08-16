#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;

my $normalized = 0;
my $rscu = 0;

GetOptions(
    "n|normalized" => \$normalized,  # -n または --normalized
    "r|rscu"       => \$rscu,        # -r または --rscu
) or die "Error in command line arguments\n";

#The first and second arguments should be the CDS-derived codon data file and the TPM file, respectively
#The data in the first columns of the input files should be consistent with each other

my %codondata;
my %tpmdata;

#if ($ARGV[0] !~ /Oryza_sativa/){die;}

#my $header = "sample\tAAC\tAAG\tACC\tACG\tAGC\tAGG\tATC\tATG\tCAC\tCAG\tCCC\tCCG\tCGC\tCGG\tCTC\tCTG\tGAC\tGAG\tGCC\tGCG\tGGC\tGGG\tGTC\tGTG\tTAC\tTAG\tTCC\tTCG\tTGC\tTGG\tTTC\tTTG\tAAA\tAAT\tACA\tACT\tAGA\tAGT\tATA\tATT\tCAA\tCAT\tCCA\tCCT\tCGA\tCGT\tCTA\tCTT\tGAA\tGAT\tGCA\tGCT\tGGA\tGGT\tGTA\tGTT\tTAA\tTAT\tTCA\tTCT\tTGA\tTGT\tTTA\tTTT\tGC3_content"; 
#print "$header\n";
#my @codons = split ("\t", $header);
#@codons = @codons[1..64];

my %synonymous = (
    'Phe' => ['TTT','TTC'],                          # フェニルアラニン
    'Leu' => ['TTA','TTG','CTT','CTC','CTA','CTG'],  # ロイシン
    'Ile' => ['ATT','ATC','ATA'],                    # イソロイシン
    'Met' => ['ATG'],                                # メチオニン (開始コドン)
    'Val' => ['GTT','GTC','GTA','GTG'],              # バリン
    'Ser' => ['TCT','TCC','TCA','TCG','AGT','AGC'],  # セリン
    'Pro' => ['CCT','CCC','CCA','CCG'],              # プロリン
    'Thr' => ['ACT','ACC','ACA','ACG'],              # スレオニン
    'Ala' => ['GCT','GCC','GCA','GCG'],              # アラニン
    'Tyr' => ['TAT','TAC'],                          # チロシン
    'His' => ['CAT','CAC'],                          # ヒスチジン
    'Gln' => ['CAA','CAG'],                          # グルタミン
    'Asn' => ['AAT','AAC'],                          # アスパラギン
    'Lys' => ['AAA','AAG'],                          # リシン
    'Asp' => ['GAT','GAC'],                          # アスパラギン酸
    'Glu' => ['GAA','GAG'],                          # グルタミン酸
    'Cys' => ['TGT','TGC'],                          # システイン
    'Trp' => ['TGG'],                                # トリプトファン
    'Arg' => ['CGT','CGC','CGA','CGG','AGA','AGG'],  # アルギニン
    'Gly' => ['GGT','GGC','GGA','GGG'],              # グリシン
    'Stop' => ['TAA','TAG','TGA'],                   # 終止コドン
);

my @codons;
for my $aa (sort keys %synonymous) {
    push @codons, @{$synonymous{$aa}};
}
my $header = join("\t", "gene", @codons, "GC3_content");
print "$header\n";

open (CF, "< $ARGV[0]") or die;
while (my $line = <CF>) {
  $line =~ s/[\r\n]+//g;
  my ($gene, $curcodondata) = split ("\t", $line, 2);
  if ($gene eq "gene" || $gene eq "all") {next;}
  $codondata{$gene} = $curcodondata;
}
close (CF);

open (TF, "< $ARGV[1]") or die;
while (my $line = <TF>) {
  $line =~ s/[\r\n]+//g;
  my ($gene, $tpms) = split ("\t", $line, 2);
  $tpmdata{$gene} = $tpms;
}
close (TF);

my @samples = split ("\t", $tpmdata{"SeqID"});
for my $i (0..$#samples) {
  #my $filename = $samples[$i];
  #$filename =~ s/_gene.tpm//;
  #$filename = $filename."_codon_demand.txt";
  #open (OF, "> $filename");
  my %tpm; #keys should be gene IDs
  my %demand; #keys should be codons
  for my $j (@codons) {$demand{$j} = 0;}
  my %codonidx; #keys should be codons
  for my $j (0..$#codons){$codonidx{$codons[$j]} = $j;}
  for my $j (keys %codondata) {
    if (!$tpmdata{$j}){next;}
    my @tpms = split ("\t", $tpmdata{$j});
    $tpm{$j} = $tpms[$i];
    my @codonnums = split ("\t", $codondata{$j});
    for my $k (keys %codonidx) {
      if ($tpm{$j}) {$demand{$k} += $tpm{$j} * ($codonnums[$codonidx{$k}] || 0);}
    }
  }
  $samples[$i] =~ s/_gene\.tpm//;
  print "$samples[$i]";
  my ($gc3demand, $at3demand) = (0, 0);
  for my $j (@codons) {
    #print "\t$demand{$j}";
    if ($j =~ /[ACGTU][ACGTU][CG]/){$gc3demand = $gc3demand + $demand{$j};}
    if ($j =~ /[ACGTU][ACGTU][ATU]/){$at3demand = $at3demand + $demand{$j};}
  }
  my $gc3content = 0;
  my $total_demand = $gc3demand + $at3demand;
  if ($total_demand > 0) {
    $gc3content = $gc3demand / $total_demand;
  }
  if (!$normalized && !$rscu) {
    for my $j (@codons) {print "\t$demand{$j}";}
  } elsif ($normalized) {
    for my $j (@codons) {
      my $normalized_demand = 0;
      if ($total_demand > 0) {$normalized_demand = $demand{$j} / $total_demand;}
      print "\t$normalized_demand";
    }
  } elsif ($rscu) {
    for my $aa (sort keys %synonymous) {
      my @codons_for_aa = @{$synonymous{$aa}};
      my $sum = 0;
      $sum += $demand{$_} for @codons_for_aa;
      my $n = scalar @codons_for_aa;
      for my $codon (@codons_for_aa) {
        my $val = ($sum > 0) ? ($demand{$codon} * $n / $sum) : 0;
        #printf "%.3f\t", $val;
        print "\t$val";
      }
    }
  }
  print "\t$gc3content\n";
  #if ($i == 0) {last;}
}

