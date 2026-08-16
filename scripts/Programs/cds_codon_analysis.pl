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

$/ = ">";

#my $pattern = qr/\.\d+$/;
#if ($ARGV[0] =~ /^Populus_/){$pattern = qr/\.\d+?\.v4\.1$/;}
#if ($ARGV[0] =~ /^Saccharomyces_/){$pattern = qr/_mRNA$/;}
#if ($ARGV[0] =~ /^Vitis_/){$pattern = qr/\.t\d+$/;}
#if ($ARGV[0] =~ /^Zea_/){$pattern = qr/_T\d+$/;}

#my $header = "gene\tAAC\tAAG\tACC\tACG\tAGC\tAGG\tATC\tATG\tCAC\tCAG\tCCC\tCCG\tCGC\tCGG\tCTC\tCTG\tGAC\tGAG\tGCC\tGCG\tGGC\tGGG\tGTC\tGTG\tTAC\tTAG\tTCC\tTCG\tTGC\tTGG\tTTC\tTTG\tAAA\tAAT\tACA\tACT\tAGA\tAGT\tATA\tATT\tCAA\tCAT\tCCA\tCCT\tCGA\tCGT\tCTA\tCTT\tGAA\tGAT\tGCA\tGCT\tGGA\tGGT\tGTA\tGTT\tTAA\tTAT\tTCA\tTCT\tTGA\tTGT\tTTA\tTTT\tGC_content\tGC3_content"; 
#if (!$rscu) {print "$header\n";} #A different header is used for the rscu option
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
my $header = join("\t", "gene", @codons, "GC_content", "GC3_content");
print "$header\n";

my $all_seqlen = 0;
my $all_gc_count = 0;
my %all_codon_count;
my $all_gc3count = 0;
my $all_at3count = 0;

my %seqflag;

while (my $seqdata = <>) {
  my ($seqid, $seq) = split(/[\r\n]+/, $seqdata, 2);
  if(!$seq){next;}
  $seq =~ s/[>\s\r\n]//g;
  $seq = uc($seq);  # 大文字に統一
  if ($seq !~ /[ACGTU]/ || length($seq) % 3 != 0) {next;}
  $seqid =~ s/.* (gene:[^ ]+) .*/$1/g;
  $seqid =~ s/ .*//g;
  if ($seqflag{$seqid}) {next;}
  $seqflag{$seqid}++;
  #$seqid =~ s/$pattern//;
  my $gc_count = ($seq =~ tr/GC//);
  my $gc_content = $gc_count / length($seq);
  $all_seqlen = $all_seqlen + length($seq);
  $all_gc_count = $all_gc_count + $gc_count;
  my %codon_count;
  for my $i (@codons) {$codon_count{$i} = 0;}
  my $gc3count = 0;
  my $at3count = 0;
  my $gc3content = 0;
  # 3塩基ごとに分割
  for (my $i = 0; $i < length($seq) - 2; $i += 3) {
    my $codon = substr($seq, $i, 3);
    next if $codon =~ /[^ACGTU]/;  # 不明塩基はスキップ
    $codon_count{$codon}++;
    $all_codon_count{$codon}++;
    if ($codon =~ /[ACGTU][ACGTU][GC]/){$gc3count++; $all_gc3count++;}
    if ($codon =~ /[ACGTU][ACGTU][ATU]/){$at3count++; $all_at3count++;}
  }
  my $total_codons = $gc3count + $at3count;
  if ($total_codons > 0){
    $gc3content = $gc3count / $total_codons;
  } else {next;}

  # 出力
  print "$seqid";
  if (!$normalized && !$rscu) {
    for my $codon (@codons) {print "\t$codon_count{$codon}";}
  } elsif ($normalized) {
    # normalized frequency 出力
    #my $total_codons = 0;
    #$total_codons += $codon_count{$_} for @codons;
    for my $codon (@codons) {
      my $freq = ($total_codons > 0) ? $codon_count{$codon} / $total_codons : 0;
      #printf "%.4f\t", $freq;
      print "\t$freq";
    }
  } elsif ($rscu) {
    # RSCU 出力
    for my $aa (sort keys %synonymous) {
      my @codons_for_aa = @{$synonymous{$aa}};
      my $sum = 0;
      $sum += $codon_count{$_} for @codons_for_aa;
      my $n = scalar @codons_for_aa;
      for my $codon (@codons_for_aa) {
        my $val = ($sum > 0) ? ($codon_count{$codon} * $n / $sum) : 0;
        #printf "%.3f\t", $val;
        print "\t$val";
      }
    }
  }
  print "\t$gc_content\t$gc3content\n";
}

my $all_gc_content = $all_gc_count / $all_seqlen;
my $all_total_codons = $all_gc3count + $all_at3count;
my $all_gc3content = $all_gc3count / $all_total_codons;

print "all";
if (!$normalized && !$rscu){
  for my $codon (@codons) {print "\t$all_codon_count{$codon}";}
} elsif ($normalized) {
  for my $codon (@codons) {
    my $freq = ($all_total_codons > 0) ? $all_codon_count{$codon} / $all_total_codons : 0;
    #printf "%.4f\t", $freq;
    print "\t$freq";
  }
} elsif ($rscu) {
  # RSCU 出力
  for my $aa (sort keys %synonymous) {
    my @codons_for_aa = @{$synonymous{$aa}};
    my $sum = 0;
    $sum += $all_codon_count{$_} for @codons_for_aa;
    my $n = scalar @codons_for_aa;
    for my $codon (@codons_for_aa) {
      my $val = ($sum > 0) ? ($all_codon_count{$codon} * $n / $sum) : 0;
      #printf "%.3f\t", $val;
      print "\t$val";
    }
  }
}
print "\t$all_gc_content\t$all_gc3content\n";
