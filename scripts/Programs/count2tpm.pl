#!/usr/bin/perl

use strict;
use warnings;
use utf8;
use Getopt::Std;
my %opt;
getopts("g:c:F:l:", \%opt);
#g is to set the gene column; c is to set the count data column; F is to set the field separator; l is to set the gene length column

#binmode STDIN, ':encoding(cp932)';
#binmode STDOUT, ':encoding(cp932)';
#binmode STDERR, ':encoding(cp932)';

#comparison file containing blast results must be provided
#comparison file should be tab-delimited and should contain arabidopsis or rice gene names

if (!$opt{g}) {
	$opt{g} = 0;
} else {
	$opt{g} = $opt{g} - 1;
}

if (!$opt{c}) {
	$opt{c} = 1;
} else {
	$opt{c} = $opt{c} - 1;
}

if (!$opt{F}) {
	$opt{F} = "\t";
}

if ($opt{g} == $opt{c}) {
	print STDERR "Options are wrong\nExitting\n";
	exit;
}



my @gene;
my @count;
my %lgth;

if ($opt{l} =~ /^\d+$/) {
	$opt{l} = $opt{l} - 1;
} elsif ($opt{l}) {
	open (FF, "< $opt{l}") or die ("No fasta file found\nExitting\n");
	$/ = ">";
	while (my $line = <FF>) {
		my @data = split(/[\r\n]/,$line,2);
		$data[0] =~ s/ .*//g;
		$data[1] =~ s/\s//g;
		$lgth{$data[0]} = length($data[1]);
	}
	close (FF);
}

$/ = "\n";

while (my $line = <>) {
	$line =~ s/[\r\n]//g;
	my @data = split($opt{F},$line);
	push (@gene, $data[$opt{g}]);
	push (@count, $data[$opt{c}]);
	if (!$lgth{$data[$opt{g}]}) {
		$lgth{$data[$opt{g}]} = $data[$opt{l}];
	}
}

my $tsum = 0;

for my $i (0..$#gene) {
	if ($lgth{$gene[$i]} =~ /^\d+$/) {
		$tsum = $tsum + (1000 * $count[$i] / $lgth{$gene[$i]});
	}
}

for my $i (0..$#gene) {
	my $tpm;
	if ($lgth{$gene[$i]} =~ /^\d+$/) {
		$tpm = 1000000 * (1000 * $count[$i] / $lgth{$gene[$i]}) / $tsum;
		print STDOUT "$gene[$i]\t$tpm\n";
	}
}

exit;
