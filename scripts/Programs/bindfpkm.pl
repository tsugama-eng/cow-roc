#!/usr/bin/perl

use strict;
use warnings;
use utf8;
use Getopt::Std;
my %opt;
getopts("g:d:p:F", \%opt);
#g is to set the gene column; d is to set the FPKM data (or data of interest) column; p is to set the pattern of file names; F is to set the field separator

#binmode STDIN, ':encoding(cp932)';
#binmode STDOUT, ':encoding(cp932)';
#binmode STDERR, ':encoding(cp932)';


if (!$opt{g}) {
	$opt{g} = 0;
} else {
	$opt{g} = $opt{g} - 1;
}

if (!$opt{d}) {
	$opt{d} = 1;
} else {
	$opt{d} = $opt{d} - 1;
}

if ($opt{g} == $opt{d}) {
	print STDERR "Options are wrong\nExitting\n";
	exit;
}

if(!$opt{p}) {
	$opt{p} = "*.fpkm";
}


if(!$opt{F}) {
	$opt{F} = "\t";
}

my %bounddata;
my $outfile = "SeqID";
my @files = glob ($opt{p});
@files = sort {$a cmp $b} @files;

for my $i (@files) {
	if ($i =~ /\.txt$/) {
		my $filename = $i;
		$filename =~ s/\.txt$//;
		$outfile = $outfile."\t".$filename;
	} else {
		$outfile = $outfile."\t".$i;
	}
	open(CF, "< $i");
	while (my $line=<CF>) {
		$line=~s/[\r\n]+//g;
		my @data = split($opt{F},$line);
		$bounddata{$data[$opt{g}]} = $bounddata{$data[$opt{g}]}."\t".$data[$opt{d}];
	}
	close (CF);
}

print STDOUT "$outfile\n";
for my $i (sort {$a cmp $b} keys %bounddata) {
	$bounddata{$i} =~ s/^\s+|\s+$//g;
	print STDOUT "$i\t$bounddata{$i}\n";
}

exit;
