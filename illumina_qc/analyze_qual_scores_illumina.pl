#!/usr/bin/env perl

use strict;
use warnings;
use IO::File;

die "$0: require [qualfile] [outfile prefix] [subset size]  as arguments\n"
  if (scalar(@ARGV) != 3);

my $max_seqs    = `grep -c ">" $ARGV[0]`;
my $qual_fh     = IO::File->new($ARGV[0], "r");
my $subset_size = $ARGV[2];
my %seqs_to_sample;
my $count = 0;

if ($subset_size > $max_seqs) {
    print "Sample larger than file\n";
    exit;
}

while ($count < $subset_size) {
    my $random_pos = int(rand($max_seqs));

    if (!exists $seqs_to_sample{$random_pos}) {
        $seqs_to_sample{$random_pos} = 1;
        $count++;
    }
}

if (-e $ARGV[1]) {
    die "File " . ($ARGV[1]) . " already exists\n";
}

my $avg_fh       = IO::File->new($ARGV[1] . ".avg_qual",      "w");
my $qual_base_fh = IO::File->new($ARGV[1] . ".qual_per_base", "w");
my $id;

my ($lane, $x, $y, $z, $readInPair);

my $seq_count = 0;

while (my $line = <$qual_fh>) {
    chomp($line);
    if ($line =~ /^>/gi) {
        ($lane, $x, $y, $z, $readInPair) = split("-", $line);

        $lane =~ s/>//gi;
    }
    elsif ($line !~ /^>/gi) {

        my @scores = split(/-/, $line);

        if (exists $seqs_to_sample{$seq_count}) {
            print $qual_base_fh join("\t", @scores), "\n";

        }
        my $sum = 0;
        foreach my $score (@scores) {
            $sum += $score;
        }

        my $avg = $sum / (scalar @scores);
        print $avg_fh join("\t", ($x, $y, $z, $readInPair, $avg)), "\n";

        $seq_count++;
    }
}
