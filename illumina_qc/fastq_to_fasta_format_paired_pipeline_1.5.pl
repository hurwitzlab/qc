#!/usr/bin/env perl 

use strict;
use warnings;
use IO::File;

die "$0: require [reads forward] [reads reverse] [outfile] as argument\n"
  if (scalar(@ARGV) != 3);

my $file_one = IO::File->new($ARGV[0], "r");
my $file_two = IO::File->new($ARGV[1], "r");
my $output_prefix = $ARGV[2];

my $done = 0;

my $outfile_fasta = IO::File->new($ARGV[2] . ".fasta", "w");
my $outfile_qual  = IO::File->new($ARGV[2] . ".qual",  "w");

my $identifier;
my $sequence;
my $qual;

while (!$done) {
    #alternate between the files, so the pairs are consecutive
    for my $file ($file_one, $file_two) {
        for (my $i = 0 ; $i < 4 ; $i++) {
            my $line = <$file>;
            chomp($line);
            if (    $i == 0
                and $line =~
                /^\@[^:]+:([0-8]):([0-9]+):([0-9]+):([0-9]+)#0\/([1-2])/gi)
            {
                my ($flowcell_lane, $tile, $x, $y, $readInPair) =
                  split(/:/, $_);

         #my($flowcell_lane, $tile ,$x, $y, $readInPair) = ($1, $2, $3, $4, $5);
                $identifier =
                  join("-", $flowcell_lane, $tile, $x, $y, $readInPair);
                print $outfile_fasta ">" . $identifier . "\n";
                print $outfile_qual ">" . $identifier . "\n";
            }
            elsif ($i == 1) {
                if ($line !~ /^[ATGCN]*$/gi) {
                    warn
"$identifier appears to have odd characters in sequence\n";
                }

                $sequence = $line;
                print $outfile_fasta $sequence . "\n";
            }
            elsif ($i == 2 and $line =~ /^\+/gi) {
                ;
            }
            elsif ($i == 3) {
                my @qual_chars = split(//, $line);
                my @qual_vals;
                for my $qual_char (@qual_chars) {
                    my $quality_score = toQual($qual_char);

                    if ($quality_score < 0 || $quality_score > 40) {
                        warn "Quality score out of range\n";
                    }
                    push(@qual_vals, $quality_score);
                }

                my $qual_string = join("-", @qual_vals);
                print $outfile_qual $qual_string . "\n";
            }
            else {

                warn "unexpected line in file: <$line>\n";
                die "Sequence: $identifier did not parse properly\n";
            }
        }
    }
    if (!(eof($file_one) and eof($file_two))) {
        $done = 0;
    }
    elsif (eof($file_one) and eof($file_two)) {
        $done = 1;
    }
    else {
        warn "Reached EOF for one file but not the other\n";
        $done = 1;
    }
}

sub toQual {
    my $quality_char = shift;

    my $quality_score = ord($quality_char) - 64;
    return $quality_score;
}
