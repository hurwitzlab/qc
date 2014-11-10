#!/usr/bin/env perl 

use strict;
use warnings;
use IO::File;

die "Require file as argument\n" if (scalar(@ARGV) != 1);

my $infile          = IO::File->new($ARGV[0], "r");
my @fields          = ();
my $prev_identifier = "";

#
# using the output of the analyse Illumina quality scores script, 
# the average quality data for mate-pair reads are collapsed into
# one line.
#

while (my $line = <$infile>) {
    chomp($line);
    my ($a, $b, $c, $d, $avg_quality) = split(/\t/, $line);

    my $identifier = join("-", ($a, $b, $c, $d));

    if (scalar(@fields) == 2) {
        my @id_fields = split("-", $prev_identifier);
        print join("\t", @id_fields), "\t", join("\t", @fields), "\n";

        @fields = ();
    }
    $prev_identifier = $identifier;
    push(@fields, $avg_quality);
}

if (scalar(@fields) == 2) {
    my @id_fields = split("-", $prev_identifier);
    print join("\t", @id_fields), "\t", join("\t", @fields), "\n";

}
