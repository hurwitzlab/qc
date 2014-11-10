#!/usr/bin/env perl 

use strict;
use warnings;
use IO::File;
use Getopt::Std;
use Readonly;
use vars qw/ %opt /;

Readonly my $DEFAULT_WINDOW_SIZE     => 5;
Readonly my $DEFAULT_WINDOW_MIN_QUAL => 20;
Readonly my $DEFAULT_WINDOW_AVG_QUAL => 20;
Readonly my $DEFAULT_MIN_LENGTH      => 0;
#special case, 2 means no observation!
Readonly my $QUAL_UNOBSERVED_NUCLEOTIDE => 2; 

my $count_avg_qual_low    = 0;
my $count_trimmed_to_zero = 0;
my $count_contained_ns    = 0;
my $too_short             = 0;

#
# This script performs QC on Illumina paired-reads using the quality
# scores (assuming reads were pre-processed by an in-house
# fastq->fasta script). # Reads are initially end-trimmed using
# sliding window of size -s, with windows trimmed that are less than
# -r in avg quality. 			    # Finally, the average quality is
# calculated over the trimmed read, and the read is thrown out if the
# average quality is less than -a.		    # The trimmed reads are
# then filtered by length. Output is a mate-pair and singleton file.
# 
# 

sub usage() {
    print
"Usage: $0 -f [fasta file] -q [qual file] -l [min length] -s [window size] -r [window min qual] -a [min avg qual] -p [outfile prefix]\n";
    print
"-f fasta file: contains the paired-reads (assuming consecutive in file forward followed by reverse) in fasta format \n";
    print
"-q quality file: contains the quality scores for the paired-reads in Phred format (separated by dashes)\n";
    print
"-l integer <optional>: all trimmeds less than this length will be thrown out. Default 0.\n";
    print
"-s window size <optional>: the length of the window to use for 'end-trimming'. Default 5.\n";
    print
"-r quality score <optional>: windows with an average quality score less than this will be trimmed (0-34). Default 20.\n";
    print
"-a quality score <optional>: trimmed reads with an avg. quality of less than this will be throw out (0-34). Default 20.\n";
    print
"-p string: the path and file prefix for your output files eg. /home/blah/quality_trimmed_data\n";
    print "-h: this message\n";
    exit;
}

sub init() {
    my $opt_string = 'f:q:l:s:r:a:p:h';
    (@ARGV and getopts("$opt_string", \%opt)) or usage();
    usage() if ($opt{'h'} or not($opt{'f'} and $opt{'q'}));
}

MAIN:
{
    init();

    my $window_size = (defined $opt{'s'}) ? $opt{'s'} : $DEFAULT_WINDOW_SIZE;
    my $window_quality_threshold =
      (defined $opt{'r'}) ? $opt{'r'} : $DEFAULT_WINDOW_MIN_QUAL;
    my $min_avg_quality =
      (defined $opt{'a'}) ? $opt{'a'} : $DEFAULT_WINDOW_AVG_QUAL;
    my $min_length = (defined $opt{'l'}) ? $opt{'l'} : $DEFAULT_MIN_LENGTH;

    #statistics

    #reads a pair at a time
    my $reads = IO::File->new($opt{'f'}, "r");
    my $quals = IO::File->new($opt{'q'}, "r");

    my $prefix = $opt{'p'};

    my $outfile_trimmed_pair_reads =
      IO::File->new($prefix . "_trimmed_mate_pairs.fasta", "w");
    my $outfile_trimmed_pair_quals =
      IO::File->new($prefix . "_trimmed_mate_pairs.qual", "w");
    my $outfile_trimmed_solitary_reads =
      IO::File->new($prefix . "_trimmed_solitary_reads.fasta", "w");
    my $outfile_trimmed_solitary_quals =
      IO::File->new($prefix . "_trimmed_solitary_reads.qual", "w");

    my $finished_first = 0;
    my %id_removed;

    my %obj_first;
    my %obj_second;

    my %file_obj = ('reads' => $reads, 'quals' => $quals);
    my $done = 0;

    while (!$done) {

        #reads and trims the forward and reverse reads
        my $ret_val1 = trim(\%file_obj, $window_size, $window_quality_threshold,
            $min_avg_quality);
        my $ret_val2 = trim(\%file_obj, $window_size, $window_quality_threshold,
            $min_avg_quality);

#if the returned objects are undefined altogether, there are no more reads to process.
        if (!(defined $ret_val1 or defined $ret_val2)) {
            $done = 1;
            next;
        }

        my $length1 =
          (defined($ret_val1->{'trimmed_seq'}))
          ? length($ret_val1->{'trimmed_seq'})
          : 0;
        my $length2 =
          (defined($ret_val2->{'trimmed_seq'}))
          ? length($ret_val2->{'trimmed_seq'})
          : 0;

        if (!defined $ret_val1 or !defined $ret_val2) {
            warn "There was an error : uneven number of sequences\n";
        }

        if (    defined $ret_val1->{'trimmed_seq'}
            and defined $ret_val2->{'trimmed_seq'}
            and $length1 >= $min_length
            and $length2 >= $min_length)
        {
            print $outfile_trimmed_pair_reads ($ret_val1->{'id'}) . "\n"
              . ($ret_val1->{'trimmed_seq'}) . "\n";
            print $outfile_trimmed_pair_quals ($ret_val1->{'id'}) . "\n"
              . ($ret_val1->{'trimmed_qual'}) . "\n";
            print $outfile_trimmed_pair_reads ($ret_val2->{'id'}) . "\n"
              . ($ret_val2->{'trimmed_seq'}) . "\n";
            print $outfile_trimmed_pair_quals ($ret_val2->{'id'}) . "\n"
              . ($ret_val2->{'trimmed_qual'}) . "\n";
        }
        elsif (defined $ret_val2->{'trimmed_seq'} and $length2 >= $min_length) {
            print $outfile_trimmed_solitary_reads ($ret_val2->{'id'}) . "\n"
              . ($ret_val2->{'trimmed_seq'}) . "\n";
            print $outfile_trimmed_solitary_quals ($ret_val2->{'id'}) . "\n"
              . ($ret_val2->{'trimmed_qual'}) . "\n";
        }
        elsif (defined $ret_val1->{'trimmed_seq'} and $length1 >= $min_length) {
            print $outfile_trimmed_solitary_reads ($ret_val1->{'id'}) . "\n"
              . ($ret_val1->{'trimmed_seq'}) . "\n";
            print $outfile_trimmed_solitary_quals ($ret_val1->{'id'}) . "\n"
              . ($ret_val1->{'trimmed_qual'}) . "\n";
        }

        if ($length1 < $min_length) {
            $too_short++;
        }
        if ($length2 < $min_length) {
            $too_short++;
        }
    }

    print
"Completed:  $count_trimmed_to_zero sequences were end-trimmed all the way, $count_avg_qual_low sequences had too low an average quality and ";
    print "$count_contained_ns sequences contained N's\n";
    print
"And there were $too_short sequences that were less than the length threshold post trimming\n";
}

#trims the next read in the file
sub trim {
    my ($file_pair, $window_size, $window_quality_threshold, $min_avg_quality)
      = @_;

    my $reads_fh = $file_pair->{'reads'};
    my $quals_fh = $file_pair->{'quals'};

    my $id_line_qual  = <$quals_fh>;
    my $id_line_reads = <$reads_fh>;

    if (!defined $id_line_qual) {
        return undef;    #finished
    }

    if ($id_line_qual ne $id_line_reads) {
        die "Identifier lines do not match: $id_line_qual vs $id_line_reads \n";
    }

    chomp($id_line_qual);
    chomp($id_line_reads);

    my $quality_line = <$quals_fh>;
    chomp($quality_line);

    my $reads_line = <$reads_fh>;
    chomp($reads_line);

    #returned because read contained at least one N
    if ($reads_line =~ /N/gi) {
        $count_contained_ns++;
        return {
            'id'           => $id_line_reads,
            'trimmed_seq'  => undef,
            'trimmed_qual' => undef
        };
    }

    my @quality_array = split(m/\-/, $quality_line);
    my $trim_to_base = $#quality_array;

    my $min_base = 1;    #because the first position is biased.

    for (my $i = scalar(@quality_array) - 1 ; $i >= $min_base ; $i--) {
        if ($quality_array[$i] == $QUAL_UNOBSERVED_NUCLEOTIDE) {
            if (    $quality_array[$i] == $QUAL_UNOBSERVED_NUCLEOTIDE
                and $trim_to_base != $i)
            {
                warn
                  "Score of '2' seen internally within read: $id_line_reads\n"
                  ;      #we're assuming that B's are not internal!
            }
            $trim_to_base = $i - 1;    #inclusive.
        }

        # else no 'B's
        elsif ($i > ($window_size - 1)) {
            my $window_sum = 0;

            #double check that this is the correct window size....
            my $window_start = $i - $window_size + 1;
            my $window_end   = $i;

            foreach my $qual_in_window (
                @quality_array[ $window_start .. $window_end ])
            {
                $window_sum += $qual_in_window;
            }

            if (($window_sum / $window_size) < $window_quality_threshold) {
                $trim_to_base = $window_start - 1;
            }
            else {
                last;
            }
        }
    }

    if ($trim_to_base < $min_base
      ) #if we trim to base 0, we remove base 0 anyhow, as recommended by Torsten.
    {
        $count_trimmed_to_zero++;
        return {
            'id'           => $id_line_reads,
            'trimmed_seq'  => undef,
            'trimmed_qual' => undef
        };
    }

    if ($trim_to_base == (scalar(@quality_array) - 1)) {
        $trim_to_base--;    #the last base is weird according to Torsten
    }

    my $quality_sum = 0;

    for (my $i = $min_base ; $i <= $trim_to_base ; $i++) {
        $quality_sum += $quality_array[$i];
    }

    my $average = $quality_sum / ($trim_to_base - $min_base + 1);

    if ($average < $min_avg_quality) {
        $count_avg_qual_low++;
        return {
            'id'           => $id_line_reads,
            'trimmed_seq'  => undef,
            'trimmed_qual' => undef
        };
    }

    #only case when something is returned..
    my @split_seq    = split("", $reads_line);
    my @trimmed_seq  = @split_seq[ $min_base .. $trim_to_base ];
    my @trimmed_qual = @quality_array[ $min_base .. $trim_to_base ];

    my %new_obj = (
        'id'           => $id_line_reads,
        'trimmed_seq'  => join("", @trimmed_seq),
        'trimmed_qual' => join("-", @trimmed_qual)
    );

    return \%new_obj;
}

