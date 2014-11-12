#!/uaopt/perl/5.14.2/bin/perl

$| = 1;

use strict;
use warnings;
use feature 'say';
use autodie;
use Cwd 'cwd';
use File::Basename qw(basename fileparse);
use File::Path;
use File::Spec::Functions;
use Getopt::Long;
use Pod::Usage;
use Bio::SeqIO;
use Bio::Seq::Quality;

my $out_dir     = cwd();;
my $verbose     = 0;
my ($help, $man_page);

GetOptions(
    'o|out=s'       => \$out_dir,
    'v|verbose'     => \$verbose,
    'help'          => \$help,
    'man'           => \$man_page,
) or pod2usage(2);

if ($help || $man_page) {
    pod2usage({
        -exitval => 0,
        -verbose => $man_page ? 2 : 1
    });
}

if (!@ARGV) {
    pod2usage('No input files');
}

if (!-d $out_dir) {
    mkpath $out_dir;
}

my $file_num = 0;
# file is in fastq format
for my $file (@ARGV) {
   open (F, "$file.fastq.trimmed.paired") || die "Cannot open fastq file\n";
   open (TEMP, ">$out_dir/$file.temp") || die "can't open temp\n";
   my $check_next = 0;
   LINE:
   while (<F>) {
       if ($_ =~ /^\+/) {
          $check_next++;
          print TEMP $_;
          next LINE;
       }
       if ($check_next == 1) {
          my $check = $_;
          if ($check =~ /^\@/) {
             $check =~ s/^\@(.*)/$1/;
             $check = "h" . $check;
          }
          print TEMP $check;
          $check_next = 0;
          next LINE;
       }
       print TEMP $_;
   }
   close F;
   close TEMP;
    
   my $out_seq_obj = Bio::SeqIO->new(
       -file   => ">$out_dir/$file.fa",
       -format => 'fasta',
   );
   
   my $out_qual_obj = Bio::SeqIO->new(
       -file   => ">$out_dir/$file.qual",
       -format => 'qual',
   );
   
   my $in_fastq_obj = Bio::SeqIO->new(
       -file   => "$out_dir/$file.temp",
       -format => 'fastq'
   );
   
   while (my $seq = $in_fastq_obj->next_seq) {
      $out_seq_obj->write_seq($seq);
      $out_qual_obj->write_seq($seq);
   }

   unlink "$out_dir/$file.temp";
} 

