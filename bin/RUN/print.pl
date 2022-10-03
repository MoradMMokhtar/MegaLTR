use strict;
my $input = $ARGV[0];
my $under = "_";
open( PFILE, "<$input" );
while (<PFILE>) {
    chomp();
    if (/(\S+)\t(\S+)\t(\S+)\t(\S+)\t([\S|\s]+)/) {
        my $id1   = $1;
        my $id2   = $2;
        my $id3   = $3;
        my $id4   = $4;
        my $id5   = $5;
        print "$id4$under$id1$under$id2\t$id4\t$id1\t$id2\t$id3\t$id5\n";
       }
         } 
close PFILE;

