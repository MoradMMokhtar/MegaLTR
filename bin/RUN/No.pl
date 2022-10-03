use strict;
my $input = $ARGV[0];
my $processid = $ARGV[1];
my @parray;
open( PFILE, "<$input" );
while (<PFILE>) {
    chomp();
    if (/([\S|\s]+)/) {
        my $id1   = $1;
        
        print "$id1\tNo\n";
       }
         } 
close GFILE;
