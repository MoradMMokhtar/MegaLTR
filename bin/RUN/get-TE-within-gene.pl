# The perl script used to identify LTR-RTs gene chimera based on LTR-Retrotransposons localization in genome sequence. The Perl script compared the gene start and end with the start and end of LTR-RT within the genome. LTR-RT was considered a LTR-RT gene chimera if it was located within the gene start and end coordinates provided by the gene annotation in the GFF files.

use strict;
my $searchfile = $ARGV[0];
my $searchfor  = $ARGV[1];
my $processid  = $ARGV[2];
my $OUT = $ARGV[3];
my @parray;
open (OUT,">$OUT/LTR_inside_genes.table.2");
open (PFILE,"<$searchfile");
while (<PFILE>) {
    chomp();
    if (/(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+([\S|\s]+)/) {
        my $id1   = $1;
        my $id2   = $2;
        my $id3   = $3;
        my $id4   = $4;
        my $id5   = $5;
        my $id6   = $6;
        my $id7   = $7;
        my $id8   = $8;
         
        push( @parray, [ $id1, $id2, $id3, $id4, $id5, $id6, $id7, $id8] );
    }
}
close PFILE;
open (GFILE,"<$searchfor");

while (<GFILE>) {
    chomp();
    if (/(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+([\S|\s]+)/) {
        my $contig  = $1;
        my $edta  = $2;
        my $region    = $3;
        my $start    = $4;
        my $end    = $5;
        my $all    = $6;
        my $all2    = $7;
        my $all3    = $8;
        my $all4    = $9;

        
        for my $i ( 0 .. $#parray ) {

      if($contig eq $parray[$i][1]&&$start <= $parray[$i][2]&& $end >= $parray[$i][3])
         {          
       

         print "$parray[$i][0]\t$parray[$i][1]\t$parray[$i][2]\t$parray[$i][3]\t$parray[$i][4]\t$parray[$i][5]\t$parray[$i][6]\t$parray[$i][7]\tYes\t$region\t$start\t$end\t$all\t$all2\t$all3\t$all4\n";
         
        print OUT "0\t$parray[$i][0]\tLTR-gene chimera\tInside the gene\t$parray[$i][1]\t$parray[$i][2]\t$parray[$i][3]\t$parray[$i][4]\t$parray[$i][5]\t$parray[$i][6]\t$parray[$i][7]\t$region\t$start\t$end\t$all\t$all2\t$all4\n";   
         } 
        }

    }
}
close GFILE;
