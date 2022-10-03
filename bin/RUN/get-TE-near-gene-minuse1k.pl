# get LTR position in gene.
use strict;
my $searchfile = $ARGV[0];
my $searchfor  = $ARGV[1];
my $up  = $ARGV[2];
my $down  = $ARGV[3];
my @parray;
open (PFILE,"<$searchfile");
while (<PFILE>) {
    chomp();
    if (/(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+([\S|\s]+)/) {
        my $id   = $1;
        my $rnu   = $2;
        my $rnu1   = $3;
        my $rnu2   = $4;
        my $rnu3   = $5;
        my $rnu7   = $3-$up;
        my $rnu8   = $4+$down;
        push( @parray, [ $id, $rnu, $rnu1, $rnu2, $rnu3, $rnu7, $rnu8] );
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
        my $all1    = $7;
        my $all2    = $8;
        my $all3    = $9;

#        my $before  = $5+1000;
#        my $plus    = '-';
        
        for my $i ( 0 .. $#parray ) {

      if($contig eq $parray[$i][1]&& $parray[$i][5] <= $start && $parray[$i][2] >=  $end)
         {          
       
my $length  = $parray[$i][2]-$end;

         print "$parray[$i][0]\tGenes Upstream LTR\t$parray[$i][1]\t$parray[$i][2]\t$parray[$i][3]\t$parray[$i][4]\t$region\t$start\t$end\t$length\t$all1\t$all3\n";
       
         } 
        }

    }
}
close GFILE;
