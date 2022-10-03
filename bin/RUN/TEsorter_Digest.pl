use strict;
my $searchfile = $ARGV[0];
my $searchfor  = $ARGV[1];
my @parray;
open( PFILE, "<$searchfile" );
while (<PFILE>) {
    chomp();
    if (/(\S+)\t([\S|\s]+)/) {
        my $id   = $1;
        my $rnu   = $2;
        push( @parray, [ $id, $rnu] );
    }
}
close PFILE;
open( GFILE, "<$searchfor" );
while (<GFILE>) {
    chomp();
    if (/(\S+)\t([\S|\s]+)/) {
        my $contig  = $1;
        my $gstart  = $2;
        for my $i ( 0 .. $#parray ) 
         {          
       if($contig eq $parray[$i][0])    
                {   
         print "$parray[$i][0]\t$parray[$i][1]\t$gstart\n";
       }
         } 
        }
    }
close GFILE;
