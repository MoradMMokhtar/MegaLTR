#!/usr/bin/perl 
use strict;
use Bio::SeqIO;
my %ids;
my $idssandefile = $ARGV[1];
my $result = $ARGV[2];
open(IDS,"<$idssandefile");
while(<IDS>)
{
#primername Id start end
if(/(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)/)
{
   push(@{$ids{$1}},[$2,$3,$4,$5,$6,$7]);
}
}
close IDS;

foreach my $k(keys(%ids))
{
 foreach my $s(@{$ids{$k}})
 {
 #print $s->[1];
 }
}
getseqbystartend();

sub getseqbystartend
{
 
my $seqfile = Bio::SeqIO->new( -file => $ARGV[0], -format => "fasta" );
while ( my $seqobj = $seqfile->next_seq )
{

 if ($ids{$seqobj->display_id()})
   {
my @array=@{$ids{$seqobj->display_id()}};
foreach my $se(@array)
{
  my @searray=@{$se};
  my $st = $searray[0];
  my $et = $searray[1];
  my $pm = $searray[2];  
  my $sm = $searray[3];
  my $sm1 = $searray[4];
  my $sm2 = $searray[5];
  my $sm3 = $searray[6];

 my $seqout=substr ($seqobj->seq(),$st,($et-$st));
 my $seqoutseq = Bio::Seq->new(-seq=>$seqout, -format => 'fasta');
  $seqoutseq->display_id(  $seqobj->display_id() ."..$st..$et|$pm|$sm|$sm1|$sm2|$sm3" );
my $writefile = Bio::SeqIO->new(-id=>"$result",-file => ">>$result", -format => "fasta");
$writefile->write_seq($seqoutseq);
}

}
}
}

