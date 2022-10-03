my $settime = $ARGV[1];
$parameter_={};
$parameter_->{"rateOfEvolution"}="arabidopsis";

#print"$settime \n";
$s=align2K("",$ARGV[0], $parameter_, "");
print($s);

sub align2K {

my ($pairID, $alignFile, $parameter_, $MESSAGES) = @_;

# calculate distances for all pairs of LTRs

	# store all the lines in alignment file
	my @lines = ();

	open(IN, $alignFile) || die "\ncan't open alignment file $alignFile: $! \n\n********* EXITING NOW *********\n";
	push (@lines, $_) while ( defined($_ = <IN>) );
	close(IN) || print $MESSAGES "can't close file: $! \n";
	
        my $clustalwOFFSET = 0;	  # will store the number of lines in the header of clustalw .aln files (= no.lines before alignment proper)
        $clustalwOFFSET++ while ( $clustalwOFFSET <= scalar(@lines) && $lines[$clustalwOFFSET] !~ /(5'|3')/ );

	# store lines corresponding to the LTR1 and LTR2 sequences
	my @LTR1lines = grep /5'/ , @lines;
	my @LTR2lines = grep /3'/ , @lines;

	# remove 'newline' characters and sequence identifier to leave only the nucleotide sequences on each "line"
	grep { chomp; s/^\S+\s+([\w-]+)$/$1/ } @LTR1lines;
	grep { chomp; s/^\S+\s+([\w-]+)$/$1/ } @LTR2lines;

	# store the lines corresponding to site identity/difference (asterisks/spaces)
	my @consensusLines;
	for (my $lineNum=1; $lineNum <= scalar(@LTR1lines); $lineNum++) {
		push(@consensusLines, $lines[($clustalwOFFSET-2) + 4 * $lineNum]);  # offset @lines to avoid the clustalW header, then get every fourth line
		# remove 'newline' and keep only the characters that refer to the LTR sequences
		chomp($consensusLines[$lineNum-1]);
		my $numSitesOnThisLine = length($LTR1lines[$lineNum-1]);
		$consensusLines[$lineNum-1] =~ s/.*(.{$numSitesOnThisLine})$/$1/;
	}

	# Concatenate LTR and consensus lines into single strings
	my $LTR1 = join '', @LTR1lines;
	my $LTR2 = join '', @LTR2lines;
	my $consensus = join '', @consensusLines;
	
	# get distances (Kimura and Tajima&Nei, plus standard deviations) and other statistics (
	# proportion of transitions and transversions, number of sites compared between the 2 seqs,
	# number of indels, transition/transversion ratio).
	return getDistances($LTR1, $LTR2, $consensus, $parameter_);

}

sub getDistances {
	# First two args are the sequences (previously made the same length by a clustalW alignment (extra hifens added if necessary),
	# third the consensus line (clustalW format). Fourth is the rate of synonymous nucleotide substitutions per site per million years.
	# Last argument is the output file handle (key in parameter hash).
	my($seq1, $seq2, $consensus, $parameter_) = @_;
	
	my $rateOfEvolution = $parameter_->{"rateOfEvolution"};
	my $arabidopsisRateOfEvolution = $settime;  # number of synonymous nucleotide substitutions per site per million years
# my $humanLowRateOfEvolution = 1.1 * 10**(-3);   	#low HUMAN number of synonymous nucleotide substitutions per site per million years
# my $humanHighRateOfEvolution = 2.1 * 10**(-3);   	#high HUMAN number of synonymous nucleotide substitutions per site per million years
	$rateOfEvolution = $arabidopsisRateOfEvolution if ($parameter_->{"arabidopsis"});
# translate sequences into numbers: A=0, C=1, G=2, T=3 (keep gaps (hyfens) as they are)
#	$seq1 =~ tr/ACGTacgt/01230123/;
#	$seq2 =~ tr/ACGTacgt/01230123/;

	# split sequences and consensus into arrays of single characters
	my @seq1 = split(//, $seq1);
	my @seq2 = split(//, $seq2);
	my @consensus = split(//, $consensus);

	my $numComparedSites = 0;
	my $transitions = 0;
	my $transversions = 0;
	my $indels = 0;
	my @ACGT_content = (0,0,0,0);
	my %ACGT_pairs = ("AC", 0, "AG", 0, "AT", 0, "CG", 0, "CT", 0, "GT", 0);

	# find site corresponding to the beginning of the alignment in case there is an initial gap
	my $initialSite = 0;
	$initialSite++ while ( $seq1[$initialSite] eq "-" || $seq2[$initialSite] eq "-" );
	# find site corresponding to the end of the alignment in case there is a final gap
	my $finalSite = length($consensus) - 1;
	$finalSite-- while ( $seq1[$finalSite] eq "-" || $seq2[$finalSite] eq "-" );

	$numComparedSites = $finalSite - $initialSite + 1;

	# get differences between sequences
	my $basesAtThisSite;
	my $lastIndelSite = $initialSite - 1;
	for (my $site = $initialSite; $site <= $finalSite; $site++) {
		$basesAtThisSite = $seq1[$site].$seq2[$site];

		if ( $consensus[$site] eq " " ) {  # this site is different between the two sequences

			if ( $basesAtThisSite =~ /-/ ) {  			# indel!
				$indels++ if ( $lastIndelSite < $site - 2 );  # do not count 2 indels if they're separated by only one base match
				$lastIndelSite = $site;
				$numComparedSites--;  # do not include gap in site comparisons
			}
			elsif ( $basesAtThisSite =~ /AC|CA/i ) {  	# A and C at this site
				$ACGT_content[0]++; $ACGT_content[1]++;
				$ACGT_pairs{"AC"}++;
				$transversions++;
			}
			elsif ( $basesAtThisSite =~ /AG|GA/i ) {  	# A and G at this site
				$ACGT_content[0]++; $ACGT_content[2]++;
				$ACGT_pairs{"AG"}++;
				$transitions++;
			}
			elsif ( $basesAtThisSite =~ /AT|TA/i ) {  	# A and C at this site
				$ACGT_content[0]++; $ACGT_content[3]++;
				$ACGT_pairs{"AT"}++;
				$transversions++;
			}
			elsif ( $basesAtThisSite =~ /CG|GC/i ) {  	# C and G at this site
				$ACGT_content[1]++; $ACGT_content[2]++;
				$ACGT_pairs{"CG"}++;
				$transversions++;
			}
			elsif ( $basesAtThisSite =~ /CT|TC/i ) {  	# C and T at this site
				$ACGT_content[1]++; $ACGT_content[3]++;
				$ACGT_pairs{"CT"}++;
				$transitions++;
			}
			elsif ( $basesAtThisSite =~ /GT|TG/i ) {  	# G and T at this site
				$ACGT_content[2]++; $ACGT_content[3]++;
				$ACGT_pairs{"GT"}++;
				$transversions++;
			}
		}

		else {  # this site is identical between the two sequences
			$ACGT_content[0]+=2 if ( $basesAtThisSite =~ /A/i );  	# A and A at this site
			$ACGT_content[1]+=2 if ( $basesAtThisSite =~ /C/i );  	# C and C at this site
			$ACGT_content[2]+=2 if ( $basesAtThisSite =~ /G/i );  	# G and G at this site
			$ACGT_content[3]+=2 if ( $basesAtThisSite =~ /T/i );  	# T and T at this site
		}
	}

	# convert numbers of nucleotides into frquencies
	foreach my $nucleotideNum (@ACGT_content) {
		$nucleotideNum /= 2*$numComparedSites;
	}
	# convert numbers of mis-matched nucleotides into frequencies
	my @ACGT_pairs;
	foreach my $key (sort keys %ACGT_pairs) {
		push(@ACGT_pairs, $ACGT_pairs{$key}/$numComparedSites);
	}
	# get K (and its variance) using Kimura's 2-parameter method
	my($K_Kimura, $V_Kimura) = getK_Kimura($numComparedSites, $transitions, $transversions);

	# get K (and its variance) using Tajima and Nei's method
	my($K_TajimaNei, $V_TajimaNei) = getK_TajimaNei($numComparedSites, $transitions+$transversions,
										@ACGT_pairs, @ACGT_content);

	my $K_Ksd = sqrt($V_Kimura);  # standard deviation
	my $K_TNsd = sqrt($V_TajimaNei);

	# time of divergence betwen two seqs in millions years
	my $timeK = "NA";
	my $timeTN = "NA";
	my $timeKsd = "NA";
	my $rateOfEvolution = $settime;
	if ($rateOfEvolution) {
	  $timeTN = $K_TajimaNei/(2*$rateOfEvolution);
	  $timeK = $K_Kimura/(2*$rateOfEvolution);
	  # compute variance of time estimate including variance for the Poisson process of nucleotide substitutions
	  $timeKsd = sqrt( ($V_Kimura*$numComparedSites**2 + $K_Kimura*$numComparedSites)/(2*$rateOfEvolution*$numComparedSites)**2 );
	}

	$transitions_numComparedSites=($transitions/$numComparedSites);
	$transversions_numComparedSites=($transversions/$numComparedSites);
	return "$K_Kimura\t$K_Ksd\t$K_TajimaNei\t$K_TNsd\t$timeK\t$timeKsd\t $numComparedSites\t$transitions\t$numComparedSites\t$transversions\t$numComparedSites\t$timeTN\t$transitions_numComparedSites\t$transversions_numComparedSites";
	# return field values
	#return $K_Kimura, $K_Ksd, $K_TajimaNei, $K_TNsd, $timeK, $timeKsd, $numComparedSites,
	#       $transitions/$numComparedSites, $transversions/$numComparedSites, $timeTN;

}

sub getK_Kimura {
	# Last two args are the number of transitions and transversions between two sequences. First argument
	# is the number of sites compared between the two sequences.
	my($length, $transitions, $transversions) = @_;
	$transitions /= $length;		# the proportion of transitions
	$transversions /= $length;	# the proportion of transversions

	my $a = 1/(1 - 2*$transitions -$transversions);
	my $b = 1/(1 - 2*$transversions);
	my $c = ($a + $b)/2;

	my $K = log($a)/2 + log($b)/4;  	# K distance

	# approximate sampling variance
	my $V = ($transitions*$a**2 + $transversions*$c**2 - ($a*$transitions + $c*$transversions)**2)/$length;

	return ($K, $V);
}

sub getK_TajimaNei {
	# First argument is the number of sites compared between the two sequences.
	# Second is the number if different nucleotides when the two seqs are aligned.
	# Next six args are the frequencies of non-matching pairs of nucleotides in the order:
	# "AC", "AG", "AT", "CG", "CT", and "GT".
	# Last arg is a list containing the frequencies of each nucleotide in the aligned seqs.
	my($length, $diffs, $fAC, $fAG, $fAT, $fCG, $fCT, $fGT, @ACGT_f) = @_;
	$diffs /= $length; 	# the proportion of different sites

	my @ACGT_pairs = ([0,$fAC,$fAG,$fAT],[0,0,$fCG,$fCT],[0,0,0,$fGT]);

	my $b1 = 1 - $ACGT_f[0]**2 - $ACGT_f[1]**2 - $ACGT_f[2]**2 - $ACGT_f[3]**2;
	my $h = 0;
	for (my $i = 0; $i < 3; $i++) {
		for (my $j= $i + 1; $j <= 3; $j++) {
			$h += ($ACGT_f[$i] && $ACGT_f[$j])? $ACGT_pairs[$i][$j]**2/(2*$ACGT_f[$i]*$ACGT_f[$j]) : 0;
		}
	}
	my $b = ($h) ? ( $b1 + $diffs**2/$h )/2 : 1;  # K will be zero if no diffs between seqs

	# K distance
	my $K = (-1) * $b * log(1 - $diffs/$b);

	# approximate sampling variance
	my $V = $b**2 * $diffs * (1-$diffs)/( ($b-$diffs)**2 * $length );

	return ($K, $V);
}
