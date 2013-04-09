#!/usr/bin/env perl

use strict;
my $samtools = "/home/uec-00/shared/production/software/samtools/samtools";
my $PRINTLANEIDS = 1;
my $DOFASTQ = 1;
my $DOREPEAT = 1;
my $DOALIGNEDCOUNTS = 1;
my $DODEPTHWINDOWS = 1;
my $DOCONVERSION = 1;
my $DOCONTAMALIGNTEST = 1;

my $INTERMEDIATE_DIRS = 0;

# State vars
my $headersPrinted = 0;

foreach my $dir (@ARGV)
{
    print STDERR "Searching DIR $dir\n";
	
    foreach my $laneNum (1..8)
    {
	my @dirlist = glob("$dir/s_$laneNum*");
	@dirlist = glob("$dir/*L00$laneNum*") if !@dirlist;
        if($dirlist[0])
        {
        
	        my @flds = ();
	        my @headers = ();
			
	        # GENERAL
	        if ($PRINTLANEIDS)
	        {
	            push(@flds,$dir); push(@headers,"FlowCelln");
	            push(@flds,$laneNum); push(@headers,"laneNum");
	        }
	
	        my $prefix = ($INTERMEDIATE_DIRS) ? "/*s_${laneNum}*/*" : "/*s_${laneNum}*";
	        my @tmpCnt = glob($prefix);
	        $prefix = "*L00$laneNum*" if !@tmpCnt;
	
	        # FASTQ counts
	        if ($DOFASTQ)
	        {
	            my $nocontamN = seqCountFastqFiles($dir."${prefix}.nocontam.*");
	            push(@flds,$nocontamN); push(@headers,"nocontamSeqs");
	            my $contamN = seqCountFastqFiles($dir."${prefix}.contam.*");
	            push(@flds,$contamN); push(@headers,"contamSeqs");
	            my $contamPolyaN = seqCountFastqFiles($dir."${prefix}.contam.polya.*");
	            push(@flds,$contamPolyaN); push(@headers,"contamPolyaSeqs");
	            my $contamAdaptersN = seqCountFastqFiles($dir."${prefix}.contam.adapters.*");
	            push(@flds,$contamAdaptersN); push(@headers,"contamAdaptersSeqs");
	            my $contamAdapterTrimN = seqCountFastqFiles($dir."${prefix}.contam.adapterTrim.*");
	            push(@flds,$contamAdapterTrimN); push(@headers,"contamAdapterTrimSeqs");            
	        }
	
	        # Repeat counts
	        if ($DOREPEAT)
	        {
	            my $gaatgN = patternCountFiles($dir."${prefix}.nocontam.*", "GAATGGAATG");
	            push(@flds,$gaatgN); push(@headers,"GAATGGAATG");
	            my $tatttN = patternCountFiles($dir."${prefix}.nocontam.*", "TATTTTATTT");
	            push(@flds,$tatttN); push(@headers,"TATTTTATTT");
	            my $cattcN = patternCountFiles($dir."${prefix}.nocontam.*", "CATTCCATTC");
	            push(@flds,$cattcN); push(@headers,"CATTCCATTC");
	        }
	
	        if ($DOALIGNEDCOUNTS)
	        {
	            #my ($fullReads, $sampleReads, $dups) = alignedCounts($dir."/ResultCount_*_${laneNum}*map.q30.txt",$dir."/ReadCounts_*_${laneNum}_maq.csv");
	            #my ($fullReads, $sampleReads, $dups) = alignedCounts($dir."/ResultCount_*_${laneNum}.bam",$dir."/ReadCounts_*_${laneNum}_maq.csv");
	            #push(@flds, $fullReads); push(@headers,"AlignedReads");
	            #push(@flds, $sampleReads); push(@headers,"SampledAlignedReads");
	            #push(@flds, $dups); push(@headers,"AlignedDuplicateReads");
	            #push(@flds, ($sampleReads>0) ? ($dups/$sampleReads) : 0); push(@headers,"AlignedDuplicateFraction");

	            if(glob($dir."/ResultCount_*_${laneNum}.bam"))
	            {
			print STDERR "checking aln...\n";
	            	my ($fullReads, $sampleReads, $dups) = alignedCounts($dir."/ResultCount_*_${laneNum}.bam",$dir."/ReadCounts_*_${laneNum}_maq.csv");
	            	push(@flds, $fullReads); push(@headers,"AlignedReads");
	            	push(@flds, $sampleReads); push(@headers,"SampledAlignedReads");
	            	push(@flds, $dups); push(@headers,"AlignedDuplicateReads");
	            	push(@flds, ($sampleReads>0) ? ($dups/$sampleReads) : 0); push(@headers,"AlignedDuplicateFraction");
	            }
	            elsif(glob($dir."/s_${laneNum}_*tophat_hits.bam"))
	            {
			print STDERR "checking tophat aln...\n";
	            	my ($fullReads, $sampleReads, $dups) = alignedCounts($dir."/s_${laneNum}_*tophat_hits.bam",$dir."/ReadCounts_*_${laneNum}_maq.csv");
	            	push(@flds, $fullReads); push(@headers,"AlignedReads");
	            	push(@flds, $sampleReads); push(@headers,"SampledAlignedReads");
	            	push(@flds, $dups); push(@headers,"AlignedDuplicateReads");
	            	push(@flds, ($sampleReads>0) ? ($dups/$sampleReads) : 0); push(@headers,"AlignedDuplicateFraction");
	            }	            
	            elsif(glob($dir."/ResultCount_*_${laneNum}_*.bam"))
	            {
			print STDERR "checking aln...\n";
	            	my ($fullReads, $sampleReads, $dups) = alignedCounts($dir."/ResultCount_*_${laneNum}_*.bam",$dir."/ReadCounts_*_${laneNum}_maq.csv");
	            	push(@flds, $fullReads); push(@headers,"AlignedReads");
	            	push(@flds, $sampleReads); push(@headers,"SampledAlignedReads");
	            	push(@flds, $dups); push(@headers,"AlignedDuplicateReads");
	            	push(@flds, ($sampleReads>0) ? ($dups/$sampleReads) : 0); push(@headers,"AlignedDuplicateFraction");
	            }	            
	        }
	
	        if ($DODEPTHWINDOWS)
	        {
	            my ($windsFw, $windsRev, $windsTotal,$chroms) = depthWindows($dir."/ReadDepths_maxIden0_*_${laneNum}_*wind100000.csv");
	            my $numWinds = scalar(@$windsFw);
	#           push(@flds,@$windsFw); push(@headers,1..$numWinds);
	#           push(@flds,@$windsRev); push(@headers,1..$numWinds);
	            print STDERR "Num headers: " . scalar(@$chroms) . "\n";
	            push(@flds,@$windsTotal); push(@headers,@$chroms);
	
	            ($windsFw, $windsRev, $windsTotal,$chroms) = depthWindows($dir."/ReadDepths_maxIden1_*_${laneNum}_*wind100000.csv");
	            $numWinds = scalar(@$windsFw);
	#           push(@flds,@$windsFw); push(@headers,1..$numWinds);
	#           push(@flds,@$windsRev); push(@headers,1..$numWinds);
	            push(@flds,@$windsTotal); push(@headers,@$chroms);
	        }

		if ($DOCONVERSION)
	        { 
	            if(glob($dir."/*${laneNum}.pileup_cg_dinucleotide.csv"))
	            {
	            	my ($convCpG) = conversionFrac($dir."/*${laneNum}.pileup_cg_dinucleotide.csv");
	            	push(@flds, $convCpG); push(@headers,"ConversionCpG");
	            }
	            if(glob($dir."/*${laneNum}.pileup_cg_dinucleotide.csv"))
	            {
	            	my ($convCpH) = conversionFrac($dir."/*${laneNum}.pileup_ch_dinucleotide.csv");
	            	push(@flds, $convCpH); push(@headers,"ConversionCpH");
	            }
	        }	
	
            if($DOCONTAMALIGNTEST)
             {
                        if(glob("$dir/aligntest_s_$laneNum*"))
                        {
                                my @testAligns = glob("$dir/aligntest_s_$laneNum*");
                                my $chunkSize = 4;
                                if(-e "$dir/s_$laneNum\_1_*sequence.1.txt")
                                {
                                	$chunkSize = `wc -l $dir/s_$laneNum\_1_*sequence.1.txt`;
                                	$chunkSize =~ /^(\d+)\s/;
                                	$chunkSize = $1;	
                                }
                                elsif(-e "$dir/s_$laneNum\_*sequence.1.txt")
                                {
                                	$chunkSize = `wc -l $dir/s_$laneNum\_*sequence.1.txt`;
                                	$chunkSize =~ /^(\d+)\s/;
                                	$chunkSize = $1;	
                                }
                                
                                $chunkSize = $chunkSize / 4;
                                foreach my $testAlign (@testAligns)
                                {
                                        $testAlign =~ /aligntest_s_$laneNum\_(.+)\.map/;
                                        my $genome = $1;
                                        my $numberAligned = `maq mapview $testAlign | wc -l`;
                                        chomp $numberAligned;
                                        my $ratioAligned = $numberAligned / $chunkSize;
                                        push(@flds, $ratioAligned); push(@headers,"$genome" . "_q0");
                                        #print "$genome: $testAlign\t$numberAligned/$chunkSize = $ratioAligned\n";
                                }
                        }
              }
	
	        # Now print
	        if (!$headersPrinted && (grep {$_} @headers))
	        {
	            print join(",",@headers)."\n";
	            $headersPrinted = 1;
	        }
	
	        # First two are ids
	        if (scalar(grep {$_} @flds) > 2)
	        {
	            print join(",",@flds)."\n" ;
	        }
	
	    }
    }
}



sub seqCountFastqFiles
{
    my ($globPat) = @_;
    return (lineCountFiles($globPat)/4);
}

sub lineCountFiles
{
    my ($globPat) = @_;

    my @files = glob($globPat);
    my $total = 0;
    foreach my $f (@files)
    {
        $total += `wc -l $f`;
        print STDERR "\t\twc -l $f\ttotal=$total\n";
    }
    
    return $total;
}

sub patternCountFiles
{
    my ($globPat, $pat) = @_;

    my @files = glob($globPat);
    my $total = 0;
    foreach my $f (@files)
    {
        my $cmd = "grep -i \"$pat\" $f | wc -l";
        print STDERR "\t\t$cmd\ttotal=$total\n";
        $total += `$cmd`;
    }
    
    return $total;
}

sub alignedCounts
{
    my ($alignmentGlobPath, $depthsGlobPat) = @_;

    my @files = glob($alignmentGlobPath);
     
    my $alignedReads = 0;
    foreach my $f (@files)
    {
        my $samtoolsOutput = `$samtools view -q 30 $f | wc -l` ; chomp $samtoolsOutput;
	print STDERR "counting $f ...\n";
        $alignedReads += $samtoolsOutput;
#               die "Can't read file $f\n" unless open(F,$f);
#               while (my $line=<F>)
#               {
#                   $alignedReads++;
#               }
#               close(F);
    }

    @files = glob($depthsGlobPat);
    my $totalReads = 0;
    my $totalDups = 0;
    foreach my $f (@files)
    {
        die "Can't read file $f\n" unless open(F,$f);

        while (my $line=<F>)
        {
            chomp $line;
            my ($dummy, $dups, $count) = split(",",$line);
            $dups = (-1 * $dups) if ($dups < 0); # They can be negative for rev strand
            
            $totalReads += $dups * $count;
            $totalDups += (($dups - 1) * $count) if ($dups >= 2);
        }

        close(F);
        print STDERR "\t\t$f\ttotalReads=$totalReads\ttotalDups=$totalDups\n";
    }
    
    return ($alignedReads, $totalReads, $totalDups);
}

# Returns ([FWwinds],[REVwinds],[TOTALwinds],[chromNums])
sub depthWindows
{
    my ($globPat) = @_;

    my @files = glob($globPat);
    my @windCounts = ([],[],[],[]); # $windCounts[0]=FW, $windCounts[1]=REV, $windCounts[2]=TOTAL, CHROM_NUMS
    foreach my $f (@files)
    {
        die "Can't read file $f\n" unless open(F,$f);

        my @windsSeen = (0,0);
        my $onLine = 0;
        while (my $line=<F>)
#       while (($onLine<2000) && (my $line=<F>))
        {
            chomp $line;
            my @flds = split(",",$line);
            my $strandInd = ($flds[5] == -1) ? 1 : 0;
            my $count = $flds[8];

            my $l = $windCounts[$strandInd];
            @{$l}[$windsSeen[$strandInd]] += $count;
            $l = $windCounts[2]; # Total
            @{$l}[$windsSeen[$strandInd]] += $count;

            # Chrom number
            my $chr = chromToNum($flds[4]);
            $chr += ($flds[6] / (1E9)); # The coords
            $l = $windCounts[3];
            @{$l}[$windsSeen[$strandInd]] = $chr;

            $windsSeen[$strandInd]++;
            $onLine++;
        }

        close(F);
    }
    
    return @windCounts;
}

sub chromToNum
{
    my ($chr) = @_;

    $chr =~ s/^c(hr)?//;

    if ($chr =~ /^\d+$/)
    {
    }
    elsif ($chr eq 'x')
    {
        $chr = 25;
    }
    elsif ($chr eq 'y')
    {
        $chr = 26;
    }
    elsif ($chr eq 'm')
    {
        $chr = 27;
    }

    return $chr;
}

sub conversionFrac
{
    my ($globPat) = @_;

    my @files = glob($globPat);
#   print STDERR "($globPat) Files: ".join(", ",@files)."\n";
    my $totals = {};
    foreach my $f (@files)
    {
        die "Can't read file $f\n" unless open(F,$f);

        while (my $line=<F>)
        {
            chomp $line;
            my ($typeId, $nuc, $strand, $cycle, $quality, $count) = split(",",$line);

            $totals->{lc($nuc)} += $count;
        }
        close(F);

        my $totalCT = ( $totals->{t} + $totals->{c} );
        my $conv = ($totalCT>0) ? ($totals->{t} / $totalCT) : 0;

        print STDERR "\t\t$f\ttotalT=$totals->{t}\ttotalC=$totals->{c}\tconv=$conv\n";
        return $conv;
    }
}
