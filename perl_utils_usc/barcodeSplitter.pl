#!/usr/bin/perl
use File::Basename;
#####SET Consts
$perfectQualityCutoff = 25;
$singleQualityErrorCutoff = 25;
@availableBarcodes = ("ACG","CAT","GTA","TGC");
####################

open($infh, $infile = $ARGV[0]) || die "couldnt open $infile";
while(!eof($infh))
{
        $line[0] = readline($infh); $line[1] = readline($infh); $line[2] = readline($infh); $line[3] = readline($infh);
        my $code = substr($line[1],0,length($availableBarcodes[0]));
        my $qual = substr($line[3],0,length($availableBarcodes[0]));
        my $resultFile = "unknown_barcodes";
        my $resultLines = $line[0] . $line[1] . $line[2] . $line[3];
        $barcodeStats{$code}++;
        for my $availableBarcode (@availableBarcodes)
        {
                if ($h = &hammingdist($code, $availableBarcode) <= 1)
                {
                        if(!($hasPerfectQuality{$qual}))
                        {
                                $hasPerfectQuality{$qual} += (ord($_) < ($perfectQualityCutoff + 64) ? 0 : 1) for split(//,$qual);
                                $hasSingleQualityError{$qual} += (ord($_) < ($singleQualityErrorCutoff + 64) ? 0 : 1) for split(//,$qual);
                        }
                        $resultFile = (($h == 0 && $hasPerfectQuality{$qual} == length($qual)) || ($h == 1 && $hasSingleQualityError{$qual} >=  length($qual)-1)) ? $availableBarcode : "lowQual_barcodes";
                        $resultLines = $line[0] . substr($line[1],5) . $line[2] . substr($line[3],5) unless $resultFile eq "lowQual_barcodes";
                }
        }
        open($resultFiles{$resultFile},">$resultFile" . "_" . basename($infile)) unless $resultFiles{$resultFile};
        my $fh =  $resultFiles{$resultFile} || die "couldn't write";
        print $fh $resultLines;
}
open($summary,">summary_" . basename($infile));
print $summary "$_ $barcodeStats{$_}\n" for (sort {$barcodeStats{$b} <=> $barcodeStats{$a} } keys %barcodeStats);
sub hammingdist{ length( $_[ 0 ] ) - ( ( $_[ 0 ] ^ $_[ 1 ] ) =~ tr[\0][\0] ) }
