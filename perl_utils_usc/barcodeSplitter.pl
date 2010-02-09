#!/usr/bin/perl
use File::Basename;
#####SET Consts
$threeQualMatch = 25;
$twoQualMatch = 25; 
@availableBarcodes = ("ACG","CAT","GTA","TGC");
####################

$threeQualMatch += 64;
$twoQualMatch += 64; 

$infile = $ARGV[0] || die "specify a file";
open($infh, "<$infile") || die "couldnt open $infile";
$infile = basename($infile);
while(!eof($infh))
{
	$line[0] = readline($infh); $line[1] = readline($infh); $line[2] = readline($infh); $line[3] = readline($infh);
	my $code = substr($line[1],0,3);
        my $qual = substr($line[3],0,3);
	my $failed = "unknown_barcodes";
	$barcodeStats{$code}++;
        for my $availableBarcode (@availableBarcodes)
        {
                if ($h = &hammingdist($code, $availableBarcode) <= 1)
                {
                        if(!($qualityThreeMatch{$qual}))
                        {
                                my ($a,$b,$c) = split(//,$qual);
                                $qualityThreeMatch{$qual} = ord($a) >= $threeQualMatch && ord($b) >= $threeQualMatch && ord($c) >= $threeQualMatch ? 1 : -1;
                                $qualityTwoMatch{$qual} = ord($a) >= $twoQualMatch && ord($b) >= $twoQualMatch || ord($a) >= $twoQualMatch && ord($c) >= $twoQualMatch || ord($b) >= $twoQualMatch && ord($c) >= $twoQualMatch ? 1 : -1;
                        }
                        if(($h == 0 && $qualityThreeMatch{$qual} == 1) || ($h == 1 && $qualityTwoMatch{$qual} == 1)) 
			{
				open($barcodes{$availableBarcode},">$availableBarcode" . "_$infile") if(!($barcodes{$availableBarcode}));
				my $fh =  $barcodes{$availableBarcode} || die "couldn't write";
				print $fh $line[0] . substr($line[1],5) . $line[2] . substr($line[3],5); 
				undef $failed;
			}
                        else 
			{
				$failed = "lowQual_barcodes";
			}
                }
        }
	if($failed)
	{
		open($barcodes{$failed},">$failed" . "_$infile") if(!($barcodes{$failed}));
		my $fh =  $barcodes{$failed} || die "couldn't write";
		print $fh $line[0] . $line[1] . $line[2] . $line[3]; 
	}
}
open($summary,">summary_$infile");
foreach $key (sort {$barcodeStats{$b} <=> $barcodeStats{$a} } keys %barcodeStats)
{
     print $summary "$key $barcodeStats{$key}\n";
}

sub hammingdist{ length( $_[ 0 ] ) - ( ( $_[ 0 ] ^ $_[ 1 ] ) =~ tr[\0][\0] ) }
