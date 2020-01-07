#!/usr/bin/env perl



#### CMMR unified microbiome pipeline ####
#### AKA PUMA ####




use warnings;
use strict;

#use DataBrowser qw(browseErr browse);
use Getopt::Long;
use Time::HiRes qw(time);
use SamCollection;
use TaxaObject;
use Set::IntervalTree;
use JSON::XS;
use Switch;
use Hash::Merge qw(merge);
use KeggObject;

my $nodeFile;
my $idCut;
my $sampleName;
my $nameFile;
my $samStr;
my $report;
my $keggFile;
my $isSorted;
my $prebuilt;
my $checkSpecies;
my $checkRest;
my $checkOutSpecies;
my $checkOutRest;
my $modules;
my $output;
my $tempDir = $ENV{'TMPDIR'};

unless($tempDir) {
	#print STDERR "tempDir = $tempDir\n";
	#die "TMPDIR not defined in ENV\n";
}

GetOptions ("samFiles=s" => \$samStr,
"nodes=s" => \$nodeFile,
"names=s" => \$nameFile,
"keggGenes=s" => \$keggFile,
"samplename=s" => \$sampleName,
"cutoff=s" => \$idCut,
"report=s" => \$report,
"isSorted" => \$isSorted,
"prebuilt=s" => \$prebuilt,
"checkpointInSpecies=s" => \$checkSpecies,
"checkpointInRest=s" => \$checkRest,
"checkpointOutSpecies=s" => \$checkOutSpecies,
"checkpointOutRest=s" => \$checkOutRest,
"modules=s" => \$modules,
"output=s" => \$output
);

my @chars = ("A".."Z", "a".."z");
my $string;
$string .= $chars[rand @chars] for 1..10;

sub parseFileStr {

	### Parses file string into list of categories and files.
	### Example = "virus=file1.sam,file2.sam,file3.sam,human=file4.sam,bacteria=file5.sam"

	my $fileStr = $_[0];
	my $reportStr = $_[1];
	my @fileStr = split ",", $fileStr;
	my @report = split ",", $reportStr;
	my $files = {};
	my $key = "";
	foreach my $part (@fileStr) {
		my $name = "";
		if ($part =~ m/(.*)=([^,]*)/) {
			$key = $1;
			$name = $2;
		} else {
			$name = $part;
		}
		unless ($key) {
			die "$fileStr not formatted correctly\n Ex: virus=file1.sam,file2.sam,human=file3.sam,file4.sam....\n";
		}
		if ($name) {
			push @{$files->{$key}}, $name;
		}
	}
	foreach my $report (@report) {
		unless (exists $files->{$report}) {
			die "report category $report defined, but no files exist with the same category\n";
		}
	}
	return($files);
}

sub revComp {
	
	#quick reverse complement

	my $seq = $_[0];
	$seq = uc($seq);
	$seq =~ tr/ACGTUYRSWKMBDHVN/TGCAARYSWMKVHDN/;
	$seq = reverse($seq);
	return ($seq);
}

sub parseLines {


	### parses sam lines from all categories, return a container of mismatches per category per reference
	
	
	my $compare = $_[0];
	my $sam = $_[1];
	my $taxa = $_[2];
	
	my $mismatchContainer= {};
	
	foreach my $category (@{$compare}) {
		my $reads = {};
		my $used = 0;
		my $checkRead = "";
		my $sequences = {};
		my $references = {};
		my $index = 0;
		my $pairs = {};
		foreach my $entry (@{$sam->{'currentLines'}->{$category}}) {
			$entry =~ m/$entryQR/;
			my $readName = $1;
			unless ($read) {
				$read = $readName;
			}
			my $code = $2;
			my $reference = $3;
			my $pos = $4;
			my $seq = $5;
			my $mismatch = $6;
			my $revSeq = revComp($seq);
			if (length($seq) < 100) {
				return 4;
			}
			if ($reference =~ m/$referenceQR/) {
				$reference = $1;
			}
			if ($report eq $category and not $taxa->getTaxaInfo($reference)) {
				next;
			}
			if (exists $sequences->{$revSeq}) {
				$seq = $revSeq;
			} elsif (not exists $sequences->{$seq}) {
				$sequences->{$seq} = 1;
			}
			my $end = $pos + length($seq);
			$references->{$category}->{$reference}->{$seq}->{$mismatch}->{$pos} = $end;
			#if ($mismatch < $mismatchMin) {
			#	$mismatchMin = $mismatch;
			#}
			unless ($used) {
				$checkRead = $readName;
				$used = 1;
			} elsif ($used and ($checkRead ne $readName)) {
				die "wrong entries in $category: $checkRead vs $readName\n";
			}
			#my $end = $pos + length($seq);
			#$mismatchContainer->{$mismatch}->{$category}->{$reference}->{$pos} = $end;
		}
		$mismatchContainer = cleanCategories($references);
	}
	return($mismatchContainer);
}

sub getBestTaxaFromBestCategory  {

	### returns from mismatchContainer the best species, or best resolved taxa if not resolved at the species level. 
	
	my $mismatchContainer = $_[0];
	my $bestCat = $_[1];
	my $taxa = $_[2];
	
	my $highestTaxa = {};
	my $spBucket = {};
	my $highestRank = 0;

	foreach my $reference (sort keys %{$mismatchContainer->{$mismatchMin}->{$bestCat}}) {
		my $refTaxID = $taxa->getTaxaId($reference);
		unless ($refTaxID) {
			next;
		}
		if ($reference =~ m/$seenNC/) {
			$seenNC = 1;
		}
		my $refRank = $taxa->highestRank($reference);
		if ($refRank > $highestRank) {
			$highestRank = $refRank;
			$highestTaxa = {};
			$spBucket = {};
			if ($taxa->{'names'}->{$refTaxID} =~ m/ sp. /) {
				$spBucket->{$refTaxID}->{$reference} = 1;
			} else {
				$highestTaxa->{$refTaxID}->{$reference} = 1;
			}
		} elsif ($refRank == $highestRank) {
			if ($taxa->{'names'}->{$refTaxID} =~ m/ sp. /) {
				$spBucket->{$refTaxID}->{$reference} = 1;
			} else {
				$highestTaxa->{$refTaxID}->{$reference} = 1;
			}
		}
	}
	return ($highestTaxa, $spBucket, $highestRank);
	
}

sub ecoliCheck {
	
	my $highestTaxa = $_[0];

	my $hasEcoli = 0;
	my %shigs;
	foreach my $taxid (keys %{$highestTaxa}) {
		if (exists $taxa->{'name'}->{$taxid} and $taxa->{'name'}->{$taxid} =~ m/Shigella/) {
			$shigs{$taxid} = 1;
		}
		if ($taxid == 562) {
			$hasEcoli = 1;
		}
	}
	if ($hasEcoli and scalar keys %shigs) {
		foreach my $id (keys %shigs) {
			delete $highestTaxa->{$id};
		}
	}
	
	return ($highestTaxa);

}

sub compareReads {

	### compare alignments per cluster across all sam files. 
	
	my $compare = $_[0];
	my $sam = $_[1];
	my $report = $_[2];
	my $entryQR = $_[3];
	my $referenceQR = $_[4];
	my $taxa = $_[5];
	my $bugs = $_[6];
	my $rest = $_[7];
	my $return = "";
	my $mismatchMin = "inf";
	my $goodRef = ();
	my $goodCategory = ();
	my $data = ();
	
	my $read = "";
	
	
	my $mismatchContainer = parseLines($compare, $sam, $taxa);
	my $localTaxaIDs;

	#these need to be refactored TODO
	if (scalar keys(%{$mismatchContainer}) == 0) {
		return 5;
	}
	if ($report and not exists $mismatchContainer->{$mismatchMin}->{$report}) {
		return 3
	}
	my @goodCat = keys %{$mismatchContainer->{$mismatchMin}};
	if (scalar @goodCat > 1) {
		return 2
	}
	my $forLookup = {};
	my $seenNC = 0;
	my ($highestTaxa, $spBucket, $highestRank) = getBestTaxaFromBestCategory($mismatchContainer, $goodCat[0], $taxa);
	
	
	
	if (scalar keys %{$spBucket} and not scalar keys %{$highestTaxa}) {
		$highestTaxa = $spBucket;
	}
	$highestTaxa = ecoliCheck($highestTaxa);
	my @keys = keys %{$highestTaxa};
	
	unless (scalar @keys == 1) {
		#TODO Handle non-species hit in a better way than this
		return(1);
		my $scores = {};
		my $rankScores = {};
		my $referenceCount = 0;
		my $toPutIn = {};
		foreach my $taxaID (@keys) {
			foreach my $reference (keys (%{$highestTaxa->{$taxaID}})) {
				$toPutIn->{$reference} = 1;
				$referenceCount++;
				my $trace = $taxa->traceTaxonomy($reference);
				unless ($trace) {
					$referenceCount--;
					delete $toPutIn->{$reference};
					next;
				}
				my @namedTrace;
				foreach my $id (@{$trace}) {
					my $name = $taxa->{'names'}->{$id};
					my $rankLevel = $taxa->{'acceptedRanks'}->{$taxa->{'ranks'}->{$id}};
					$rankScores->{$id}++;
					$scores->{$id} += $rankLevel;
					push @namedTrace, $name;
				}
				my $traceStr = join "->", @namedTrace;
			}
		}
		foreach my $rank (keys %{$rankScores}) {
			unless ($rankScores->{$rank} == $referenceCount) {
				delete $scores->{$rank};
			}
		}
		my @highestRank = sort {$scores->{$b} <=> $scores->{$a}} keys %{$scores};
		unless (scalar @highestRank) {
			return 1;
		}
		my $highestName = $taxa->{'names'}->{$highestRank[0]};
		foreach my $reference (keys %{$toPutIn}) {
			foreach my $start (keys %{$mismatchContainer->{$mismatchMin}->{$goodCat[0]}->{$reference}}) {
				$rest->{$highestRank[0]}->{$read}->{$reference}->{$start} = $mismatchContainer->{$mismatchMin}->{$goodCat[0]}->{$reference}->{$start};
			}
		}
	} elsif ($highestRank != 7) {
		#TODO Same as above
		return(1);
		my $taxaID = $keys[0];
		foreach my $reference (keys %{$highestTaxa->{$taxaID}}) {
			foreach my $start (keys %{$mismatchContainer->{$mismatchMin}->{$goodCat[0]}->{$reference}}) {
				$rest->{$taxaID}->{$read}->{$reference}->{$start} = $mismatchContainer->{$mismatchMin}->{$goodCat[0]}->{$reference}->{$start};
			}

		}

	} elsif(scalar(@keys) == 1) {
		my $taxaID = $keys[0];
		my $ncQR = qr/^NC/;
		foreach my $reference (keys %{$highestTaxa->{$taxaID}}) {
			if ($seenNC and not ($reference =~ m/$ncQR/)) {
				next;
			}
			foreach my $start (keys %{$mismatchContainer->{$mismatchMin}->{$goodCat[0]}->{$reference}}) {
				$bugs->{$taxaID}->{$reference}->{$read}->{$start} = $mismatchContainer->{$mismatchMin}->{$goodCat[0]}->{$reference}->{$start};
			}
		}
	} else {
		return 1;
	}
	return 0;

}

sub getNextRead {
		my $fh = $_[0];
		my $category = $_[1];
		my $currentHeader = $_[2];
		my $prevLine = $_[3];
		my $currentLines = $_[4];
		my $finished = $_[5];
		my $headQR = $_[6];
		my $line;
		my $read;
		while (1) {
			if (exists $prevLine->{$category}) {
				$line = $prevLine->{$category};
				$line =~ m/$headQR/;
				$read = $1;
				$currentHeader->{$category} = $read;
				delete $currentLines->{$category};
				delete $prevLine->{$category};
			} elsif (eof($fh->{$category})) {
				$finished->{$category} = 1;
				last
			} else {
				$line = $fh->{$category}->getline;
				chomp $line;
				$line =~ m/$headQR/;
				$read = $1;
			}
			unless (exists $currentHeader->{$category}) {
				 $currentHeader->{$category} = $read;
			}
			if ($read ne $currentHeader->{$category}) {
				$prevLine->{$category} = $line;
				last;
			} else {
				push @{$currentLines->{$category}}, $line;
				redo
			}
		}
		return($currentHeader, $prevLine, $currentLines, $finished);

}


sub mergeKeggList {
	my $keggList = shift;
	my $mergeIn = shift;
	foreach my $ko (keys %{$mergeIn}) {
		$keggList->{$ko} += $mergeIn->{$ko};
	}
}

sub readsToKeggList {
	my $keggHits = shift;
	my $keggList = {};
	foreach my $read (keys %{$keggHits}) {
		foreach my $ko (@{$keggHits->{$read}}) {
			$keggList->{$ko}++;
		}
	}
	return $keggList;
}

sub mergeKegg {
	my $keggHits = $_[0];
	my $mergeIn = $_[1];
	foreach my $read (keys %{$mergeIn}) {
		if (exists $keggHits->{$read}) {
			push @{$keggHits->{$read}}, @{$mergeIn->{$read}};
		} else {
			$keggHits->{$read} = $mergeIn->{$read};
		}
	}
}

sub calculateCoverageAndKeggHits {
	my $reads = $_[0];
	my $length = $_[1];
	my $reference = $_[2];
	my $kegg = $_[3];
	my $used = [];
	my $letters = 0;
	my $currEnd;
	my $firstStart = 0;
	my $gapAmount = 0;
	my $keggHits = {};
	$kegg->createTree($reference);
	foreach my $read (sort { (sort {$a <=> $b} keys %{$reads->{$a}})[0] <=> (sort {$a <=> $b} keys %{$reads->{$b}})[0] } keys %{$reads}) {
		push @{$used}, $read;
		my @starts = sort {$a <=> $b} keys %{$reads->{$read}};
		my $start = $starts[0];
		my $stop = $reads->{$read}->{$start};
		my $keggHit = $kegg->query($reference, $start, $stop);
		if ($keggHit != 0 and scalar(@{$keggHit})) {
			$keggHits->{$read} = $keggHit;
		}
		my $dist = $stop - $start;
		$letters += $dist;

		unless ($firstStart) {
			$firstStart = $start;
			$currEnd = $stop;
			next;
		}
		if ($start > $currEnd) {
			my $gap = $start - $currEnd;
			$gapAmount += $gap;
			$currEnd = $stop;
			next
		}
		if ($stop > $currEnd) {
			$currEnd = $stop;
		}

	}
	$kegg->removeTree($reference);
	my $covered = $currEnd - ($firstStart + $gapAmount);
	return ($covered, $letters, $used, $keggHits);

}
sub collapseHits {
	my $array = $_[0];
	my $hash = {};
	foreach my $read (@{$array}) {
		$hash->{$read} = 1;
	}
	return (scalar(keys%{$hash}));
}

##### START #####

### Set up file handles ###
my $entryQR = qr/([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t[^\t]+\t[^\t]+\t[^\t]+\t[^\t]+\t[^\t]+\t([^\t]+)\t.*NM:i:(\d+)/;
my $referenceQR = qr/.*\|([^|]+)\|/;
unless ($report) {
	$report = "";
}


my $taxa = TaxaObject->new({"nodeFile" => $nodeFile, "nameFile" => $nameFile, "prebuilt" => $prebuilt});
#TODO Ewww Goto
if ($checkSpecies and $checkRest) {
	goto CHECKPOINT;
}

my $samFiles = parseFileStr($samStr, $report);
my $sam = SamCollection->new({"isSorted" => $isSorted});
$sam->openSamFiles($samFiles);


my $bugsByTaxaIDSpecies = {};
my $bugsByTaxaIDRest = {};
my $readCount = 0;
my $goodCount = 0;
my $lowResCount = 0;
my $notReport = 0;
my $contamination = 0;
my $tooShort = 0;
my $allGone = 0;


while (1) {
	$readCount++;
	if ($readCount % 100000 == 0) {
		print STDERR "processed $readCount Reads\n";
		my $ratio = $goodCount / $readCount;
		print STDERR "ratio of good to rejected: $goodCount / $readCount = $ratio\n";
		print STDERR "Rest of the reads\n";
		print STDERR "Contamination: $contamination\n";
		print STDERR "best hit not microbial: $notReport\n";
		print STDERR "Does not hit species level: $lowResCount\n";
		print STDERR "All Reads Thrown Out: $allGone\n";
		print STDERR "Read too short: $tooShort\n\n";
	}
	my @lowestToHighest = sort {$sam->{'currentHeader'}->{$a} cmp $sam->{'currentHeader'}->{$b}} keys (%{$sam->{'currentHeader'}});
	unless (scalar(@lowestToHighest)) {
		last;
	}
	my $lowest = $sam->{'currentHeader'}->{$lowestToHighest[0]};
	my $compare = [];
	foreach my $category (@lowestToHighest) {
		if ($sam->{'currentHeader'}->{$category} eq $lowest) {
			push @{$compare}, $category;
		}
	}
	my $status = compareReads($compare, $sam, $report, $entryQR, $referenceQR, $taxa, $bugsByTaxaIDSpecies, $bugsByTaxaIDRest);
	unless ($status) {
		$goodCount++;
	} else {
		switch($status) {
			case 1 {$lowResCount++}
			case 2 {$contamination++}
			case 3 {$notReport++}
			case 4 {$tooShort++}
			case 5 {$allGone++}
		}

	}
	#advance the compared files
	foreach my $category (@{$compare}) {
		if (exists $sam->{'finished'}->{$category}) {
			delete $sam->{'currentHeader'}->{$category};
			delete $sam->{'finished'}->{$category};
			$sam->closeFile($category);
			next;
		}
		$sam->getNextRead($category);
	}
}
$sam->closeFiles();
if ($checkOutSpecies and $checkOutRest) {
	my $bugsSpecies = encode_json($bugsByTaxaIDSpecies);
	my $bugsRest = encode_json($bugsByTaxaIDRest);
	open BUGS, ">$checkOutSpecies";
	print BUGS "$bugsSpecies";
	close BUGS;
	open BUGS, ">$checkOutRest";
	print BUGS "$bugsRest";
	close BUGS;
}


# TODO Don't cheat with a GOTO
CHECKPOINT:
if ($checkSpecies and $checkRest) {
	open IN, "$checkSpecies";
	while (my $line = <IN>) {
		chomp $line;
		$bugsByTaxaIDSpecies = decode_json($line);
	}
	close IN;
	open IN, "$checkRest";
	while (my $line = <IN>) {
		chomp $line;
		$bugsByTaxaIDRest = decode_json($line);
	}
	close IN;
}

my $kegg = KeggObject->new({'keggAnnotate' => $keggFile});



my $taxaAmount = scalar(keys(%{$bugsByTaxaIDSpecies}));
my $taxaDone = 0;
my $vector;
my $collapsedTaxaInformation = {};
my $metagenomeKeggList = {};
my $ko2taxa = {};
my $taxaName = {};
foreach my $taxaID (keys %{$bugsByTaxaIDSpecies}) {
	$taxaDone++;
	my $name = $taxa->{'names'}->{$taxaID};
	$collapsedTaxaInformation->{$taxaID}->{'name'} = $name;
	$taxaName->{$taxaID} = $name;

	my $assemblies = {};
	my $assemblyLength = {};
	my $assmUsed = {};
	my $overallHits = [];
	my $referenceAmount = scalar(keys %{$bugsByTaxaIDSpecies->{$taxaID}});
	my $referenceCount = 0;
	my $avgLen = 0;
	my $overallKeggHits = {};
	my $goodCount = 0;
	foreach my $reference (keys %{$bugsByTaxaIDSpecies->{$taxaID}}) {
		unless (exists $taxa->{'knownTaxa'}->{$reference}->{'type'}) {
			next;
		} elsif (exists($taxa->{'bad'}->{$reference})) {
			next;
		}
		$referenceCount++;
		if ($taxa->{'knownTaxa'}->{$reference}->{'type'} eq 'wgs') {
			my $assemblyAccession = $taxa->{'knownTaxa'}->{$reference}->{'wgs'};
			unless (exists $assemblyLength->{$assemblyAccession} and $assemblyLength->{$assemblyAccession}) {
				$assemblyLength->{$assemblyAccession} = $taxa->{'wgs'}->{$assemblyAccession}->{'overallLength'};
			}
			my ($covered, $letters, $used, $keggHits) = calculateCoverageAndKeggHits($bugsByTaxaIDSpecies->{$taxaID}->{$reference}, $assemblyLength->{$assemblyAccession}, $reference, $kegg);
			unless (exists $assmUsed->{$assemblyAccession}) {
				$assmUsed->{$assemblyAccession} = [];
			}
			$assemblies->{$assemblyAccession}->{$reference}->{'keggHits'} = $keggHits;
			push @{$assmUsed->{$assemblyAccession}}, @{$used};
			push @{$overallHits}, @{$used};
			$assemblies->{$assemblyAccession}->{$reference}->{'stats'} = [$covered, $letters];
		} else {
			$goodCount++;
			my $length = $taxa->{'knownTaxa'}->{$reference}->{'length'};
			my ($covered, $letters, $used, $keggHits) = calculateCoverageAndKeggHits($bugsByTaxaIDSpecies->{$taxaID}->{$reference}, $length, $reference, $kegg);
			push @{$overallHits}, @{$used};
			my $coverageAmount = $covered/$length;
			my $DOC = $letters/$covered;
			my $readHits = scalar(@{$used});
			$avgLen += $length;
			mergeKegg($overallKeggHits, $keggHits);
			$collapsedTaxaInformation->{$taxaID}->{'references'}->{$reference}->{'depth'} = $DOC;
			$collapsedTaxaInformation->{$taxaID}->{'references'}->{$reference}->{'coverage'} = $coverageAmount;
		}
	}
	unless($referenceCount) {
		delete ($collapsedTaxaInformation->{$taxaID});
		next;
	}
	my $totalHits = collapseHits($overallHits);
	foreach my $assembly (keys %{$assemblies}) {
		$goodCount++;
		my $totalCoverage = 0;
		my $totalLetters = 0;
		my $overallLength = $assemblyLength->{$assembly};
		my $keggHits = {};
		foreach my $contig (keys %{$assemblies->{$assembly}}) {
			$totalCoverage += $assemblies->{$assembly}->{$contig}->{'stats'}->[0];
			$totalLetters += $assemblies->{$assembly}->{$contig}->{'stats'}->[1];
			mergeKegg($overallKeggHits, $assemblies->{$assembly}->{$contig}->{'keggHits'});
		}
		$avgLen += $overallLength;
		my $coverageAmount = $totalCoverage/$overallLength;
		my $DOC = $totalLetters/$totalCoverage;
		$collapsedTaxaInformation->{$taxaID}->{'references'}->{$assembly}->{'depth'} = $DOC;
		$collapsedTaxaInformation->{$taxaID}->{'references'}->{$assembly}->{'coverage'} = $coverageAmount;
	}
	$collapsedTaxaInformation->{$taxaID}->{'totalHits'} = $totalHits;
	$collapsedTaxaInformation->{$taxaID}->{'averageGenomeLength'} = int($avgLen / $goodCount);
	my $keggList = readsToKeggList($overallKeggHits);
	foreach my $ko (keys %{$keggList}) {
		$ko2taxa->{$ko}->{$taxaID} = $keggList->{$ko};
	}
	$collapsedTaxaInformation->{$taxaID}->{'keggList'} = $keggList;
	mergeKeggList($metagenomeKeggList, $keggList);
}
undef($taxa);
$kegg->parseModules($modules);
foreach my $taxaID (keys %{$collapsedTaxaInformation}) {
	my $keggList = $collapsedTaxaInformation->{$taxaID}->{'keggList'};
	if (scalar(%{$keggList})) {
		$collapsedTaxaInformation->{$taxaID}->{'kegg'} = $kegg->getModuleInfo($keggList);
	} else {
		$collapsedTaxaInformation->{$taxaID}->{'kegg'} = 0;
	}
}
my $metagenomeModules = 0;
if (scalar(%{$metagenomeKeggList})) {
	$metagenomeModules = $kegg->getModuleInfo($metagenomeKeggList);
}
my $out = {};
$out->{'metagenomeModules'} = $metagenomeModules;
$out->{'taxaInfo'} = $collapsedTaxaInformation;
$out->{'ko2taxa'} = $ko2taxa;

#browseErr($out);
my $outStr = encode_json($out);
open OUT, ">$output";
print OUT "$outStr";
close OUT;





















__END__
