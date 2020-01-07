package TaxaObject;

#use DataBrowser qw(browse browseErr);
use XML::Hash;
use JSON::XS;
use Try::Tiny;
use LWP::UserAgent;
use Time::HiRes qw(ualarm);;
use warnings;
use strict;

$SIG{ALRM} = sub { die };

sub new {
	my $class = shift;
	my $settings = {};
	if (scalar (@_)) {
		$settings = shift;
	}
	my $self = {};
	$self->{'knownTaxa'} = {};
	$self->{'names'} = {};
	$self->{'parents'} = {};
	$self->{'ranks'} = {};
	$self->{'xml'} = XML::Hash->new();
	$self->{'taxaTrace'} = {};
	$self->{'wgs'} = {};
	$self->{'genomes'} = {};
	$self->{'nullArray'} = ();
	$self->{'nullHash'} = {};
	$self->{'bad'} = {};
	$self->{'ua'} = LWP::UserAgent->new;
	$self->{'ua'}->timeout(12);
	$self->{'acceptedRanks'} = {
	"kingdom" => 1,
	"phylum" => 2,
	"class" => 3,
	"order" => 4,
	"family" => 5,
	"genus" => 6,
	"species" => 7
	};
	if (exists($settings->{'prebuilt'})) {
		open PREBUILT, "$settings->{'prebuilt'}" or die "can't open prebuilt file\n";
		while (my $line = <PREBUILT>) {
			chomp $line;
			my @parts = split /\t/, $line;
			$self->{$parts[0]} = decode_json($parts[1]);
		}
	}
	foreach my $known (keys %{$self->{'knownTaxa'}}) {
			delete $self->{'knownTaxa'}->{'taxaTrace'};
	}
	if (exists ($settings->{'nodeFile'})) {
		createNodes($self, $settings->{'nodeFile'});
	}
	if (exists($settings->{'nameFile'})) {
		makeNames($self, $settings->{'nameFile'});
	}
	bless $self, $class;
	return $self;

}
sub clearBad {
	my $self = shift;
	$self->{'bad'} = {};
}

sub highestRank {
	my $self = shift;
	my $refName = shift;
	my $status = $self->getTaxaInfo($refName);
	$self->traceTaxonomy($refName);
	return ($self->{'knownTaxa'}->{$refName}->{'highestRank'});
}

sub getLength {
	my $self = shift;
	my $refName = shift;
	my $status = $self->getTaxaInfo($refName);
	unless ($status) {
		return 0;
	}
	return ($self->{'knownTaxa'}->{$refName}->{'length'});
}

sub getTaxaId {
	my $self = shift;
	my $refName = shift;
	my $taxaTrace = $self->traceTaxonomy($refName);
	if ($taxaTrace) {
		return $taxaTrace->[-1];
	} else {
		return 0;
	}
}


sub createNodes {
	my $self = shift;
	my $nodesFile = shift;
	open NODES, "$nodesFile" or die "Can't open Nodes File";
	my $nodeQR = qr/([^\t]*)\t[^\t]*\t([^\t]*)\t[^\t]*\t([^\t]*)\t/;
	while (my $line = <NODES>) {
		$line =~ m/$nodeQR/;
		my $node = $1;
		my $parent = $2;
		my $rank = $3;
		$self->{'parents'}->{$node} = $parent;
		$self->{'ranks'}->{$node} = $rank;
	}
	close NODES;
}

sub makeNames {
	my $self = shift;
	my $namesFile = shift;
	open NAMES, "$namesFile";
	my $nameQR = qr/([^\t]*)\t[^\t]*\t([^\t]*)\t[^\t]*\t[^\t]*\t[^\t]*([^\t]*)\t/;
	while (my $line = <NAMES>) {
		chomp $line;
		my @parts = split "\t", $line;
		my $nodeId = $parts[0];
		my $name = $parts[2];
		my $type = $parts[6];
		unless ($type eq "scientific name") {
			next;
		}
		$self->{'names'}->{$nodeId} = $name;
	}
	close NAMES;
}


sub getTaxaInfo {
	my $self = shift;
	my $refName = shift;
	if (exists $self->{'bad'}->{$refName}) {
		return 0;
	}
	if (exists $self->{'knownTaxa'}->{$refName} and scalar($self->{'knownTaxa'}->{$refName})) {
		return 1;
	}
	$self->lookUpTaxa($refName);
	return 1;
}


sub lookUpTaxa {
	my $self = shift;
	my $refName = $_[0];
	my $inTry = shift;
	my $jsonHash = $self->ncbiFetchSummary($refName);
	unless ($jsonHash) {
		return 1;
	}
	unless (exists $jsonHash->{'result'}->{'uids'} and defined $jsonHash->{'result'}->{'uids'}) {
		return 1;
	}
	if (scalar (@{$jsonHash->{'result'}->{'uids'}}) > 1) {
		print STDERR "TAXA Object Error: $refName had more than 1 result:\n";
		print STDERR "Using first ID only: $jsonHash->{'result'}->{'uids'}->[0]\n";
		print STDERR "END TAXA Object Error\n\n";
	}
	unless (scalar(@{$jsonHash->{'result'}->{'uids'}})) {
		unless (defined($inTry)) {
			$inTry = 11;
		} elsif ($inTry == 0) {
			return 1;
		}
		$inTry--;
		sleep (12 - $inTry);
		return($self->lookUpTaxa($refName, $inTry));
	}
	my $id = $jsonHash->{'result'}->{'uids'}->[0];
	my $result = $jsonHash->{'result'}->{$id};
	my $type = "";
	if ($result->{'genome'}) {
		$type = $result->{'genome'};
	}
	if ($type =~ m/plasmid/) {
		$self->{'knownTaxa'}->{$refName}->{'type'} = "plasmid";
		return 1;
	}
	my $completeness = $result->{'completeness'};
	if ($completeness ne "complete" and $refName =~ m/NZ_[A-Z]{4}/) {
		my $assemblyAccession = $result->{'assemblyacc'};
		unless ($assemblyAccession) {
			$refName =~ m/NZ_(.*)\./;
			$assemblyAccession = $1;
		}
		unless ($assemblyAccession) {
			$self->{'bad'}->{$refName} = 1;
			return 1;
		}
		$assemblyAccession =~ s/\d\d\d\d\d\d$/000000/g;
		$self->fillInWgs($refName, $assemblyAccession);
	} elsif ($completeness ne "complete" and $refName =~ m/NZ_[A-Z]{2}/ and $type ne "chromosome") {
		my $extra = $result->{'extra'};
		my $wgsProj;
		my $assemblyAccession;
		if ($extra =~ m/WGS:([^\|]*)\|/) {
			$wgsProj = $1;
			$wgsProj =~ s/^.._//g;
			$wgsProj =~ s/\..*$//g;
		}
		if ($wgsProj) {
			if (length($wgsProj) == 4) {
				my $try = 1;
				my $found = 0;
				while ($try < 100) {
					my $num = sprintf("%02d", $try);
					my $tryAssm = $wgsProj . $num . "000000";
					my $tryResult = $self->ncbiSearch($tryAssm);
					if ($tryResult and not $found) {
						$assemblyAccession = $tryAssm;
						$found = 1;
					} elsif ($found and not $tryResult) {
						last
					} elsif ($found and $tryResult) {
						$assemblyAccession = $tryAssm;
					}
					$try++;
				}
			} elsif (length($wgsProj) == 6) {
				$assemblyAccession = "$wgsProj" . "000000";
			}
			if ($assemblyAccession) {
				$self->fillInWgs($refName, $assemblyAccession);
			} else {
				$self->{'bad'}->{$refName} = 1;
				print STDERR "unable to find assembly accession: $refName\n";
				print STDERR "Extra for $refName: $extra\n";
			}
		}
		if (not $assemblyAccession and  $result->{'assemblyacc'}) {
			$assemblyAccession = $result->{'assemblyacc'};
			$self->fillInWgs($refName, $assemblyAccession);
		}
	} else {
		if ($completeness eq "complete" or $refName =~ m/^NC/ or $type eq "chromosome") {
			$self->{'knownTaxa'}->{$refName}->{'type'} = "genome";
			$self->{'knownTaxa'}->{$refName}->{'length'} = $result->{'slen'};
			$self->{'knownTaxa'}->{$refName}->{'taxid'} = $result->{'taxid'};
		}
	}
	return 1;
}

sub fillInWgs {
	my $self = shift;
	my $refName = shift;
	my $assemblyAccession = shift;
	my $contigNameCanonical = $refName;
	$contigNameCanonical =~ s/^.._//g;
	if (exists $self->{'wgs'}->{$assemblyAccession} and not (exists $self->{'knownTaxa'}->{$refName})) {
		$self->{'knownTaxa'}->{$refName}->{'type'} = 'wgs';
		$self->{'knownTaxa'}->{$refName}->{'wgs'} = $assemblyAccession;
		$self->{'knownTaxa'}->{$refName}->{'length'} = $self->{'wgs'}->{$assemblyAccession}->{'overallLength'};
		$self->{'knownTaxa'}->{$refName}->{'taxid'} = $self->{'wgs'}->{$assemblyAccession}->{'taxid'};
		return 1;
	}
	my $assmHash = $self->ncbiFetchSummary($assemblyAccession);
	unless ($assmHash) {
		$self->{'bad'}->{$refName} = 1;
		return 0;
	}
	$self->{'knownTaxa'}->{$refName}->{'type'} = "wgs";
	$self->{'wgs'}->{$assemblyAccession} = {};
	$self->{'knownTaxa'}->{$refName}->{'wgs'} = $assemblyAccession;
	my $assmResult = $assmHash->{'result'};
	my $contigCount = scalar(@{$assmResult->{'uids'}});
	foreach my $contig (@{$assmResult->{'uids'}}) {
		my $accession = $assmResult->{$contig}->{'accessionversion'};
		my $length = $assmResult->{$contig}->{'slen'};
		my $taxid = $assmResult->{$contig}->{'taxid'};
		$self->{'knownTaxa'}->{$refName}->{'length'} += $length;
		$self->{'wgs'}->{$assemblyAccession}->{'contigs'}->{$accession}->{'length'} = $length;
		$self->{'wgs'}->{$assemblyAccession}->{'taxid'} = $taxid;
		$self->{'knownTaxa'}->{$refName}->{'taxid'} = $taxid;
	}
	$self->{'wgs'}->{$assemblyAccession}->{'overallLength'} = $self->{'knownTaxa'}->{$refName}->{'length'};
	return 1;
}

sub ncbiFetchSummary {
	my $self = shift;
	my $id = shift;
	my $list = $self->ncbiSearch($id);
	unless($list) {
		return 0;
	}
	return $self->ncbiGrab($list);

}

sub ncbiSearch {
	my $self = shift;
	my $id = shift;
	my $url = "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=nucleotide&term=$id&usehistory=y&retmax=100000";
	my $try = 10;
	my $single = 0;
	my @list;
	while ($try) {
		$try--;
		my $sleeptime = 11 - (11 - $try);
		my $response;
		$response = try {
			eval {
				alarm(15);
				return $self->{'ua'}->get($url);
				alarm(0);
			}
		} catch {0};
		alarm(0);
		my $returnXml;
		if ($response and $response->is_success) {
			$returnXml = $self->{'xml'}->fromXMLStringtoHash($response->decoded_content);  # or whatever
		} else {
			my $sleeptime = 11 - (11 - $try);
			sleep ($sleeptime);
		}
		unless (ref($returnXml) eq "HASH") {
			sleep ($sleeptime);
		}
		if (not exists($returnXml->{'eSearchResult'}->{'IdList'}->{'Id'})) {
			return 0;
		}
		if (ref($returnXml->{'eSearchResult'}->{'IdList'}->{'Id'}) eq "ARRAY") {
			foreach my $accession (@{$returnXml->{'eSearchResult'}->{'IdList'}->{'Id'}}) {
				push @list, $accession->{'text'};
			}
			last;
		} else {
			push @list, $returnXml->{'eSearchResult'}->{'IdList'}->{'Id'}->{'text'};
			last;
		}
	}
	unless (scalar(@list)) {
		return 0
	} else {
		return \@list;
	}

}

sub ncbiGrab {
	my $self = shift;
	my $list = shift;
	my @list = @{$list};
	my $returnHash;
	if (scalar(@list) < 10) {
		my $idStr = join ",", @list;
		$returnHash = $self->ncbiGet($idStr);
	} else {
		my $idStr = join ",", @list;
		my @localList;
		my $count = 1;
		$returnHash = {};
		$returnHash->{'result'}->{'uids'} = ();
		while (scalar(@list) > 0) {
			my $item = pop @list;
			push @localList, $item;
			if ($count == 160) {
				my $idStr = join ",", @localList;
				my $localReturn = $self->ncbiPost($idStr);
				$self->mergeResultsList($returnHash, $localReturn);
				@localList = ();
				$count = 0;
				sleep 2;
			}
			$count++;
		}
		if (scalar(@localList)) {
			push @list, @localList;
		}
		if (scalar(@list)) {
			my $idStr = join ",", @localList;
			my $localReturn = $self->ncbiPost($idStr);
			$self->mergeResultsList($returnHash, $localReturn);
		}
	}
	return $returnHash;
}

sub mergeResultsList {
	my $self = shift;
	my $returnHash = shift;
	my $local = shift;
	unless ($local) {
		return
	}
	push @{$returnHash->{'result'}->{'uids'}}, @{$local->{'result'}->{'uids'}};
	foreach my $uid (@{$local->{'result'}->{'uids'}}) {
		$returnHash->{'result'}->{$uid} = $local->{'result'}->{$uid};
	}
}


sub parseJson {
	my $self = shift;
	my $return = shift;
	unless ($return) {
		return(0);
	}
	$return =~ s/ true,/ "true",/g;
	$return =~ s/ false,/ "false",/g;
	my $returnHash = try {decode_json($return)} catch {0};
	return ($returnHash);
}

sub ncbiPost {
	my $self = shift;
	my $idStr = shift;
	my $url = "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=nuccore&rettype=gb&retmode=json";
	my $req = new HTTP::Request POST => "$url";
	$req->content_type('application/x-www-form-urlencoded');
	$req->content("&id=$idStr");
	my $try = 10;
	my $returnHash;
	while ($try) {
		$try--;
		my $response;
		$response = try{eval {
			alarm(15);
			return($self->{'ua'}->request($req));
			alarm(0);
		}} catch {0};
		alarm(0);
		my $return;
		if ($response and $response->is_success) {
			$return = $response->content;
		} else {
			my $code = $response->code();
			if ($code == 502 or $code == 500) {
				my @list = split ",", $idStr;
				my $listLength = scalar(@list);
				my @secondHalf = splice( @list, int($listLength/2));
				unless (scalar(@list) and scalar(@secondHalf)) {
					return(0);
				}
				my $newLen1 = scalar(@list);
				my $newLen2 = scalar(@secondHalf);
				my $idStr1 = join ",", @list;
				my $idStr2 = join ",", @secondHalf;
				my $return = $self->ncbiPost($idStr1);
				my $mergeIn = $self->ncbiPost($idStr2);
				$self->mergeResultsList($return, $mergeIn);
				return($return);
			}
		}
		$returnHash = try {$self->parseJson($return)} catch {0};
		if (ref($returnHash) eq "HASH") {
			last;
		} else {
			my $sleeptime = 1 * (11 - $try);
			sleep($sleeptime);
		}
	}
	return($returnHash);
}

sub ncbiGet {
	my $self = shift;
	my $id = shift;
	my $url = "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=nuccore&id=$id&rettype=gb&retmode=json";
	my $try = 10;
	my $returnHash;
	while ($try) {
		$try--;
		my $response;
		$response = try{eval {
			alarm(15);
			return($self->{'ua'}->get($url));
			alarm(0);
		}} catch {0};
		alarm(0);
		my $return;
		if ($response and $response->is_success) {
			$return = $response->decoded_content;
		}
		$returnHash = try {$self->parseJson($return)} catch {0};
		if (ref($returnHash) eq "HASH") {
			last;
		} else {
			my $sleeptime = 1 * (11 - $try);
			sleep $sleeptime;
		}
	}
	return($returnHash);
}

sub traceTaxonomy {
	my $self = shift;
	my $refName = shift;
	unless (exists $self->{'knownTaxa'}->{$refName}) {
		my $status = $self->getTaxaInfo($refName);
		unless ($status) {
			return (0);
		}
	}
	unless (exists $self->{'knownTaxa'}->{$refName}->{'taxid'} and $self->{'knownTaxa'}->{$refName}->{'taxid'}) {
		$self->{'knownTaxa'}->{$refName}->{'taxaTrace'} = $self->{'nullArray'};
		return ($self->{'knownTaxa'}->{$refName}->{'taxaTrace'});
	}
	if (exists $self->{'knownTaxa'}->{$refName}->{'taxaTrace'}) {
		return ($self->{'knownTaxa'}->{$refName}->{'taxaTrace'});
	}
	my $taxid = $self->{'knownTaxa'}->{$refName}->{'taxid'};
	my $currentNode = $taxid;
	my $currentRank = $self->{'ranks'}->{$taxid};
	unless ($currentNode and $currentRank) {
		$self->{'knownTaxa'}->{$refName}->{'taxaTrace'} = $self->{'nullArray'};
		return ($self->{'knownTaxa'}->{$refName}->{'taxaTrace'});
	}
	my @trace;
	my $highestRank;
	my $highestRankID;
	if (exists($self->{'acceptedRanks'}->{$currentRank})) {
		unless($highestRank) {
			$highestRank = $self->{'acceptedRanks'}->{$currentRank};
			$highestRankID = $currentNode;
			$self->{'knownTaxa'}->{$refName}->{'highestRank'} = $highestRank;
			if (exists $self->{'taxaTrace'}->{$currentNode}) {
				$self->{'knownTaxa'}->{$refName}->{'taxaTrace'} = $self->{'taxaTrace'}->{$currentNode};
				return ($self->{'knownTaxa'}->{$refName}->{'taxaTrace'});
			}
		}
		push @trace, $currentNode;
	}
	my $tryCount = 100;
	while ($currentNode != 1) {
		$tryCount--;
		if ($tryCount == 0) {
			$self->{'knownTaxa'}->{$refName}->{'taxaTrace'} = $self->{'nullArray'};
			return ($self->{'knownTaxa'}->{$refName}->{'taxaTrace'});
		}
		$currentNode = $self->{'parents'}->{$currentNode};
		$currentRank = $self->{'ranks'}->{$currentNode};
		if (exists($self->{'acceptedRanks'}->{$currentRank}) and $currentNode != 1) {
			unless($highestRank) {
				$highestRank = $self->{'acceptedRanks'}->{$currentRank};
				$highestRankID = $currentNode;
				$self->{'knownTaxa'}->{$refName}->{'highestRank'} = $highestRank;
				if (exists $self->{'taxaTrace'}->{$currentNode}) {
					$self->{'knownTaxa'}->{$refName}->{'taxaTrace'} = $self->{'taxaTrace'}->{$currentNode};
					return ($self->{'knownTaxa'}->{$refName}->{'taxaTrace'});
				}
			}
			push @trace, $currentNode;
		}
	}
	@trace = reverse(@trace);
	$self->{'taxaTrace'}->{$highestRankID} = \@trace;
	$self->{'knownTaxa'}->{$refName}->{'taxaTrace'} = $self->{'taxaTrace'}->{$highestRankID};
	$self->{'knownTaxa'}->{$refName}->{'speciesTaxaID'} = $highestRankID;
	return ($self->{'knownTaxa'}->{$refName}->{'taxaTrace'});
}

return (1);
