package KeggObject;

#use DataBrowser qw(browse browseErr);
use Set::IntervalTree;
use Storable qw(dclone);
use warnings;
use strict;


sub new {
	my $class = shift;
	my $settings;
	if (scalar(@_)) {
		$settings = shift;
	}
	my $self = {};
	$self->{'raw'} = {};
	$self->{'genomes'} = {};
	my $keggFile;
	unless (exists($settings->{'keggAnnotate'})) {
		die "KEGG annotation file not found\n";
	} else {
		$keggFile = $settings->{'keggAnnotate'};
	}
	open IN, "$keggFile" or die "Can't open KEGG annotations: $keggFile\n";
	my $lineQR = qr/([^\t]+)*\t(.*)/;
	while (my $line = <IN>) {
		chomp $line;
		$line =~ m/$lineQR/;
		$self->{'raw'}->{$1} = $2;
	}
	close IN;
	bless($self, $class);
	return $self;

}

sub createTree {
	my $self = shift;
	my $genome = shift;
	unless (exists $self->{'raw'}->{$genome}) {
		return 0;
	}
	$self->{'genomes'}->{$genome} = Set::IntervalTree->new;
	my @parts = split /\|/, $self->{'raw'}->{$genome};
	my $partsQR = qr/(.*)_(.*)_(.*)/;
	foreach my $part (@parts) {
		$part =~ m/$partsQR/;
		my $KO = $3;
		my $start = $1;
		my $stop = $2;
		$self->{'genomes'}->{$genome}->insert($KO, $start, $stop);
	}
}

sub removeTree {
	my $self = shift;
	my $genome = shift;
	delete($self->{'genomes'}->{$genome});

}

sub query {
	my $self = shift;
	my $genome = shift;
	my $start = shift;
	my $stop = shift;
	unless (exists $self->{'genomes'}->{$genome}) {
		return(0)
	}
	return ($self->{'genomes'}->{$genome}->fetch($start, $stop));
}
sub getModuleInfo {
	my $self = shift;
	my $geneList = shift;
	my $totalHits = {};
	my $modules = dclone($self->{'modules'});
	foreach my $module (keys %{$modules}) {
		my $count = 0;
		my $complete = 0;
		my @stepCompleteness;
		foreach my $step (@{$modules->{$module}->{'pathway'}}) {
			$count++;
			$self->fillStep($step, $geneList, 0);
			my $subcomplete = $self->determineStepCompleteness($step);
			unless (scalar (@{$subcomplete})) {
				$count--;
				next;
			}
			my $localComplete = 0;
			foreach my $elem (@{$subcomplete}) {
				$localComplete += $elem;
			}
			my $total = scalar (@{$subcomplete});
			my $localCompleteness = $localComplete / $total;
			push @stepCompleteness, $localCompleteness;;
			$complete += $localCompleteness;
		}
		my $completeness = $complete/$count;
		if ($completeness == 0) {
			delete $modules->{$module};
			next;
		}
		$modules->{$module}->{'completeness'} = $completeness;
		$modules->{$module}->{'step_completeness'} = \@stepCompleteness;
		my $hitCount = 0;
		foreach my $member (@{$modules->{$module}->{'members'}}) {
			if (exists $geneList->{$member}) {
				$hitCount += $geneList->{$member};
			}
		}
		if ($hitCount == 0) {
			delete $modules->{$module};
			next;
		}
		delete ($modules->{$module}->{'pathway'});
		$modules->{$module}->{'hits'} = $hitCount;
		$totalHits->{$modules->{$module}->{'type'}} += $hitCount;
	}
	return($modules);

}
sub parseModules {
	my $self = shift;
	my $in = shift;
	my $modules = {};
	my $types = {};
	my $skip = 0;
	open IN, "$in" or die "Can't open KEGG modules file\n";
	my $currentName = "";
	my $entryQR = qr/^ENTRY\s+(\S+)\s+(\S+)\s+/;
	my $nameQR = qr/^NAME\s+(.*)/;
	my $defQR = qr/^DEFINITION\s+(.*)/;
	my $mQR = qr/^M/;
	my $currentType;
	while (my $line = <IN>) {
		chomp $line;
		if ($line =~ m/$entryQR/) {
			my $name = $1;
			my $type = $2;
			$currentType = $type;
			if ($name =~ m/$mQR/) {
				$skip = 0;
			} else {
				$skip = 1;
				next;
			}
			$currentName = $name;
			$modules->{$name} = {};;
			$modules->{$name}->{'type'} = $type;
			$types->{$type}->{$name} = 1;
			next;
		}
		if ($skip) {
			next;
		}
		if ($line =~ m/$nameQR/){
			$modules->{$currentName}->{'name'} = $1;
			next;
		}
		if ($line =~ m/$defQR/) {
			my $defline = $1;
			my $checkTopSplit = $self->checkForTopLevelBalancedOr($defline);
			if (scalar(@{$checkTopSplit}) == 1) {
				my ($members, $parse) = $self->doDefElement($defline);
				$modules->{$currentName}->{'pathway'} = $parse;
				$modules->{$currentName}->{'members'} = $members;
			} else {
				my $count = 0;
				my $engName = $modules->{$currentName}->{'name'};
				delete $modules->{$currentName};
				delete $types->{$currentType}->{$currentName};
				foreach my $str (@{$checkTopSplit}) {
					$str =~ s/^\(//g;
					$str =~ s/\)$//g;
					my $name = $currentName . ".$count";
					my ($members, $parse) = $self->doDefElement($str);
					$modules->{$name}->{'pathway'} = $parse;
					$modules->{$name}->{'members'} = $members;
					$modules->{$name}->{'name'} = $engName;
					$modules->{$name}->{'type'} = $currentType;
					$types->{$currentType}->{$name} = 1;
					$count++;
				}
			}
		}

	}
	close IN;
	$self->{'modules'} = $modules;
	$self->{'types'} = $types;
	#browse($self->{'modules'});
	#browse($self->{'types'});
}

sub doDefElement {
	my $self = shift;
	my $defline = shift;
	my $tempLine = $defline;
	$tempLine =~ s/[^K0-9]/ /g;
	$tempLine =~ s/^\s+//g;
	$tempLine =~ s/\s+$//g;
	my @members = split /\s+/, $tempLine;
	my $parse = $self->parseDefinitionLine($defline);
	return (\@members, $parse);
}

sub checkForTopLevelBalancedOr {
	my $self = shift;
	my $line = shift;
	my $count = 0;
	my $seenParen = 0;
	my @letters = split //, $line;
	my $levels = ();
	my $parenCount = 0;
	foreach my $letter (@letters) {
		if (not $seenParen and $letter eq "(") {
			$seenParen = 1;
			$parenCount = 1;
		} elsif ($letter eq "(") {
			$parenCount++;
		} elsif ($letter eq ")") {
			$parenCount--;
		}
		if ($letter eq "," and $parenCount == 0) {
			$count++;
			next;
		}
		$levels->[$count] .= $letter;
	}
	return($levels);
}

sub parseDefinitionLine {
	my $self = shift;
	my $line = shift;;
	my $results = {};
	my $leftCount = $line =~ tr/(//;
	my $rightCount = $line =~ tr/)//;
	my $diff = $leftCount - $rightCount;
	unless ($diff == 0) {
		if ($diff > 0) {
			for (my $i = 0; $i < $diff; $i++) {
				$line =~ s/^\(//;
			}
		} else {
			for (my $i = 0; $i < ($diff * -1); $i++) {
				$line =~ s/\}$//g;
			}
		}
	}
	my @steps = split /\s+/, $line;
	my $steps = $self->fixSteps(\@steps, " ");
	my @resolvedSteps;
	foreach my $step (@{$steps}) {
		my @temp;
		push @temp, $step;
		push @resolvedSteps, $self->parseSteps(\@temp, "");
	}
	unless(scalar(@resolvedSteps)) {
		die "couldn't parse\n";
	}
	return(\@resolvedSteps);
}

sub parseSteps {
	my $self = shift;
	my @steps = @{$_[0]};
	my $tab = $_[1];
	my $opt = $_[2];
	my @delimiters = (',', '+', ' ', '-');
	my $stepcount = 0;
	my $resolvedSteps = {};
	foreach my $step (@steps) {
		my $typeAdd = "";
		if ($opt and $stepcount > 0) {
			$typeAdd = ".OPT";
		}
		$stepcount++;
		my @substeps;
		my $overall = {};
		if ($step =~ m/^\(.*\)$/) {
			$step =~ s/^\(//g;
			$step =~ s/\)$//g;
		}
		my $leftCount = $step =~ tr/(//;
		my $rightCount = $step =~ tr/)//;
		unless ($leftCount == $rightCount) {
			die "$step is not balanced\n";
		}
		if ($leftCount) {
			$step =~ m/\(/;
			my $leftStart = $-[0];
			$step =~ m/\)/;
			my $rightStart = $-[0];
			if ($leftStart > $rightStart) {
				$step = "($step)";
			}
		}
		my $firstDelim = "";
		my @possibilities;
		my $poss;
		foreach my $delimiter (@delimiters) {
			@possibilities = split "\Q$delimiter\E", $step;
			$poss = $self->fixSteps(\@possibilities, "$delimiter");
			if (scalar(@{$poss}) > 1) {
				$firstDelim = $delimiter;
				last
			}
		}
		unless ($firstDelim) {
			if ($step =~ m/(.)(K.*)/) {
				my $type = $1;
				my $knum = $2;
				my $newtab = $tab . "\t";
				if ($type eq "-" or $type eq " ") {
					$resolvedSteps->{"OPT"} = $self->parseSteps(["$knum.OPT"], "$newtab");
				} else {
					print STDERR "$step\n";
					die "invalid prefix in pathway math\n";
				}
			} else {
				$step = $step . "$typeAdd";
				$resolvedSteps->{$step} = 0;
			}
		} else {
			my $newtab = $tab . "\t";
			my $type;
			my $optStatus = "";
			if ($firstDelim eq "-") {
				$type = "OPT";
				$optStatus = 1;
			} elsif ($firstDelim eq ",") {
				$type = "OR";
			} else {
				$type = "AND";
			}
			$type = $type . "$typeAdd";
			if (exists $resolvedSteps->{$type}) {
				my $count = 2;
				my $origType = $type;
				while (exists $resolvedSteps->{$type}) {
					$type = $origType . "$count";
					$count++;
				}
				$resolvedSteps->{$type} = $self->parseSteps($poss, "$newtab", $optStatus);
			} else {
				$resolvedSteps->{$type} = $self->parseSteps($poss, "$newtab", $optStatus);
			}


		}
	}
	return($resolvedSteps);
}

sub fixSteps {
	my $self = shift;
	my @steps = @{$_[0]};
	my $join = $_[1];
	my @origSteps;
	my @temp;
	my $parenCount = 0;
	foreach my $elem (@steps) {
		unless ($elem) {
			next;
		}
		my $leftCount = $elem =~ tr/(//;
		my $rightCount = $elem =~ tr/)//;
		my $overallCount = $leftCount - $rightCount;
		$parenCount += $overallCount;
		if ($parenCount == 0) {
			push @temp, $elem;
			my $string = join "$join", @temp;
			push @origSteps, $string;
			@temp = ();
		} else {
			push @temp, $elem;
		}
	}
	return(\@origSteps);
}

sub fillStep {
	my $self = shift;
	my $step = $_[0];
	my $geneList = $_[1];
	my $minCount = $_[2];
	my $tab = $_[3];
	unless ($tab) {
		$tab = "";
	}
	my $optQR = qr/(\S+)\.OPT/;
	foreach my $substep (keys %{$step}) {
		my $reference = ref($step->{$substep});
		if ($reference) {
			my $newtab = $tab . "\t";
			$self->fillStep($step->{$substep}, $geneList, 0, $newtab);
		} else {
			my $key = $substep;
			if ($key =~ m/$optQR/) {
				$key = $1;
			}
			if (exists $geneList->{$key}) {
				$step->{$substep} = $geneList->{$key};
			} else {
				$step->{$substep} = 0;
			}
		}
	}
}

sub optAndHandling {
	my $self = shift;
	my $subcomplete = $_[0];
	my $good = 0;
	foreach my $elem (@{$subcomplete}) {
		$good += $elem;
	}
	my $total = scalar (@{$subcomplete});
	my $localComplete = $good/$total;
	return ($localComplete);

}

sub determineStepCompleteness {
	my $self = shift;
	my $step = $_[0];
	my $count = 0;
	my $complete = 0;
	my @partCompleteness;
	foreach my $part (keys %{$step}) {
		my $reference = ref($step->{$part});
		if ($reference) {
			if ($part =~ m/OPT/) {
				my %temp = %{$step->{$part}};
				foreach my $key (keys %temp) {
					if ($key =~ m/OPT/) {
						delete $temp{$key};
					}
				}
				unless(scalar(keys %temp)) {
					next;
				}
				my $subcomplete = $self->determineStepCompleteness(\%temp);
				my $localComplete = $self->optAndHandling($subcomplete);
				push @partCompleteness, $localComplete;
			}
			if ($part =~ m/AND/) {
				my $subcomplete = $self->determineStepCompleteness($step->{$part});
				my $localComplete = $self->optAndHandling($subcomplete);
				push @partCompleteness, $localComplete;
			}
			if ($part =~ m/OR/) {
				my $subcomplete = $self->determineStepCompleteness($step->{$part});
				my $max = 0;
				foreach my $elem (@{$subcomplete}) {
					if ($elem > $max) {
						$max = $elem;
					}
				}
				push @partCompleteness, $max;
			}
		} else {
			if ($step->{$part}) {
				push @partCompleteness, 1;
			} elsif ($part eq "--") {
				push @partCompleteness, 1;
			} else {
				push @partCompleteness, 0;
			}
		}
	}
	return (\@partCompleteness);
}

return 1;
