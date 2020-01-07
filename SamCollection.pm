package SamCollection;

use warnings;
use strict;

sub new {
	my $class = shift;
	my $settings = {};
	if (scalar(@_)) {
		$settings = shift;
	}

	my $self = {};
	$self->{'currentLines'} = {};
	$self->{'prevLine'} = {};
	$self->{'currentHeader'} = {};
	$self->{'lastHeader'} = {};
	$self->{'finished'} = {};
	$self->{'fh'} = {};
	$self->{'categories'} = {};
	$self->{'headQR'} = qr/([^\t]+)\t/;
	$self->{'isSorted'} = 0;
	if (exists $settings->{'isSorted'} and $settings->{'isSorted'}) {
		$self->{'isSorted'} = 1;
	}

	bless $self, $class;

	return $self;
}


sub openSamFiles {

	### takes in a hash of arrays ###
	#HASH
	#category1 -> ARRAY
	#	file1.txt
	#	file2.txt
	#category2 -> ARRAY
	#	file3.txt
	#	file4.txt

	my $self = shift;
	my $samFiles = $_[0];
	$ENV{'LC_ALL'} = "POSIX";
	foreach my $category (keys %{$samFiles}) {
		if (exists $self->{'categories'}->{$category}) {
			die "$category file handle already open\n";
		} else {
			$self->{'categories'}->{$category} = 1;
		}
		my $files = join " ", @{$samFiles->{$category}};
		unless ($self->{'isSorted'}) {
			open $self->{'fh'}->{$category}, "bzcat $files | sort -k1,10 -t\$\'\t\' -s |" or die "can't open for $category\n";
		} else {
			open $self->{'fh'}->{$category}, "bzcat $files |" or die "can't open for $category\n";
		}
		$self->getNextRead($category);

	}
	return $self;
}

sub closeFile {
	my $self = shift;
	my $category = shift;
	close $self->{'fh'}->{$category};
	delete $self->{'fh'}->{$category};

}

sub closeFiles {
	my $self = shift;
	foreach my $fh (keys %{$self->{'fh'}}) {
		close $self->{'fh'}->{$fh};
		delete $self->{'fh'}->{$fh};
	}

}


sub getNextRead {
	my $self = shift;
	my $category = shift;
	my $headQR;
	unless (scalar @_) {
		$headQR = $self->{'headQR'};
	} else {
		$headQR = shift(@_);
	}
	my $line;
	my $read;
	while (1) {
		if (exists $self->{'prevLine'}->{$category}) {
			#print STDERR "Here in previous line\n";
			#print STDERR "LINE = $self->{'prevLine'}->{$category}\n";
			$line = $self->{'prevLine'}->{$category};
			$line =~ m/$headQR/;
			$read = $1;
			$self->{'currentHeader'}->{$category} = $read;
			#print STDERR "SET currentHeader to $read\nActual variable $self->{'currentHeader'}->{$category} \n";
			delete $self->{'currentLines'}->{$category};
			delete $self->{'prevLine'}->{$category};
		} elsif (eof($self->{'fh'}->{$category})) {
			$self->{'finished'}->{$category} = 1;
			last
		} else {
			$line = $self->{'fh'}->{$category}->getline;
			chomp $line;
			$line =~ m/$headQR/;
			$read = $1;
		}
		unless (exists $self->{'currentHeader'}->{$category}) {
			 $self->{'currentHeader'}->{$category} = $read;
		}
		if ($read ne $self->{'currentHeader'}->{$category}) {
			$self->{'prevLine'}->{$category} = $line;
			#print STDERR "set {'prevLine'}->{$category} to $line\n";
			#print STDERR "ACTUAL VARIABLE: $self->{'prevLine'}->{$category}\n";
			last;
		} else {
			push @{$self->{'currentLines'}->{$category}}, $line;
			redo;
		}
	}
	return($self);

}


return 1;
