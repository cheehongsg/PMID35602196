#!env perl

##
## Covid19 - Analysis
## Copyright (C) 2020-2021 Chee-Hong WONG (Dr Chia-Lin Wei Laboratory)
##

use strict;
use Data::Dumper;
use Getopt::Long;
use FindBin;
use lib $FindBin::Bin;
use JSON;

my $G_DEBUG=0;
my $G_MINMAPQ=0;
# my $G_MINMAPQ=255;
# my $G_SPLITPOOL_MINMAPQ=0;
my $G_SPLITPOOL_MINMAPQ=255;

# my $G_DEBUG=1;
my $G_ANNOTATION_PRIMERS = $FindBin::Bin.'/'.'../GT-Primers/primers.bed';
my $G_PRIMER_BUFFER_BP=2;

my $G_ANNOTATION_JUNCTIONS = $FindBin::Bin.'/'.'../annotation/sar-cov2-junctions.sorted.bed';
my $G_ANNOTATION_FEATURES = $FindBin::Bin.'/'.'../annotation/sar-cov2-features.bed';

my $G_SAMPLEMETAFILE = 'sampleMeta.xls';

# perl arctictools.pl <command>
my $G_USAGE = "
$0 <command> -h
<command> [splitreads,splitpool,splitpoolcov]
";

my $command = undef;
if (!defined $ARGV[0] || substr($ARGV[0],0,1) eq '-' || substr($ARGV[0],0,2) eq '--') {
	die("Please specify the command.\n",$G_USAGE);
}
$command = shift @ARGV; $command = lc $command;

# auto-flush for both stderr and stdout
select(STDERR);
$| = 1;
select(STDOUT);
$| = 1;

if ('splitreads' eq $command) {
	runSplitReads();
} elsif ('splitpool' eq $command) {
	runSplitPool();
} elsif ('splitpoolcov' eq $command) {
	runGenerateCoverage();
} elsif ('tabulatecov' eq $command) {
	runTabulateCoverage();
} elsif ('gridexpress' eq $command) {
	runGridExpression();
} else {
	print $G_USAGE;
}

exit 0;

sub loadArticPrimers {
    my ($primersRef) = @_;
    open INFILE, $G_ANNOTATION_PRIMERS || die "Fail to open $G_ANNOTATION_PRIMERS\n$!\n";
    $primersRef->{primerPairs} = ();
    $primersRef->{primers} = {};
    $primersRef->{lastPrimerId} = 0;
    $primersRef->{posToLPrimer} = ();
    $primersRef->{posToRPrimer} = ();
    $primersRef->{maxpos} = 0;
    $primersRef->{minpos} = 0;
    $primersRef->{maxposMarker} = 'E';

    while (<INFILE>) {
        chomp();
        my @bits = split(/\t/);
        my @nameBits = split(/\_/, $bits[3]);

        my %primer = (pid=>int($nameBits[1]), start=>$bits[1]+1, end=>$bits[2]);
        $primersRef->{lastPrimerId} = $primer{pid} if ($primer{pid}>$primersRef->{lastPrimerId});
        if ('LEFT' eq $nameBits[2]) {
            $primer{type} = 'L';
            $primer{id} = sprintf("%02dL", $primer{pid});
        } elsif ('RIGHT' eq $nameBits[2]) {
            $primer{type} = 'R';
            $primer{id} = sprintf("%02dR", $primer{pid});
        } else {
            die "Unknown id $bits[3]\n";
        }

        $primersRef->{primers}->{$primer{id}} = \%primer;
        if (!defined $primersRef->{primerPairs}->[$primer{pid}]) {
            $primersRef->{primerPairs}->[$primer{pid}] = {id=>$primer{pid}, L=>undef, R=>undef};
        }
        $primersRef->{primerPairs}->[$primer{pid}]->{$primer{type}} = \%primer;

        $primersRef->{maxpos} = $bits[2] if ($primersRef->{maxpos}<$bits[2]);
    }
    close INFILE;

    # set up the pos lookup
    # primers
    for(my $p=1; $p<=$primersRef->{lastPrimerId}; ++$p) {
        my $primerPairRef = $primersRef->{primerPairs}->[$p];
        my $primerLRef = $primerPairRef->{L};
        my $primerRRef = $primerPairRef->{R};
        for(my $i=$primerLRef->{start}; $i<=$primerLRef->{end}; ++$i) {
            $primersRef->{posToLPrimer}->[$i] = $primerLRef->{id};
        }
        for(my $i=1; $i<=$G_PRIMER_BUFFER_BP; ++$i) {
            $primersRef->{posToLPrimer}->[$primerLRef->{start}-$i] = sprintf("%s-%d", $primerLRef->{id}, $i);
        }
        for(my $i=$primerRRef->{start}; $i<=$primerRRef->{end}; ++$i) {
            $primersRef->{posToRPrimer}->[$i] = $primerRRef->{id};
        }
        for(my $i=1; $i<=$G_PRIMER_BUFFER_BP; ++$i) {
            $primersRef->{posToRPrimer}->[$primerRRef->{end}+$i] = sprintf("%s+%d", $primerRRef->{id}, $i);
        }
    }
    # forward
    for(my $p=1; $p<=$primersRef->{lastPrimerId}; ++$p) {
        my $primerPairRef = $primersRef->{primerPairs}->[$p];
        my $primerLRef = $primerPairRef->{L};
        my $id = sprintf("%02dM", $primerPairRef->{id});
        if ($p==$primersRef->{lastPrimerId}) {
            my $primerRRef = $primerPairRef->{R};
            for(my $i=$primerLRef->{end}+1; $i<$primerRRef->{start}-1; ++$i) {
                $primersRef->{posToLPrimer}->[$i] = $id;
            }
        } else {
            my $primerPairNextRef = $primersRef->{primerPairs}->[$p+1];
            my $primerLNextRef = $primerPairNextRef->{L};
            for(my $i=$primerLRef->{end}+1; $i<$primerLNextRef->{start}-$G_PRIMER_BUFFER_BP; ++$i) {
                $primersRef->{posToLPrimer}->[$i] = $id;
            }
        }
    }
    # reverse
    for(my $p=1; $p<=$primersRef->{lastPrimerId}; ++$p) {
        my $primerPairRef = $primersRef->{primerPairs}->[$p];
        my $primerRRef = $primerPairRef->{R};
        my $id = sprintf("%02dM", $primerPairRef->{id});
        if (1==$p) {
            my $primerLRef = $primerPairRef->{L};
            for(my $i=$primerLRef->{start}; $i<$primerRRef->{start}; ++$i) {
                $primersRef->{posToRPrimer}->[$i] = $id;
            }
        } else {
            my $primerPairPrevRef = $primersRef->{primerPairs}->[$p-1];
            my $primerRPrevRef = $primerPairPrevRef->{R};
            for(my $i=$primerRPrevRef->{end}+$G_PRIMER_BUFFER_BP+1; $i<$primerRRef->{start}; ++$i) {
                $primersRef->{posToRPrimer}->[$i] = $id;
            }
        }
    }

    # B=BEGIN, E=END
    if (defined $primersRef->{primerPairs}->[1]) {
        my $primerPairRef = $primersRef->{primerPairs}->[1];
        my $primerLRef = $primerPairRef->{L};
        $primersRef->{minpos} = $primerLRef->{start};
    }
}

sub locatePrimer {
    my ($exonsRef, $strand, $primersRef) = @_;
    # TODO: get the start and end position
    # TODO: determine the start and end position based on $strand
    if ('+' eq $strand) {
        # forward strand ; use the first, and maybe next exon?
        my $exonRef = $exonsRef->[0];
        if ($exonRef->{start}>$primersRef->{maxpos}) {
            return $primersRef->{maxposMarker} = 'E';
        } else {
            if ($exonRef->{start}<$primersRef->{minpos}) {
                if ($exonRef->{end}>$primersRef->{minpos}) {
                    my $retVal = $primersRef->{posToLPrimer}->[$primersRef->{minpos}];
                    return $retVal;
                }
            } else {
                my $retVal = $primersRef->{posToLPrimer}->[$exonRef->{start}];
                if (!defined $retVal || '' eq $retVal) {
                    $retVal = $primersRef->{posToRPrimer}->[$exonRef->{start}];
                    if (defined $retVal) {
                        $retVal = lc($retVal);
                    }
                }
                return $retVal;
            }
        }
    } else {
        # reverse strand ; use the last, and maybe before exon?
        my $exonRef = $exonsRef->[scalar(@{$exonsRef})-1];
        if ($exonRef->{end}>$primersRef->{maxpos}) {
            return $primersRef->{maxposMarker} = 'E';
        } else {
            # TODO: will start's special condition happen here?!
            my $retVal = $primersRef->{posToRPrimer}->[$exonRef->{end}];
            if (!defined $retVal || '' eq $retVal) {
                $retVal = $primersRef->{posToLPrimer}->[$exonRef->{end}];
                if (defined $retVal) {
                    $retVal = lc($retVal);
                }
            }
            return $retVal;
        }
    }
}

sub isPrimerChimeric {
    my ($exonsRef, $strand, $primersRef) = @_;
    return 0 if (scalar(@{$exonsRef})==0);

    # check if each exon is just a primer region
    if ('+' eq $strand) {
        my $prevPrimerPairRef = undef;
        foreach my $exonRef (@{$exonsRef}) {
            my @singleton = ($exonRef);
            my $zp = locatePrimer(\@singleton, $strand, $primersRef);
            my $primerRef = $primersRef->{primers}->{$zp};
            my $start = $primerRef->{start}-$G_PRIMER_BUFFER_BP;
            my $end = $primerRef->{end}+$G_PRIMER_BUFFER_BP;
            if ($start<=$exonRef->{start} && $exonRef->{end}<=$end) {
                if (defined $prevPrimerPairRef && $primerRef->{pid}!=$prevPrimerPairRef->{id} 
                    && $exonRef->{end}<=$prevPrimerPairRef->{'R'}->{end}) {
                    # still within the last fragment, although the next primer has been assigned
                } else {
                    return 1;
                }
            }
            $prevPrimerPairRef = $primersRef->{primerPairs}->[$primerRef->{pid}];
        }
    } else {
        my $prevPrimerPairRef = undef;
        foreach my $exonRef (reverse @{$exonsRef}) {
            my @singleton = ($exonRef);
            my $zp = locatePrimer(\@singleton, $strand, $primersRef);
            my $primerRef = $primersRef->{primers}->{$zp};
            my $start = $primerRef->{start}-$G_PRIMER_BUFFER_BP;
            my $end = $primerRef->{end}+$G_PRIMER_BUFFER_BP;
            if ($start<=$exonRef->{start} && $exonRef->{end}<=$end) {
                if (defined $prevPrimerPairRef && $primerRef->{pid}!=$prevPrimerPairRef->{id} 
                    && $prevPrimerPairRef->{'L'}->{start}<=$exonRef->{start}) {
                    # still within the last fragment, although the next primer has been assigned
                } else {
                    return 1;
                }
            }
            $prevPrimerPairRef = $primersRef->{primerPairs}->[$primerRef->{pid}];
        }
    }

    return 0;
}

sub replaceTag {
    my ($tag, $tagType, $tagValue, $tagsRef) = @_;

    my $found = 0;
    for(my $i=11; $i<scalar(@{$tagsRef}); ++$i) {
        my $existing = $tagsRef->[$i];
        if ($existing =~ /^$tag:$tagType:/) {
            $tagsRef->[$i] = sprintf("%s:%s:%s", $tag, $tagType, $tagValue);
            $found = 1;
        }
    }
    if (0==$found) {
        push @{$tagsRef}, sprintf("%s:%s:%s", $tag, $tagType, $tagValue);
    }
}

sub getEndSite {
    my ($startPos, $cigar) = @_;
    my @cigarOps = split(/([SHNDIM=X])/, $cigar);
    my $numOps = scalar(@cigarOps);
    my $refPos = $startPos-1;
    for(my $i=0; $i<$numOps; $i+=2) {
        my $span = $cigarOps[$i];
        my $op = $cigarOps[$i+1];
        if ('N' eq $op || 'D' eq $op || 'M' eq $op || '=' eq $op || 'X' eq $op) {
            $refPos += $span;
        }
    }
    return $refPos;
}

sub getMaxJump {
    my ($jumpsRef) = @_;
    my $seenMaxN = 0;
    foreach my $jumpRef (@{$jumpsRef}) {
        my $Ns = $jumpRef->{j3} - $jumpRef->{j5} - 1;
        $seenMaxN = $Ns if ($Ns>$seenMaxN);
    }
    return $seenMaxN;
}

# TODO: sync with getExonSites for condition
sub getSwitchingSites {
    my ($startPos, $cigar, $jumpsRef) = @_;
    @{$jumpsRef} = ();
    my @cigarOps = split(/([SHNDIM=X])/, $cigar);
    my $numOps = scalar(@cigarOps);
    my $refPos = $startPos-1;
    for(my $i=0; $i<$numOps; $i+=2) {
        my $span = $cigarOps[$i];
        my $op = $cigarOps[$i+1];
        # print $refPos;
        if ('N' eq $op) {
            push @{$jumpsRef}, {j5=>$refPos, j3=>$refPos+$span};
        } elsif ('D' eq $op && $span>=20) {
            push @{$jumpsRef}, {j5=>$refPos, j3=>$refPos+$span};
        }

        if ('N' eq $op || 'D' eq $op || 'M' eq $op || '=' eq $op || 'X' eq $op) {
            $refPos += $span;
        }
        # print " ",$refPos," ", $span,$op, "\n";
    }
}

sub getExonSites {
    my ($startPos, $cigar, $exonsRef) = @_;
    @{$exonsRef} = ();
    my @cigarOps = split(/([SHNDIM=X])/, $cigar);
    my $numOps = scalar(@cigarOps);

    # let's work around the PacBio IsoSeq alignment bug
    # 82S24=26403N319=34D370=31D358=42N2D1754=
    my @adjusted = ();
    my $adjusted = 0;
    my $prevOp = '';
    my $prevSpan = 0;
    for(my $i=0; $i<$numOps; $i+=2) {
        my $span = $cigarOps[$i];
        my $op = $cigarOps[$i+1];

        if ('N' eq $op || 'D' eq $op) {
            if ('N' eq $prevOp || 'D' eq $prevOp) {
                # collate
                # $prevOp kept as-is
                $prevSpan += $span;
                $adjusted++;
            } else {
                # flush-and-reset
                push @adjusted, $prevSpan, $prevOp if ('' ne $prevSpan);
                $prevSpan = $span;
                $prevOp = $op;
            }
        } else {
            # flush-and-reset
            push @adjusted, $prevSpan, $prevOp if ('' ne $prevSpan);
            $prevSpan = $span;
            $prevOp = $op;
        }
    }
    # flush-and-reset
    push @adjusted, $prevSpan, $prevOp if ('' ne $prevSpan);
    if ($adjusted>0) {
        # collatable N/D in PacBio IsoSeq alignment
        @cigarOps = ();
        push @cigarOps, @adjusted;
        $numOps = scalar(@cigarOps);
    }

    my $refPos = $startPos-1;
    my $start = $refPos;
    my $end = -1;
    for(my $i=0; $i<$numOps; $i+=2) {
        my $span = $cigarOps[$i];
        my $op = $cigarOps[$i+1];
        # print $refPos;
        if ('N' eq $op) {
            $end = $refPos;
            push @{$exonsRef}, {start=>$start+1, end=>$end};
            $start = $refPos+$span;
        } elsif ('D' eq $op && $span>=20) {
            $end = $refPos;
            push @{$exonsRef}, {start=>$start+1, end=>$end};
            $start = $refPos+$span;
        }

        if ('N' eq $op || 'D' eq $op || 'M' eq $op || '=' eq $op || 'X' eq $op) {
            $refPos += $span;
        }
    }
    if ($refPos>$start) {
        $end = $refPos;
        push @{$exonsRef}, {start=>$start+1, end=>$end};
    }
}

sub runSplitReads {
    my $bamFile = $ARGV[0];

    my %primers = (); loadArticPrimers(\%primers);

    # let's work thru' the bam file
    # - generate .proper.bam (encode Lxx-Rxx, Lxx-Mxx, Mxx-Rxx, Mxx-Mxx)
    # - generate .improper.bam
    my $outFile = $bamFile; $outFile=~s/\.bam$//; $outFile.='.amplicontag.bam';
    open INFILE, "samtools view -H $bamFile |" || die "Fail to open $bamFile\n$!\n";
    open OUTFILE, "| samtools view -Sb > $outFile" || die "Fail to open $outFile\n$!\n";
    while (<INFILE>) {
        print OUTFILE $_;
    }
    close INFILE;

    # TODO: decide if we want to have mapq score here!
    # open INFILE, "samtools view -q $G_MINMAPQ $bamFile |" || die "Fail to open $bamFile\n$!\n";
    open INFILE, "samtools view $bamFile |" || die "Fail to open $bamFile\n$!\n";
    printf STDERR "Processing $bamFile..";
    my $block = 10000; my $trigger = $block; my $numRecs = 0;
    while (<INFILE>) {
        chomp();
        my @bits = split(/\t/);

        my $isSequencingArtefact = 0;
        for(my $i=11; $i<scalar(@bits); ++$i) {
            if ($bits[$i] eq 'DT:Z:SQ') {
                $isSequencingArtefact = 1;
                last;
            }
        }
        next if ($isSequencingArtefact == 1);

        # let's find the maximum Ns in this alignment
        # this maxN must be within the specified bounds
        my $end = getEndSite($bits[3], $bits[5]);
        replaceTag('ZE', 'i', $end, \@bits);

        # let's find the maximum Ns in this alignment
        # this maxN must be within the specified bounds
        my @jumps = ();
        getSwitchingSites($bits[3], $bits[5], \@jumps);
        my $seenMaxN = getMaxJump(\@jumps);
        replaceTag('ZN', 'i', $seenMaxN, \@bits);
        my @values = ();
        grep { push @values, sprintf("%d-%d", $_->{j5}, $_->{j3});}  (@jumps);
        replaceTag('ZJ', 'Z', join(",", @values), \@bits);
        replaceTag('ZK', 'Z', scalar(@values).','.$seenMaxN, \@bits);

        # TODO: ZT:Z:proper/improper
        my $zt = '';
        my $flag = int($bits[1]);
        if (2==($flag & 2)) {
            $zt = 'P';
        } else {
            $zt = 'I';
        }

        # TODO: ZP:Z:Lxx,Rxx,Mxx
        # TODO: only handle both reads of read-pair mapped?
        my $zp = 'U';
        if (4==($flag & 4)) {
            # read is unmapped
        } else {
            if (256==($flag & 256)) {
                # not primary alignment
                $zp = 'NPA';
            } elsif (512==($flag & 512)) {
                # read fails platform/vendor quality checks
                $zp = 'QUAL';
            } elsif (1024==($flag & 1024)) {
                # read is PCR or optical duplicate
                $zp = 'DUP';
            } elsif (2048==($flag & 2048)) {
                # supplementary alignment
                $zp = 'SA';
            } else {
                my @exons = ();
                getExonSites($bits[3], $bits[5], \@exons);
                if (16==($flag & 16)) {
                    my $isPrimerChimeric = isPrimerChimeric(\@exons, '-', \%primers);
                    if (0==$isPrimerChimeric) {
                        # read reverse strand, needs the end
                        $zp = locatePrimer(\@exons, '-', \%primers);
                        if ($zp eq '') {
                            print STDERR "DEBUG\n";
                        }
                    } else {
                        $zp = 'PCD';
                    }
                } else {
                    my $isPrimerChimeric = isPrimerChimeric(\@exons, '+', \%primers);
                    # read forwards strand, needs the start
                    if (0==$isPrimerChimeric) {
                        $zp = locatePrimer(\@exons, '+', \%primers);
                        if ($zp eq '') {
                            print STDERR "DEBUG\n";
                        }
                    } else {
                        $zp = 'PCD';
                    }
                }
                # TODO: attempt to conclude the properness of pairing?
                #       or can we wait for second pass?
                if (64==($flag & 64)) {
                    # first in pair
                } elsif (128==($flag & 128)) {
                    # second in pair
                } else {
                    # TODO: not paired-end?
                }
            }
        }
        replaceTag('ZP', 'Z', $zp, \@bits);
        replaceTag('ZT', 'Z', $zt, \@bits);

        # TODO: ZQ:Z:Lxx-Rxx, Lxx-Mxx, Mxx-Rxx, Mxx-Mxx

        # TODO: determine the primer information

        # TODO: inject the primer information

        # write to the appropriate bam file!
        print OUTFILE join("\t", @bits), "\n";
        # print join("\t", @bits), "\n" if ($zp =~ /M$/ || $zp eq '');
        print join("\t", @bits), "\n" if ($zp eq '');

        $numRecs++;
        if ($numRecs>=$trigger) {
            $trigger += $block;
            printf STDERR " %.2fM..", $numRecs / 1e6;
        }
    }
    close OUTFILE;
    close INFILE;
    printf STDERR " %d.. done!\n", $numRecs;
    # let's index the file
    my @commands = ('samtools', 'index', $outFile);
    my $commands = join(' ', @commands);
    system($commands);
}

sub getReadPairOutType {
    my ($qname, $readPairsRef) = @_;

    return 'ignore' if (!exists $readPairsRef->{$qname});

    my $readPairRef = $readPairsRef->{$qname};
    # 1) single read only
    return 'ignore' if (!exists $readPairRef->{1});
    return 'ignore' if (!exists $readPairRef->{2});

    # 2) does not contain marker 'M' in any of the read of pair
    my $value = lc(substr($readPairRef->{1}->{ZP}, -1));
    return 'ignore' if ('m' eq $value);
    $value = lc(substr($readPairRef->{2}->{ZP}, -1));
    return 'ignore' if ('m' eq $value);

    # 3) R/1 and R/2 must be of different strand
    return 'ignore' if ($readPairRef->{1}->{strand} eq $readPairRef->{2}->{strand});

    # 4) the primer-pair has to be of exclusively odd or even pool!
    return 'pcd' if ('PCD' eq $readPairRef->{1}->{ZP} || 'PCD' eq $readPairRef->{2}->{ZP});

    my $r1odd = (1==(int(substr($readPairRef->{1}->{ZP}, 0, -1)) % 2)) ? 1 : 0;
    my $r2odd = (1==(int(substr($readPairRef->{2}->{ZP}, 0, -1)) % 2)) ? 1 : 0;

    # 5) must be in-ward pointing orientation
    if ('+' eq $readPairRef->{1}->{strand}) {
        # if ($readPairRef->{2}->{end}<$readPairRef->{1}->{start}) {
        if ($readPairRef->{2}->{start}<$readPairRef->{1}->{start} || $readPairRef->{1}->{end}>$readPairRef->{2}->{end}) {
            return 'ignore';
        }
    } else {
        # if ($readPairRef->{1}->{end}<$readPairRef->{2}->{start}) {
        if ($readPairRef->{1}->{start}<$readPairRef->{2}->{start} || $readPairRef->{2}->{end}>$readPairRef->{1}->{end}) {
            return 'ignore';
        }
    }
    if ($r1odd != $r2odd) {
        return 'mixed';
    } elsif (1==$r1odd) {
        return 'odd';
    } else {
        return 'even';
    }
}

sub runSplitPool {
    my $bamFile = $ARGV[0];

    # TODO: let's count the mapq filter separately!

    # let's summarize each read-pair so that we can gather statistics
    # open INFILE, "samtools view -q $G_MINMAPQ $bamFile |" || die "Fail to open $bamFile\n$!\n";
    open INFILE, "samtools view $bamFile |" || die "Fail to open $bamFile\n$!\n";
    printf STDERR "Processing $bamFile..";
    my %readPairs = ();
    my $block = 10000; my $trigger = $block; my $numRecs = 0;
    my $mapqFailed = 0;
    while (<INFILE>) {
        chomp();
        my @bits = split(/\t/);

        my $isSequencingArtefact = 0;
        for(my $i=11; $i<scalar(@bits); ++$i) {
            if ($bits[$i] eq 'DT:Z:SQ') {
                $isSequencingArtefact = 1;
                last;
            }
        }
        next if ($isSequencingArtefact == 1);

        # we skip some classes of read
        my $flag = int($bits[1]);
        next if (4==($flag & 4)); # read is unmapped
        next if (256==($flag & 256)); # not primary alignment
        next if (512==($flag & 512)); # read fails platform/vendor quality checks
        next if (1024==($flag & 1024)); # read is PCR or optical duplicate
        next if (2048==($flag & 2048)); # # supplementary alignment

        # let's get the tags and gather the statistics
        # ZE, ZN, ZJ, ZK
        # ZT
        # ZP
        my %tags = ();
        for(my $i=11; $i<scalar(@bits); ++$i) {
            my ($tag, $tagType, $tagValue) = split(/\:/, $bits[$i], 3);
            $tags{$tag} = $tagValue;
        }

        my $qname = $bits[0];
        if (!exists $readPairs{$qname}) {
            # '1'=>{}, '2'=>{}
            my %item = ();
            $readPairs{$qname} = \%item;
        }

        my $readPairRef = $readPairs{$qname};
        $readPairRef->{ZT} = $tags{ZT};
        my $read1or2 = 0;
        if (64==($flag & 64)) {
            $read1or2 = 1;
        } elsif (128==($flag & 128)) {
            $read1or2 = 2;
        } else {
            die "ERROR: Unknown read (1 or 2?)\n$_\n";
        }
        if (!exists $readPairs{$qname}->{$read1or2}) {
            my %item = ();
            $readPairs{$qname}->{$read1or2} = \%item;
        }
        my $readPairReadRef = $readPairs{$qname}->{$read1or2};
        my $strand = (16==($flag & 16)) ? '-' : '+';
        $readPairReadRef->{strand} = $strand;
        $readPairReadRef->{ZP} = $tags{ZP};
        $readPairReadRef->{start} = int($bits[3]);
        $readPairReadRef->{end} = int($tags{ZE});
        $readPairReadRef->{mapq} = int($bits[4]);

        if ($bits[4]<$G_SPLITPOOL_MINMAPQ) {
            $mapqFailed++;
        }

        $numRecs++;
        if ($numRecs>=$trigger) {
            $trigger += $block;
            printf STDERR " %.2fM..", $numRecs / 1e6;
        }
    }
    close INFILE;
    printf STDERR " %d.. (%d mapq<%d) done!\n", $numRecs, $mapqFailed, $G_SPLITPOOL_MINMAPQ;

    # output the read-pairs information
    my $outFile = $bamFile; $outFile=~s/\.bam$//; $outFile.='.xls.gz';
    printf STDERR "Writing $outFile..";
    open OUTFILE, "| gzip -c - > $outFile" || die "Fail to open $outFile\n$!\n";
    print OUTFILE '#', $bamFile, "\n";
    my @cols = ('readId', 'presence', 'R1_ZP', 'R2_ZP', 'R1_strand', 'R2_strand', 'R1_mapq', 'R2_mapq');
    print OUTFILE '#', join("\t", @cols), "\n";
    my @ids = sort keys %readPairs;
    foreach my $id (@ids) {
        my $readPairRef = $readPairs{$id};
        @cols = ($id);
        my $presence = 0;
        $presence |= 1 if (exists $readPairRef->{1});
        $presence |= 2 if (exists $readPairRef->{2});
        push @cols, $presence;
        push @cols, (exists $readPairRef->{1}) ? $readPairRef->{1}->{ZP} : '.';
        push @cols, (exists $readPairRef->{2}) ? $readPairRef->{2}->{ZP} : '.';
        push @cols, (exists $readPairRef->{1}) ? $readPairRef->{1}->{strand} : '.';
        push @cols, (exists $readPairRef->{2}) ? $readPairRef->{2}->{strand} : '.';
        push @cols, (exists $readPairRef->{1}) ? $readPairRef->{1}->{mapq} : '.';
        push @cols, (exists $readPairRef->{2}) ? $readPairRef->{2}->{mapq} : '.';

        print OUTFILE join("\t", @cols), "\n";
    }
    close OUTFILE;
    printf STDERR " done!\n";

    # output the cleaned read-pairs
    # we are going to exclude:
    # 1) single read only
    # 2) does not contain marker 'M' in any of the read of pair
    # 3) R/1 and R/2 must be of different strand
    # 4) the primer-pair has to be of exclusively odd or even pool!
    # the output groups:
    # A) <prefix>.oddpool.bam
    # B) <prefix>.evenpool.bam
    # C) <prefix>.mixedpool.bam

    my $outPrefix = $bamFile; $outPrefix=~s/\.bam$//;
    my $oddFile = sprintf("%s.oddpool.bam", $outPrefix);
    my $evenFile = sprintf("%s.evenpool.bam", $outPrefix);
    my $mixedFile = sprintf("%s.mixedpool.bam", $outPrefix);
    my $pcdFile = sprintf("%s.primerchimeric.bam", $outPrefix);
    open INFILE, "samtools view -H $bamFile |" || die "Fail to open $bamFile\n$!\n";
    open OUTOFILE, "| samtools view -Sb > $oddFile" || die "Fail to open $oddFile\n$!\n";
    open OUTEFILE, "| samtools view -Sb > $evenFile" || die "Fail to open $evenFile\n$!\n";
    open OUTMFILE, "| samtools view -Sb > $mixedFile" || die "Fail to open $mixedFile\n$!\n";
    open OUTPFILE, "| samtools view -Sb > $pcdFile" || die "Fail to open $pcdFile\n$!\n";
    while (<INFILE>) {
        print OUTOFILE $_;
        print OUTEFILE $_;
        print OUTMFILE $_;
        print OUTPFILE $_;
    }
    close INFILE;

    open INFILE, "samtools view -q $G_MINMAPQ $bamFile |" || die "Fail to open $bamFile\n$!\n";
    printf STDERR "Processing $bamFile..";
    my $block = 10000; my $trigger = $block; my $numRecs = 0;
    my $oddCount = 0; my $evenCount = 0; my $mixedCount = 0; my $pcdCount = 0; my $ignoreCount = 0;
    my $oddMapqCount = 0; my $evenMapqCount = 0; my $mixedMapqCount = 0; my $pcdMapqCount = 0; my $ignoreMapqCount = 0;
    while (<INFILE>) {
        chomp();
        my @bits = split(/\t/);

        my $isSequencingArtefact = 0;
        for(my $i=11; $i<scalar(@bits); ++$i) {
            if ($bits[$i] eq 'DT:Z:SQ') {
                $isSequencingArtefact = 1;
                last;
            }
        }
        next if ($isSequencingArtefact == 1);

        # we skip some classes of read
        my $flag = int($bits[1]);
        next if (4==($flag & 4)); # read is unmapped
        next if (256==($flag & 256)); # not primary alignment
        next if (512==($flag & 512)); # read fails platform/vendor quality checks
        next if (1024==($flag & 1024)); # read is PCR or optical duplicate
        next if (2048==($flag & 2048)); # # supplementary alignment

        my $outType = getReadPairOutType($bits[0], \%readPairs);
        if ('odd' eq $outType) {
            $oddCount++;
            $oddMapqCount++ if ($bits[4]<$G_SPLITPOOL_MINMAPQ);
            print OUTOFILE join("\t", @bits), "\n";
        } elsif ('even' eq $outType) {
            $evenCount++;
            $evenMapqCount++ if ($bits[4]<$G_SPLITPOOL_MINMAPQ);
            print OUTEFILE join("\t", @bits), "\n";
        } elsif ('mixed' eq $outType) {
            $mixedCount++;
            $mixedMapqCount++ if ($bits[4]<$G_SPLITPOOL_MINMAPQ);
            print OUTMFILE join("\t", @bits), "\n";
        } elsif ('pcd' eq $outType) {
            $pcdCount++;
            $pcdMapqCount++ if ($bits[4]<$G_SPLITPOOL_MINMAPQ);
            print OUTPFILE join("\t", @bits), "\n";
        } else {
            $ignoreCount++;
            $ignoreMapqCount++ if ($bits[4]<$G_SPLITPOOL_MINMAPQ);
            # ignore!
        }

        $numRecs++;
        if ($numRecs>=$trigger) {
            $trigger += $block;
            printf STDERR " %.2fM..", $numRecs / 1e6;
        }
    }
    close OUTOFILE;
    close OUTEFILE;
    close OUTMFILE;
    close INFILE;
    printf STDERR " %d.. odd(%d,%d mapq<%d), even(%d,%d mapq<%d), mixed(%d,%d mapq<%d), pcd(%d,%d mapq<%d), ignore(%d,%d mapq<%d) done!\n", $numRecs, 
        $oddCount, $oddMapqCount, $G_SPLITPOOL_MINMAPQ,
        $evenCount, $evenMapqCount, $G_SPLITPOOL_MINMAPQ,
        $mixedCount, $mixedMapqCount, $G_SPLITPOOL_MINMAPQ,
        $pcdCount, $pcdMapqCount, $G_SPLITPOOL_MINMAPQ,
        $ignoreCount, $ignoreMapqCount, $G_SPLITPOOL_MINMAPQ;

    # let's index the file
    foreach my $file ($oddFile, $evenFile, $mixedFile, $pcdFile) {
        printf STDERR "Indexing %s..", $file;
        my @commands = ('samtools', 'index', $file);
        my $commands = join(' ', @commands);
        system($commands);
        printf STDERR " done!\n";
    }
}

sub runGenerateCoverage {
    my $bamFile = $ARGV[0];

    # let's generate the coverage for the splitted pool
    my $outPrefix = $bamFile; $outPrefix=~s/\.bam$//;
    my $oddFile = sprintf("%s.oddpool.bam", $outPrefix);
    my $evenFile = sprintf("%s.evenpool.bam", $outPrefix);
    my $mixedFile = sprintf("%s.mixedpool.bam", $outPrefix);
    foreach my $file ($oddFile, $evenFile, $mixedFile) {
        open INFILE, "samtools view -q $G_MINMAPQ $file |" || die "Fail to open $file\n$!\n";
        printf STDERR "Processing $file..";

        my $outfile = $file; $outfile =~ s/\.bam$//; $outfile.= '.bdg';
        open OUTFILE, ">$outfile" || die "Fail to open $outfile\n$!\n";

        my %readPairs = ();
        my $block = 10000; my $trigger = $block; my $numRecs = 0;
        while (<INFILE>) {
            chomp();
            my @bits = split(/\t/);

            # we skip some classes of read
            my $flag = int($bits[1]);
            next if (4==($flag & 4)); # read is unmapped
            next if (256==($flag & 256)); # not primary alignment
            next if (512==($flag & 512)); # read fails platform/vendor quality checks
            next if (1024==($flag & 1024)); # read is PCR or optical duplicate
            next if (2048==($flag & 2048)); # # supplementary alignment

            my $qname = $bits[0];
            if (!exists $readPairs{$qname}) {
                # '1'=>{}, '2'=>{}
                my %item = ();
                $readPairs{$qname} = \%item;
            }

            my $readPairRef = $readPairs{$qname};
            my $read1or2 = 0;
            if (64==($flag & 64)) {
                $read1or2 = 1;
            } elsif (128==($flag & 128)) {
                $read1or2 = 2;
            } else {
                die "ERROR: Unknown read (1 or 2?)\n$_\n";
            }
            if (!exists $readPairs{$qname}->{$read1or2}) {
                my %item = ();
                $readPairs{$qname}->{$read1or2} = \%item;
            }
            my $readPairReadRef = $readPairs{$qname}->{$read1or2};
            my @exons = ();
            # keep both exons set and work on it
            getExonSites($bits[3], $bits[5], \@exons);
            $readPairReadRef->{exons} = \@exons;

            $numRecs++;
            if ($numRecs>=$trigger) {
                $trigger += $block;
                printf STDERR " %.2fM..", $numRecs / 1e6;
            }
        }
        close INFILE;
        printf STDERR " %d.. done!\n", $numRecs;

        # set up the genome backbone
        open INFILE, "samtools view -H $file |" || die "Fail to open $file\n$!\n";
        my $genomeName = 'unnamed';
        my $genomeSize = -1;
        while (<INFILE>) {
            next if (!/^\@SQ/);
            chomp();
            die "Does not support multiple reference sequences!\n" if ($genomeSize!=-1);
            my @bits = split(/\t/);
            my @nameBits = split(/\:/, $bits[1]);
            $genomeName = $nameBits[1];
            my @lenBits = split(/\:/, $bits[2]);
            $genomeSize = int($lenBits[1]);
        }
        close INFILE;

        # update the count
        printf STDERR "Accumulating reads..";
        $block = 10000; $trigger = $block; $numRecs = 0;
        my @genomeReadPairCounts = (0) x ($genomeSize+1);
        while (my ($rpid, $readPairRef) = each %readPairs) {
            next if (!exists $readPairRef->{1});
            next if (!exists $readPairRef->{2});
            my %coverage = ();
            foreach my $readRef ($readPairRef->{1}, $readPairRef->{2}) {
                foreach my $exonRef (@{$readRef->{exons}}) {
                    for(my $i=$exonRef->{start}; $i<=$exonRef->{end}; ++$i) {
                        $coverage{$i}++;
                    }
                }
            }
            foreach my $gpos (keys %coverage) {
                $genomeReadPairCounts[$gpos]++;
            }

            $numRecs++;
            if ($numRecs>=$trigger) {
                $trigger += $block;
                printf STDERR " %.2fM..", $numRecs / 1e6;
            }
        }
        printf STDERR " %d.. done!\n", $numRecs;

        # write the bedgraph
        printf STDERR "Processing $outfile..";
        for(my $i=1; $i<=$genomeSize; ++$i) {
            next if (0==$genomeReadPairCounts[$i]);
            my @cols = ($genomeName, $i-1, $i, $genomeReadPairCounts[$i]);
            print OUTFILE join("\t", @cols), "\n";
        }
        close OUTFILE;
        printf STDERR " done!\n";
    }
}

sub loadSampleInfo {
  my ($file, $rowsRef) = @_;
  open INFILE, $file || die "Fail to open $file\n$!\n";
  my $headers = <INFILE>; chomp($headers);
  my @headers = split(/\t/, $headers);
  my $numHeaders = scalar(@headers);
  while (<INFILE>) {
    chomp();
    my @bits = split(/\t/);
    my %item = ();
    for(my $i=0; $i<$numHeaders; ++$i) {
      $item{$headers[$i]} = $bits[$i];
    }
    my @parts = split(/\./, $item{Sample});
    $item{id} = $parts[0];
    $rowsRef->{$item{id}} = \%item;
  }
  close INFILE;
}

sub loadPoolCoverage {
    my ($bdgFile, $covsRef) = @_;
    @{$covsRef} = ();
    open INFILE, $bdgFile || die "Fail to open $bdgFile\n$!\n";
    while (<INFILE>) {
        chomp();
        my @bits = split(/\t/);
        $covsRef->[int($bits[2])] = int($bits[3]);
    }
    close INFILE;
}

sub tabulatePoolCoverage {
    my ($oddBdg, $evenBdg, $primersRef, $rowsRef) = @_;

    # to calculate max
    # average and median does not make sense for artic
    my @cov = (); loadPoolCoverage($oddBdg, \@cov);
    for(my $i=1; $i<=$primersRef->{lastPrimerId}; $i+=2) {
        my $primerPairRef = $primersRef->{primerPairs}->[$i];
        my $max = 0; my $total = 0;
        my @counts = ();
        for(my $j=$primerPairRef->{L}->{start}; $j<=$primerPairRef->{R}->{end}; ++$j) {
            $max = $cov[$j] if ($cov[$j]>$max);
            $total += $cov[$j];
            push @counts, $cov[$j];
        }
        my $span = $primerPairRef->{R}->{end} - $primerPairRef->{L}->{start} + 1;
        my $average = $total / $span;
        @counts = sort { int($a)<=>int($b) } @counts;
        my $median = $counts[int($span/2)];
        if (0==($span % 2)) {
            $median += $counts[int($span/2)-1];
            $median /= 2;
        }
        $rowsRef->{ampliconCov} = { lastPrimerId=>$primersRef->{lastPrimerId}, coverage=>[]} if (!exists $rowsRef->{ampliconCov});
        $rowsRef->{ampliconCov}->{coverage}->[$i] = {id=>$i, max=>$max, total=>$total, span=>$span, median=>$median};
    }
    @cov = (); loadPoolCoverage($evenBdg, \@cov);
    for(my $i=2; $i<=$primersRef->{lastPrimerId}; $i+=2) {
        my $primerPairRef = $primersRef->{primerPairs}->[$i];
        my $max = 0; my $total = 0;
        my @counts = ();
        for(my $j=$primerPairRef->{L}->{start}; $j<=$primerPairRef->{R}->{end}; ++$j) {
            $max = $cov[$j] if ($cov[$j]>$max);
            $total += $cov[$j];
            push @counts, $cov[$j];
        }
        my $span = $primerPairRef->{R}->{end} - $primerPairRef->{L}->{start} + 1;
        my $average = $total / $span;
        @counts = sort { int($a)<=>int($b) } @counts;
        my $median = $counts[int($span/2)];
        if (0==($span % 2)) {
            $median += $counts[int($span/2)-1];
            $median /= 2;
        }
        $rowsRef->{ampliconCov} = { lastPrimerId=>$primersRef->{lastPrimerId}, coverage=>[]} if (!exists $rowsRef->{ampliconCov});
        $rowsRef->{ampliconCov}->{coverage}->[$i] = {id=>$i, max=>$max, total=>$total, span=>$span, median=>$median};
    }
}

sub tabulateCountAndJunction {
    my ($bamFile, $primersRef, $rowsRef, $prefixLabel) = @_;

    # let's work thru' the bam file
    # 1) what's the total count
    # 2) what are the junctions and what's the background?
    $rowsRef->{jumps} = {} if (!exists $rowsRef->{jumps});
    # TODO: do we filter!!!
    open INFILE, "samtools view -q $G_MINMAPQ $bamFile |" || die "Fail to open $bamFile\n$!\n";
    printf STDERR "Processing $bamFile..";
    my $block = 10000; my $trigger = $block; my $numRecs = 0;
    my %ids = ();
    while (<INFILE>) {
        chomp();
        my @bits = split(/\t/);

        # TODO: ZT:Z:proper/improper
        my $zt = '';
        my $flag = int($bits[1]);
        if (2==($flag & 2)) {
            $zt = 'P';
        } else {
            $zt = 'I';
        }
        die Dumper(\@bits) if ('P' ne $zt);

        # let's find the maximum Ns in this alignment
        # this maxN must be within the specified bounds
        my @jumps = ();
        getSwitchingSites($bits[3], $bits[5], \@jumps);

        my $zp = 'U';
        if (4==($flag & 4)) {
            # read is unmapped
        } else {
            if (256==($flag & 256)) {
                # not primary alignment
                $zp = 'NPA';
            } elsif (512==($flag & 512)) {
                # read fails platform/vendor quality checks
                $zp = 'QUAL';
            } elsif (1024==($flag & 1024)) {
                # read is PCR or optical duplicate
                $zp = 'DUP';
            } elsif (2048==($flag & 2048)) {
                # supplementary alignment
                $zp = 'SA';
            } else {
                # TODO: have to process each junction
                # my @exons = ();
                # getExonSites($bits[3], $bits[5], \@exons);
                # (16==($flag & 16))
                # $zp = locatePrimer(\@exons, '-', \%primers);
                # $zp = locatePrimer(\@exons, '+', \%primers);

                foreach my $jumpRef (@jumps) {
                    my $key = sprintf('%d-%d', $jumpRef->{j5}, $jumpRef->{j3});
                    if (!exists $rowsRef->{jumps}->{$key}) {
                        # TODO: report on the amplicon id for j5 and j3?
                        my @exons = ({start=>$jumpRef->{j5}-1, end=>$jumpRef->{j5}}, {start=>$jumpRef->{j3}-1, end=>$jumpRef->{j3}});
                        my $zpJ5 = locatePrimer(\@exons, '+', $primersRef);
                        my $zpJ3 = locatePrimer(\@exons, '-', $primersRef);
                        $rowsRef->{jumps}->{$key} = {j5=>$jumpRef->{j5}, j3=>$jumpRef->{j3}, j5AmpId=>$zpJ5, j3AmpId=>$zpJ3, reads=>{}};
                    }
                    $rowsRef->{jumps}->{$key}->{reads}->{$bits[0]}++;
                }
            }
        }

        if ($zp ne '') {
            $ids{$bits[0]}++;
        }

        $numRecs++;
        if ($numRecs>=$trigger) {
            $trigger += $block;
            printf STDERR " %.2fM..", $numRecs / 1e6;
        }
    }
    close OUTFILE;
    close INFILE;
    printf STDERR " %d.. done!\n", $numRecs;

    $rowsRef->{$prefixLabel.'rp'} = scalar(keys %ids);
}

sub runTabulateCoverage {
    my $bamFile = $ARGV[0];

    # primers
    my %primers = (); loadArticPrimers(\%primers);

    my %readPairs = ();
    my $outPrefix = $bamFile; $outPrefix=~s/\.bam$//;

    # 1) we will tally the max, median, average base coverage for each amplicon
    my $oddBdgFile = sprintf("%s.oddpool.bdg", $outPrefix);
    my $evenBdgFile = sprintf("%s.evenpool.bdg", $outPrefix);
    tabulatePoolCoverage($oddBdgFile, $evenBdgFile, \%primers, \%readPairs);

    # TODO: count the junction too?!
    # 2) we are going to count the number of read-pairs in odd & even pool
    # odd pool
    if (0==0) {
        my $file = sprintf("%s.oddpool.bam", $outPrefix);
        tabulateCountAndJunction ($file, \%primers, \%readPairs, 'odd-');
    }
    # even pool
    if (0==0) {
        my $file = sprintf("%s.evenpool.bam", $outPrefix);
        tabulateCountAndJunction ($file, \%primers, \%readPairs, 'even-');
    }
    $readPairs{'total-rp'} = $readPairs{'odd-rp'} + $readPairs{'even-rp'};

    # reporting!
    if (0==0) {
        my $numTotalRP = $readPairs{'total-rp'};
        if (0==$numTotalRP) {
            print "WARNING: total read-pair has zero entry, resetting to 1 for division\n";
            $numTotalRP = 1;
        }
        my $outfile = sprintf("%s.coverage.statistics", $outPrefix);
        open OUTFILE, ">$outfile" || die "Fail to open $outfile\n$!\n";
        my @cols = ('#', 'data_file', $bamFile);
        print OUTFILE join("\t", @cols), "\n";
        @cols = ('#', 'total_read_pair', $readPairs{'total-rp'});
        print OUTFILE join("\t", @cols), "\n";
        @cols = ('#', 'odd_read_pair', $readPairs{'odd-rp'}, sprintf("%.6f", $readPairs{'odd-rp'} / $numTotalRP));
        print OUTFILE join("\t", @cols), "\n";
        @cols = ('#', 'even_read_pair', $readPairs{'even-rp'}, sprintf("%.6f", $readPairs{'even-rp'} / $numTotalRP));
        print OUTFILE join("\t", @cols), "\n";

        print OUTFILE "#\n";
        @cols = ('#amplicon_id', 'span', 'max', 'median', 'total', 'average');
        print OUTFILE join("\t", @cols), "\n";

        my $coveragesRef = $readPairs{'ampliconCov'}->{'coverage'};
        for(my $i=1; $i<=$readPairs{'ampliconCov'}->{'lastPrimerId'}; ++$i) {
            my $coverageRef = $coveragesRef->[$i];
            @cols = ();
            grep { push @cols, (exists $coverageRef->{$_} && '' ne $coverageRef->{$_}) ? $coverageRef->{$_} : 0; } ('id', 'span', 'max', 'median', 'total');
            push @cols, sprintf("%.2f", $coverageRef->{'total'} / $coverageRef->{'span'});
            print OUTFILE join("\t", @cols), "\n";
        }
        close OUTFILE;

        my @jumps = sort { $a->{j5}<=>$b->{j5} || $a->{j3}<=>$b->{j3} } values %{$readPairs{jumps}};
        my $outfile = sprintf("%s.coverage.junctions", $outPrefix);
        open OUTFILE, ">$outfile" || die "Fail to open $outfile\n$!\n";
        my $G_SEQDEPTH = 1000000; # one million
        # my $G_SEQDEPTH = 100000; # hundred thousand
        # my $G_AMPLICONDEPTH = 1000; # one thousand
        # my $G_AMPLICONDEPTH = 100000; # hundred thousand
        my $G_AMPLICONDEPTH = 1000000; # one million
        @cols = ('#j5', 'j3', 'num_rp', 'j5AmpId', 'j3AmpId', 'j5AmpDepth', 'j3AmpDepth', 'RPM_mapped', 'RPM_pool', 'RPMPM_mapped_ampliconJ5', 'RPMPM_mapped_ampliconJ3', 'RPMPM_pool_ampliconJ5', 'RPMPM_pool_ampliconJ3');
        #@cols = ('#j5', 'j3', 'num_rp', 'j5AmpId', 'j3AmpId', 'j5AmpDepth', 'j3AmpDepth', 'RPHT_mapped', 'RPHT_pool', 'RPHTPHT_mapped_ampliconJ5', 'RPHTPHT_mapped_ampliconJ3', 'RPHTPHT_pool_ampliconJ5', 'RPHTPHT_pool_ampliconJ3');
        print OUTFILE join("\t", @cols), "\n";
        my $numEvenRP = $readPairs{'even-rp'};
        if (0==$numEvenRP) {
            print "WARNING: total even read-pair has zero entry, resetting to 1 for division\n";
            $numEvenRP = 1;
        }
        my $numOddRP = $readPairs{'odd-rp'};
        if (0==$numOddRP) {
            print "WARNING: total odd read-pair has zero entry, resetting to 1 for division\n";
            $numOddRP = 1;
        }
        foreach my $jumpRef (@jumps) {
            my $numJuncs = scalar(keys %{$jumpRef->{reads}});
            @cols = ($jumpRef->{j5}, $jumpRef->{j3}, $numJuncs, $jumpRef->{j5AmpId}, $jumpRef->{j3AmpId});
            my $j5Id = int($jumpRef->{j5AmpId});
            my $j5Depth = $coveragesRef->[$j5Id]->{'max'};
            my $j3Id = int($jumpRef->{j3AmpId});
            my $j3Depth = $coveragesRef->[$j3Id]->{'max'};
            push @cols, $j5Depth, $j3Depth;
            if (0==$j5Id) {
                print "WARNING: Skipping read which starts with non-amplicon segment\n", Dumper($jumpRef), "\n";
                next;
            }
            if (0==$j3Id) {
                print "WARNING: Skipping read which starts with non-amplicon segment\n", Dumper($jumpRef), "\n";
                next;
            }
            if (0==$j5Depth) {
                print "WARNING: j5 $jumpRef->{j5AmpId} has zero depth, resetting to 1 for division\n";
                $j5Depth = 1;
            }
            if (0==$j3Depth) {
                print "WARNING: j5 $jumpRef->{j5AmpId} has zero depth, resetting to 1 for division\n";
                $j3Depth = 1;
            }

            # junc / mapped: RPHT
            my $normDepth1 = $numJuncs / $readPairs{'total-rp'} * $G_SEQDEPTH;
            push @cols, sprintf("%.2f", $normDepth1);
            # junc / pool: RPHT
            my $normDepth2 = $numJuncs / ((0==($j5Id % 2)) ? $numEvenRP : $numOddRP) * $G_SEQDEPTH;
            push @cols, sprintf("%.2f", $normDepth2);
            # junc5 / amplicon: RPHTPT
            push @cols, sprintf("%.2f", $normDepth1 / $j5Depth * $G_AMPLICONDEPTH);
            # junc3 / amplicon: RPHTPT
            push @cols, sprintf("%.2f", $normDepth1 / $j3Depth * $G_AMPLICONDEPTH);
            # junc5 / amplicon: RPHTPT
            push @cols, sprintf("%.2f", $normDepth2 / $j5Depth * $G_AMPLICONDEPTH);
            # junc3 / amplicon: RPHTPT
            push @cols, sprintf("%.2f", $normDepth2 / $j3Depth * $G_AMPLICONDEPTH);

            print OUTFILE join("\t", @cols), "\n";
        }
        close OUTFILE;
    }

}

sub getSample {
    my ($file) = @_;
    my @bits = split(/\./, $file);
    return ($bits[0], $bits[1]);
}

sub readJumpFile {
    my ($bamFile, $id, $jumpsRef) = @_;

    # 1) COV19-040.GT20-12810.sarcov2.amplicontag.coverage.statistics
    #    total, pool,ults amplicon statistics
    # 2) COV19-040.GT20-12810.sarcov2.amplicontag.coverage.junctions
    #    junctions' result

    my $outPrefix = $bamFile; $outPrefix=~s/\.bam$//;
    my $file = sprintf("%s.coverage.statistics", $outPrefix);
    open INFILE, $file || die "Fail to open $file\n$!\n";
    if (!exists $jumpsRef->{samples}->{keyed}->{$id}) {
        my %sample = (id=>$id, file=>$file, idx=>scalar(@{$jumpsRef->{samples}->{ordered}}));
        $jumpsRef->{samples}->{keyed}->{$id} = \%sample;
        push @{$jumpsRef->{samples}->{ordered}}, \%sample;
    }
    my $sampleRef = $jumpsRef->{samples}->{keyed}->{$id};
    # my $sampleIdx = $jumpsRef->{samples}->{keyed}->{$id}->{idx};
    my $sampleIdx = $sampleRef->{idx};
    while (<INFILE>) {
        chomp();
        my @bits = split(/\t/);
        if (/^#/) {
            last if ('#' ne $bits[0]);
            if (scalar(@bits)>=3) {
                next if ('data_file' eq $bits[1]);
                $sampleRef->{$bits[1]} = $bits[2];
            }
        } else {
            # we skip the amplicon coverage for now
            last;
        }
    }
    close INFILE;


    $file = sprintf("%s.coverage.junctions", $outPrefix);
    open INFILE, $file || die "Fail to open $file\n$!\n";
    my $header = <INFILE>; chomp($header); $header =~ s/^\#//;
    my @headers = split(/\t/, $header);
    my $numHeaders = scalar(@headers);
    while (<INFILE>) {
        chomp();
        my @bits = split(/\t/);
        my %item = ();
        for(my $i=0; $i<$numHeaders; ++$i){
            $item{$headers[$i]} = $bits[$i];
        }
        my $jumpId = sprintf("j5-%d-j3-%d", $item{j5}, $item{j3});
        if (!exists $jumpsRef->{keyed}->{$jumpId}) {
            my %jump = (j5=>$item{j5}, j3=>$item{j3}, counts=>[], j5AmpId=>$item{j5AmpId}, j3AmpId=>$item{j3AmpId});
            $jumpsRef->{keyed}->{$jumpId} = \%jump;
        }
        if (defined $jumpsRef->{keyed}->{$jumpId}->{counts}->[$sampleIdx]) {
            die "Jump $jumpId already present for sample $id\n$_\n",Dumper($jumpsRef->{keyed}->{$jumpId}->{counts}->[$sampleIdx]);
        }
        $jumpsRef->{keyed}->{$jumpId}->{counts}->[$sampleIdx] = \%item;
    }
    close INFILE;
}

sub countSampleTotalJunctionsReads {
    my ($jumpsRef, $minReads) = @_;
    my $orderedSamplesRef = $jumpsRef->{samples}->{ordered};
    foreach my $sampleRef (@{$orderedSamplesRef}) {
        $sampleRef->{totalJuns} = 0;
        $sampleRef->{totalJunsPassed} = 0;
    }
    while (my ($jKey, $jDataRef) = each %{$jumpsRef->{keyed}}) {
        my $countsRef = $jDataRef->{counts};
        for(my $i=0; $i<scalar(@{$countsRef}); ++$i) {
            next if (!defined $countsRef->[$i]);
            my $iSampleRef = $orderedSamplesRef->[$i];
            my $count = $countsRef->[$i]->{num_rp};
            $iSampleRef->{totalJuns} += $count;
            $iSampleRef->{totalJunsPassed} += $count if ($count>=$minReads);
        }
    }
}

sub writeJumpSites {
    my ($file, $jumpSitesRef) = @_;
    open OUTFILE, ">$file" || die "Fail to open $file\n$!\n";
    foreach my $site (sort { $a<=>$b} keys %{$jumpSitesRef}) {
        my $siteRef = $jumpSitesRef->{$site};
        my @cols = ('MN908947.3', $site-1, $site, $site, 0, '+');
        print OUTFILE join("\t", @cols), "\n";
    }
    close OUTFILE;
}

sub annotateJumpsWithJunctions {
    my ($jumpSitesRef) = @_;
    # my $G_ANNOTATION_JUNCTIONS = '../../../sar-cov2-junctions.sorted.bed';
    #    gzip -dc KEL3805A2733.GT20-08106.sarcov2.dedup.jumps.txt.gz \
    #    | awk -v FS="\t" -v OFS="\t" '{print "MN908947.3",$1-1,$1,"j5-"NR"-"$1"-"$2"-c"$3,$3,"+";}' \
    #    | closestBed -D a -a - -b ../../../sar-cov2-junctions.sorted.bed
    my $tmpFile = 'tmp.bed';
    writeJumpSites($tmpFile, $jumpSitesRef);

    my @commands = ('closestBed', '-D', 'a', '-a', $tmpFile, '-b', $G_ANNOTATION_JUNCTIONS);
    push @commands, '|';
    my $command = join(' ', @commands);
    open INFILE, $command || die "Fail to $command\n$!\n";
    while (<INFILE>) {
        chomp();
        my @bits = split(/\t/);
        my $site = $bits[3];
        my $junction = $bits[9];
        my $distance = $bits[12];
        die "Jump site $site is non-existent!\n" if (!exists $jumpSitesRef->{$site});
        $jumpSitesRef->{$site}->{closest} = {} if ($jumpSitesRef->{$site}->{closest});
        $jumpSitesRef->{$site}->{closest}->{$junction} = $distance;
    }
    close INFILE;
}

sub annotateJumpsWithFeatures {
    my ($jumpSitesRef) = @_;
    # my $G_ANNOTATION_FEATURES = '../../../sar-cov2-features';
    #    gzip -dc KEL3805A2733.GT20-08106.sarcov2.dedup.jumps.txt.gz \
    #    | awk -v FS="\t" -v OFS="\t" '{print "MN908947.3",$1-1,$1,"j5-"NR"-"$1"-"$2"-c"$3,$3,"+";}' \
    #    | intersectBed -wo -a - -b ../../../sar-cov2-features.bed
    my $tmpFile = 'tmp.bed';
    writeJumpSites($tmpFile, $jumpSitesRef);

    my @commands = ('intersectBed', '-wo', '-a', $tmpFile, '-b', $G_ANNOTATION_FEATURES);
    push @commands, '|';
    my $command = join(' ', @commands);
    open INFILE, $command || die "Fail to $command\n$!\n";
    while (<INFILE>) {
        chomp();
        my @bits = split(/\t/);
        my $site = $bits[3];
        my $feature = $bits[9];
        die "Jump site $site is non-existent!\n" if (!exists $jumpSitesRef->{$site});
        $jumpSitesRef->{$site}->{overlap} = {} if ($jumpSitesRef->{$site}->{overlap});
        $jumpSitesRef->{$site}->{overlap}->{$feature} = 1;
    }
    close INFILE;
}

sub annotateJumps {
    my ($jumpsRef) = @_;

    # TODO:
    # 1) generate the unique j5 list
    # 2) annotate with junction + features
    #    gzip -dc KEL3805A2733.GT20-08106.sarcov2.dedup.jumps.txt.gz \
    #    | awk -v FS="\t" -v OFS="\t" '{print "MN908947.3",$1-1,$1,"j5-"NR"-"$1"-"$2"-c"$3,$3,"+";}' \
    #    | closestBed -D a -a - -b ../../../sar-cov2-junctions.sorted.bed
    #    gzip -dc KEL3805A2733.GT20-08106.sarcov2.dedup.jumps.txt.gz \
    #    | awk -v FS="\t" -v OFS="\t" '{print "MN908947.3",$1-1,$1,"j5-"NR"-"$1"-"$2"-c"$3,$3,"+";}' \
    #    | intersectBed -wo -a - -b ../../../sar-cov2-features.bed
    # 3) generate the unique j3 list
    # 4) annotate with junction + features
    #    gzip -dc KEL3805A2733.GT20-08106.sarcov2.dedup.jumps.txt.gz \
    #    | cut -f 2 | sort -k1,1n | uniq \
    #    | awk -v FS="\t" -v OFS="\t" '{print "MN908947.3",$1-1,$1,"j3-"$1,0,"+";}' \
    #    | closestBed -D a -a - -b ../../../sar-cov2-junctions.sorted.bed
    #    gzip -dc KEL3805A2733.GT20-08106.sarcov2.dedup.jumps.txt.gz \
    #    | cut -f 2 | sort -k1,1n | uniq \
    #    | awk -v FS="\t" -v OFS="\t" '{print "MN908947.3",$1-1,$1,"j3-"$1,0,"+";}' \
    #    | intersectBed -wo -a - -b ../../../sar-cov2-features.bed

    # 1) generate the unique j5 list
    my %j5s = ();
    while (my ($jKey, $jDataRef) = each %{$jumpsRef->{keyed}}) {
        if (!exists $j5s{$jDataRef->{j5}}) {
            $j5s{$jDataRef->{j5}} = {j5=>$jDataRef->{j5}, refs=>[]};
        }
        push @{$j5s{$jDataRef->{j5}}->{refs}}, $jDataRef;
        $jDataRef->{j5annoRef} = $j5s{$jDataRef->{j5}};
    }
    # 2) annotate with junction + features
    annotateJumpsWithJunctions(\%j5s);
    annotateJumpsWithFeatures(\%j5s);

    # 3) generate the unique j3 list
    my %j3s = ();
    while (my ($jKey, $jDataRef) = each %{$jumpsRef->{keyed}}) {
        if (!exists $j3s{$jDataRef->{j3}}) {
            $j3s{$jDataRef->{j3}} = {j3=>$jDataRef->{j3}, refs=>[]};
        }
        push @{$j3s{$jDataRef->{j3}}->{refs}}, $jDataRef;
        $jDataRef->{j3annoRef} = $j3s{$jDataRef->{j3}};
    }
    # 4) annotate with junction + features
    annotateJumpsWithJunctions(\%j3s);
    annotateJumpsWithFeatures(\%j3s);
}

sub reportJumpsGrid {
    my ($jumpsRef, $minReads, $minSampleSupport, $reportAuxResult, $reportPrefix) = @_;


    my $file = sprintf("%s.grid.absolute.expression.ge%dreadsupport.ge%dsample.xls", $reportPrefix, $minReads, $minSampleSupport);
    open OUTFILE, ">$file" || die "Fail to open $file\n$!\n";

    # report the number of proper pairs?!
    # report a grid table ; append samples at the right columns
    my @cols = ();
    push @cols, ('.') x 11, 'total_read_pair';
    grep { push @cols, $_->{total_read_pair}; } @{$jumpsRef->{samples}->{ordered}};
    print OUTFILE join("\t", @cols), "\n";

    @cols = ();
    push @cols, ('.') x 11, 'odd_read_pair';
    grep { push @cols, $_->{odd_read_pair}; } @{$jumpsRef->{samples}->{ordered}};
    print OUTFILE join("\t", @cols), "\n";

    @cols = ();
    push @cols, ('.') x 11, 'even_read_pair';
    grep { push @cols, $_->{even_read_pair}; } @{$jumpsRef->{samples}->{ordered}};
    print OUTFILE join("\t", @cols), "\n";

    @cols = ();
    push @cols, ('.') x 11, 'total';
    grep { push @cols, $_->{totalJuns}; } @{$jumpsRef->{samples}->{ordered}};
    print OUTFILE join("\t", @cols), "\n";

    @cols = ();
    push @cols, ('.') x 11, 'passed';
    grep { push @cols, $_->{totalJunsPassed}; } @{$jumpsRef->{samples}->{ordered}};
    print OUTFILE join("\t", @cols), "\n";

    # headers
    @cols = ('j5', 'j3', 'j5closest', 'j5closestdist', 'j5overlap', 'j3closest', 'j3closestdist', 'j3overlap', 'j5AmpId', 'j3AmpId', '#sample', '#passedSample');
    grep { push @cols, $_->{id}; } @{$jumpsRef->{samples}->{ordered}};
    grep { push @cols, "rpm_".$_->{id}; } @{$jumpsRef->{samples}->{ordered}};
    if (0!=$reportAuxResult) {
        grep { push @cols, "j5rpmpma_".$_->{id}; } @{$jumpsRef->{samples}->{ordered}};
        grep { push @cols, "j3rpmpma_".$_->{id}; } @{$jumpsRef->{samples}->{ordered}};
    }
    print OUTFILE join("\t", @cols), "\n";

    # data rows! orderd by the junctions
    my @ordered = values %{$jumpsRef->{keyed}};
    @ordered = sort { $a->{j5}<=>$b->{j5} || $a->{j3}<=>$b->{j3} } @ordered;
    my $numColumns = scalar(@{$jumpsRef->{samples}->{ordered}});
    foreach my $juncRef (@ordered) {
        #@cols = ('j5', 'j3', 'j5closest', 'j5overlap', 'j3closest', 'j3overlap', '#sample');
        @cols = ($juncRef->{j5}, $juncRef->{j3});

        # j5 closest, overlap
        my @tags = keys %{$juncRef->{j5annoRef}->{closest}};
        if (scalar(@tags)>0) {
            my @values = (); grep { push @values, $juncRef->{j5annoRef}->{closest}->{$_}; } @tags;
            push @cols, join(";", @tags), join(";", @values);
        } else {
            push @cols, ('.') x 2;
        }
        @tags = keys %{$juncRef->{j5annoRef}->{overlap}};
        push @cols, (scalar(@tags)>0) ? join(";", @tags) : '.';

        # j3 closest, overlap
        @tags = keys %{$juncRef->{j3annoRef}->{closest}};
        if (scalar(@tags)>0) {
            my @values = (); grep { push @values, $juncRef->{j3annoRef}->{closest}->{$_}; } @tags;
            push @cols, join(";", @tags), join(";", @values);
        } else {
            push @cols, ('.') x 2;
        }
        @tags = keys %{$juncRef->{j3annoRef}->{overlap}};
        push @cols, (scalar(@tags)>0) ? join(";", @tags) : '.';

        # let's report the amplicon id for j5 and j3
        push @cols, $juncRef->{j5AmpId}, $juncRef->{j3AmpId};

        # TODO: numSamples .. count the numbers
        # TODO: samples details
        my @values = ();
        my @auxValues = ();
        my @j5Values = ();
        my @j3Values = ();
        my $numSamples = 0;
        my $numPassedSamples = 0;
        my $countsRef = $juncRef->{counts};
        for(my $i=0; $i<$numColumns; ++$i) {
            if (!defined $countsRef->[$i]) {
                push @values, '.';
                push @auxValues, '.';
                push @j5Values, '.';
                push @j3Values, '.';
            } else {
                $numSamples++;
                my $count = $countsRef->[$i]->{'num_rp'};
                if ($count>=$minReads) {
                    die "Zero count col",$i," junc",$juncRef->{j5},"-",$juncRef->{j3},"!\n"  if (0==$count);
                    push @values, $count;
                    push @auxValues, $countsRef->[$i]->{'RPM_mapped'};
                    push @j5Values, $countsRef->[$i]->{'RPMPM_mapped_ampliconJ5'};
                    push @j3Values, $countsRef->[$i]->{'RPMPM_mapped_ampliconJ3'};
                    $numPassedSamples++ ;
                } else {
                    push @values, '.';
                    push @auxValues, '.';
                    push @j5Values, '.';
                    push @j3Values, '.';
                }
            }
        }
        push @cols, $numSamples, $numPassedSamples;
        push @cols, @values;
        push @cols, @auxValues;
        if (0!=$reportAuxResult) {
            push @cols, @j5Values;
            push @cols, @j3Values;
        }

        print OUTFILE join("\t", @cols), "\n";
    }
    close OUTFILE;
}

sub runGridExpression {
    # we now have 
    # 1) COV19-040.GT20-12810.sarcov2.amplicontag.coverage.statistics
    #    total, pool,ults amplicon statistics
    # 2) COV19-040.GT20-12810.sarcov2.amplicontag.coverage.junctions
    #    junctions' result
    #
    # let's generate the junctions (we have read support and number of samples filters)
    #
    # a) junction, j5id, j3id, j5depth, j3depth, summary statistics, samples(read_support)
    # b) absolute expr: junction, j5id, j3id, j5depth, j3depth, summary statistics, samples(RPM_mapped)
    # c) relative expr: junction, j5id, j3id, j5depth, j3depth, summary statistics, samples(RPMPM_mapped_ampliconJ5), samples(RPMPM_mapped_ampliconJ3)
    #
    # read support, ge1 (no filter), ge2 (min. filter)
    #

    my $minReadSupport = shift @ARGV; $minReadSupport = int($minReadSupport);
    my $minSampleSupport = shift @ARGV; $minSampleSupport = int($minSampleSupport);
    my $reportAuxResult = shift @ARGV;
    my $reportPrefix = shift @ARGV;

    # my @jumpsFiles = glob $ARGV[0];

    my %jumps = (keyed=>{}, samples=>{ordered=>[], keyed=>{}});
    for(my $i=0; $i<scalar(@ARGV); ++$i) {
        my $jumpsFile = $ARGV[$i];
        my ($sampleId, $libId) = getSample($jumpsFile);
        my $id = sprintf("%s.%s", $sampleId, $libId);
        readJumpFile($jumpsFile, $id, \%jumps);
    }

    # tally the total counts
    # without and with filter
    countSampleTotalJunctionsReads(\%jumps, $minReadSupport);

    # annotate the jump
    annotateJumps(\%jumps);

    # report the full table
    # TODO: generate report A) - standard
    # TODO: generate report B) - normalized expression
    # TODO: generate report C) - normalized expression normalized amplicon
    reportJumpsGrid(\%jumps, $minReadSupport, $minSampleSupport, $reportAuxResult, $reportPrefix);
}
