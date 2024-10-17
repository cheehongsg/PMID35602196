#!env perl

##
## Covid19 - Analysis
## Copyright (C) 2020-2021 Chee-Hong WONG (Dr Chia-Lin Wei Laboratory)
##

$|++;
use strict;
use Data::Dumper;
use FindBin;

# take a series of jumps and translate the product
# take the longest version? or keep the different possiblities


# https://en.wikipedia.org/wiki/DNA_codon_table
my %G_AA = (
    'F' => {aa=>'F', aa3=>'Phe', property=>'Nonpolar'}
    , 'L' => {aa=>'L', aa3=>'Leu', property=>'Nonpolar'}
    , 'I' => {aa=>'I', aa3=>'Ile', property=>'Nonpolar'}
    , 'M' => {aa=>'M', aa3=>'Met', property=>'Nonpolar'}
    , 'V' => {aa=>'V', aa3=>'Val', property=>'Nonpolar'}
    , 'S' => {aa=>'S', aa3=>'Ser', property=>'Polar'}
    , 'P' => {aa=>'P', aa3=>'Pro', property=>'Nonpolar'}
    , 'T' => {aa=>'T', aa3=>'Thr', property=>'Polar'}
    , 'A' => {aa=>'A', aa3=>'Ala', property=>'Nonpolar'}
    , 'Y' => {aa=>'Y', aa3=>'Tyr', property=>'Polar'}
    , '*' => {aa=>'*', aa3=>'***', property=>'Stop codon'}
    , 'H' => {aa=>'H', aa3=>'His', property=>'Basic'}
    , 'Q' => {aa=>'Q', aa3=>'Gln', property=>'Polar'}
    , 'N' => {aa=>'N', aa3=>'Asn', property=>'Polar'}
    , 'K' => {aa=>'K', aa3=>'Lys', property=>'Basic'}
    , 'D' => {aa=>'D', aa3=>'Asp', property=>'Acidic'}
    , 'E' => {aa=>'E', aa3=>'Glu', property=>'Acidic'}
    , 'C' => {aa=>'C', aa3=>'Cys', property=>'Polar'}
    , 'W' => {aa=>'W', aa3=>'Trp', property=>'Nonpolar'}
    , 'R' => {aa=>'R', aa3=>'Arg', property=>'Basic'}
    , 'G' => {aa=>'G', aa3=>'Gly', property=>'Nonpolar'}
);
my %G_STD_CODON = (
    #
    'TTT' => 'F', 'TTC' => 'F', 'TTA' => 'L', 'TTG' => 'L'
    , 'CTT' => 'L', 'CTC' => 'L', 'CTA' => 'L', 'CTG' => 'L'
    , 'ATT' => 'I', 'ATC' => 'I', 'ATA' => 'I', 'ATG' => 'M'
    , 'GTT' => 'V', 'GTC' => 'V', 'GTA' => 'V', 'GTG' => 'V'
    #
    , 'TCT' => 'S', 'TCC' => 'S', 'TCA' => 'S', 'TCG' => 'S'
    , 'CCT' => 'P', 'CCC' => 'P', 'CCA' => 'P', 'CCG' => 'P'
    , 'ACT' => 'T', 'ACC' => 'T', 'ACA' => 'T', 'ACG' => 'T'
    , 'GCT' => 'A', 'GCC' => 'A', 'GCA' => 'A', 'GCG' => 'A'
    #
    , 'TAT' => 'Y', 'TAC' => 'Y', 'TAA' => '*', 'TAG' => '*'
    , 'CAT' => 'H', 'CAC' => 'H', 'CAA' => 'Q', 'CAG' => 'Q'
    , 'AAT' => 'N', 'AAC' => 'N', 'AAA' => 'K', 'AAG' => 'K'
    , 'GAT' => 'D', 'GAC' => 'D', 'GAA' => 'E', 'GAG' => 'E'
    #
    , 'TGT' => 'C', 'TGC' => 'C', 'TGA' => '*', 'TGG' => 'W'
    , 'CGT' => 'R', 'CGC' => 'R', 'CGA' => 'R', 'CGG' => 'R'
    , 'AGT' => 'S', 'AGC' => 'S', 'AGA' => 'R', 'AGG' => 'R'
    , 'GGT' => 'G', 'GGC' => 'G', 'GGA' => 'G', 'GGG' => 'G'
);
my %G_5PrimeLeader = (start=>1, end=>85, name=>'TRS-L');


my $genomeFile = $FindBin::Bin.'/../annotation/'.'MN908947.3.fasta';
my $annotationFile = $FindBin::Bin.'/../annotation/'.'SARS-CoV-2-annotations.bed';

my $genome = loadGenome($genomeFile);
my %annotations = ();
loadAnnotations($annotationFile, \%annotations);

# calculate the standard ORFs
loadStandardORFs($genome, \%annotations);


my $jumpsFile = $ARGV[0];
reportTranslation($jumpsFile, $genome, \%annotations);

exit 0;

sub loadStandardORFs {
    my ($genome, $annotationsRef) = @_;
    # let's load the sequences of CDS
    foreach my $annotationRef (@{$annotationsRef->{byCoord}}) {
        my $isNSP = 0;
        if ($annotationRef->{name}=~/^nsp\d+/) {
            $annotationRef->{aaLen} = -1;
            $annotationRef->{seq} = '';
            $isNSP = 1;
        }

        my $finalSequence = substr($genome, $annotationRef->{start}-1, $annotationRef->{end}-$annotationRef->{start}+1);
        # attempt to translate!
        my $seqLen = length($finalSequence);
        my @frames = ();
        my @marks = ([],[],[]);
        my $endFrame = (0==$isNSP) ? 2 : 0;
        for(my $frame=0; $frame<=$endFrame; ++$frame) {
            my @aas = ();
            for(my $i=$frame; $i<$seqLen; $i+=3) {
                my $codon = substr($finalSequence, $i, 3);
                last if (length($codon)<3);
                my $aa = $G_STD_CODON{$codon};
                push @aas, $aa;
                if ('M' eq $aa) {
                    push @{$marks[$frame]}, {frame=>$frame, ntpos=>$i, aapos=>$#aas, aa=>$aa};
                } elsif ('*' eq $aa) {
                    push @{$marks[$frame]}, {frame=>$frame, ntpos=>$i, aapos=>$#aas, aa=>$aa};
                }
            }
            $frames[$frame] = \@aas;
        }

        if (0==$isNSP) {
            # extract the longest?
            # extract the first 2 MET?
            my %longest = (aaLen=>0, frame=>-1, ntstart=>-1, ntend=>-1, aastart=>-1, aaend=>-1, seq=>'');
            for(my $frame=0; $frame<=2; ++$frame) {
                my $marksRef = $marks[$frame];
                my $metMark = -1;
                for(my $i=0; $i<scalar(@{$marksRef}); ++$i) {
                    if ('M' eq $marksRef->[$i]->{aa}) {
                        $metMark = $i;
                        last;
                    }
                }
                if (-1!=$metMark) {
                    my $stopMark = -1;
                    for(my $i=$metMark+1; $i<scalar(@{$marksRef}); ++$i) {
                        if ('*' eq $marksRef->[$i]->{aa}) {
                            $stopMark = $i;
                            last;
                        }
                    }

                    if (-1!=$stopMark) {
                        if (-1==$longest{ntstart} || $longest{ntstart}>$marksRef->[$metMark]->{ntpos}) {
                            my $aaSpan = $marksRef->[$stopMark]->{aapos} - $marksRef->[$metMark]->{aapos} + 1;
                            $longest{aaLen} = $aaSpan;
                            $longest{frame} = $frame;
                            $longest{ntstart} = $marksRef->[$metMark]->{ntpos};
                            $longest{ntend} = $marksRef->[$stopMark]->{ntpos};
                            $longest{aastart} = $marksRef->[$metMark]->{aapos};
                            $longest{aaend} = $marksRef->[$stopMark]->{aapos};
                            $longest{seq} = '';
                            for(my $i=$longest{aastart}; $i<=$longest{aaend}; ++$i) {
                                $longest{seq} .= $frames[$frame]->[$i];
                            }
                        }
                    }
                }
            }
            $annotationRef->{aaLen} = $longest{aaLen};
            $annotationRef->{seq} = $longest{seq};
        } else {
            # nsp, just take the first frame!
            $annotationRef->{aaLen} = scalar(@{$frames[0]});
            $annotationRef->{seq} = join('', @{$frames[0]});
        }
    }
}

# TODO: options
# 1) we might have to get the first or second MET
# 2) we might have to get the longest possible
sub getCdsOverlap {
    my ($site, $annotationsRef) = @_;

    my $siteCdsRef = undef;
    if ($site<=85) {
        # we assume that this is from TRS-L, some suboptimal alignment
        $siteCdsRef = \%G_5PrimeLeader;
    } else {
        # find the ORF which the site is in and assume it is affected
        my $jump5 = $site;
        my $cdssRef = $annotationsRef->{byCoord};
        for(my $i=0; $i<scalar(@{$cdssRef}); ++$i) {
            my $cdsRef = $cdssRef->[$i];
            if (($cdsRef->{start}-1)<=$jump5 && $jump5<=$cdsRef->{end}) {
                # this is the ORF that we are interested
                $siteCdsRef = $cdsRef;
                last;
            }
        }
    }
    return $siteCdsRef;
}

sub getProducts {
    my ($genome, $sitesRef, $annotationsRef, $productsRef) = @_;
    @{$productsRef} = ();

    return if (0==scalar(@{$sitesRef}));

    # TODO: record the j5 and j3 hits so that we can predict the kind of impact
    # TODO: if the data is from pacbio, we will not need to guess, thus optional
    my $assumeStart = 0; # start is 0-based
    my $j5CdsRef = getCdsOverlap($sitesRef->[0], $annotationsRef);
    my $j3CdsRef = getCdsOverlap($sitesRef->[1], $annotationsRef);
    if (defined $j5CdsRef) {
        $assumeStart = $j5CdsRef->{start}-1;
    } else {
        if (0==$assumeStart) {
            my $jump5 = $sitesRef->[0];
            print STDERR "\nWARNING: Fail to locate a good start (overlap with cds) for $jump5..$sitesRef->[1] (line $.)!\n";
            $assumeStart = $jump5 - 75; # we have read lenght of 150, and we assume half is covered!
        }
    }


    my @positions = ($assumeStart);    # start is 0-based
    push @positions, @{$sitesRef};
    push @positions, length($genome);
    
    my @exons = ();
    for(my $i=0; $i<scalar(@positions); $i+=2) {
        push @exons, {start=>$positions[$i]+(0==$i?0:-1), end=>$positions[$i+1]};
    }

    my $finalSequence = '';
    for(my $i=0; $i<scalar(@exons); ++$i) {
        my $exonRef = $exons[$i];
        my $exonSeq = substr($genome, $exonRef->{start}, $exonRef->{end}-$exonRef->{start});
        # print STDERR $exonRef->{start}, "\t", $exonRef->{end}, "\t", $exonSeq, "\n";
        $finalSequence .= $exonSeq;
    }

    # attempt to translate!
    my $seqLen = length($finalSequence);
    my @frames = ();
    my @marks = ([],[],[]);
    for(my $frame=0; $frame<=2; ++$frame) {
        my @aas = ();
        for(my $i=$frame; $i<$seqLen; $i+=3) {
            my $codon = substr($finalSequence, $i, 3);
            last if (length($codon)<3);
            my $aa = $G_STD_CODON{$codon};
            push @aas, $aa;
            if ('M' eq $aa) {
                push @{$marks[$frame]}, {frame=>$frame, ntpos=>$i, aapos=>$#aas, aa=>$aa};
            } elsif ('*' eq $aa) {
                push @{$marks[$frame]}, {frame=>$frame, ntpos=>$i, aapos=>$#aas, aa=>$aa};
            }
        }
        $frames[$frame] = \@aas;
    }

    # extract the longest?
    # extract the first 2 MET?
    my %longest = (aaLen=>0, frame=>-1, ntstart=>-1, ntend=>-1, aastart=>-1, aaend=>-1, seq=>'');
    $longest{j5CdsRef}=$j5CdsRef;
    $longest{j3CdsRef}=$j3CdsRef;
    for(my $frame=0; $frame<=2; ++$frame) {
        my $marksRef = $marks[$frame];
        my $metMark = -1;
        for(my $i=0; $i<scalar(@{$marksRef}); ++$i) {
            if ('M' eq $marksRef->[$i]->{aa}) {
                $metMark = $i;
                last;
            }
        }
        if (-1!=$metMark) {
            my $stopMark = -1;
            for(my $i=$metMark+1; $i<scalar(@{$marksRef}); ++$i) {
                if ('*' eq $marksRef->[$i]->{aa}) {
                    $stopMark = $i;
                    last;
                }
            }

            if (-1!=$stopMark) {
                if (-1==$longest{ntstart} || $longest{ntstart}>$marksRef->[$metMark]->{ntpos}) {
                    my $aaSpan = $marksRef->[$stopMark]->{aapos} - $marksRef->[$metMark]->{aapos} + 1;
                    $longest{aaLen} = $aaSpan;
                    $longest{frame} = $frame;
                    $longest{ntstart} = $marksRef->[$metMark]->{ntpos};
                    $longest{ntend} = $marksRef->[$stopMark]->{ntpos};
                    $longest{aastart} = $marksRef->[$metMark]->{aapos};
                    $longest{aaend} = $marksRef->[$stopMark]->{aapos};
                    $longest{seq} = '';
                    for(my $i=$longest{aastart}; $i<=$longest{aaend}; ++$i) {
                        $longest{seq} .= $frames[$frame]->[$i];
                    }
                }
            }
        }
    }
    push @{$productsRef}, \%longest;
}

sub reportTranslation {
    my ($file, $genome, $annotationsRef) = @_;

    my $fh = undef;
    if (!defined $file || '-' eq $file) {
        $fh = *STDIN;
    } else {
        open $fh, $file || die "Fail to open $file\n$!\n";
    }

    my @cols = ('#lineId');
    grep { push @cols, $_; } ('aaLen', 'frame', 'ntstart', 'ntend', 'aastart', 'aaend', 'j5p', 'j3p', 'match', 'seq');
    print join("\t", @cols), "\n";

    my $block = 10;
    my $trigger = $block;
    my $line = 0;
    while (<$fh>) {
        chomp();
        my @bits = split(/[\t,]/);
        my @products = ();
        getProducts($genome, \@bits, $annotationsRef, \@products);

        @cols = ($.);
        grep { push @cols, $products[0]->{$_}; } ('aaLen', 'frame', 'ntstart', 'ntend', 'aastart', 'aaend');
        push @cols, (defined $products[0]->{j5CdsRef}) ? $products[0]->{j5CdsRef}->{name} : '.';
        push @cols, (defined $products[0]->{j3CdsRef}) ? $products[0]->{j3CdsRef}->{name} : '.';
        # TODO: brute force right now
        my $productLabel = '.';
        foreach my $annotationRef (@{$annotationsRef->{byCoord}}) {
            if ($annotationRef->{seq} eq $products[0]->{'seq'}) {
                $productLabel = $annotationRef->{name};
                last;
            }
        }
        push @cols, $productLabel;
        push @cols, $products[0]->{'seq'};
        print join("\t", @cols), "\n";

        $line++;
        if ($line>=$trigger) {
            print STDERR $line,'.. ';
            $trigger += $block;
        }
    }
    if (!defined $file || '-' eq $file) {
    } else {
        close INFILE;
    }
    print STDERR $line,'.. Done.', "\n";
}

sub loadAnnotations {
    my ($file, $annotationsRef) = @_;
    open INFILE, $file || die "Fail to open $file\n$!\n";
    # exclude TRS-*, CDS-frameshift
    my @annotations = ();
    my %annotations = ();
    while (<INFILE>) {
        chomp();
        my @bits = split(/\t/);
        next if ($bits[3]=~/^TRS-/);
        next if ($bits[3]=~/frameshift/);
        
        my %item = (start=>$bits[1]+1, end=>$bits[2], name=>$bits[3]);
        $item{name} =~ s/^CDS-//;
        push @annotations, \%item;
        $annotations{$item{name}} = \%item;
    }
    close INFILE;
    @annotations = sort { $a->{start}<=>$b->{start} || $a->{end}<=>$b->{end}} @annotations;
    $annotationsRef->{byCoord} = \@annotations;
    $annotationsRef->{byName} = \%annotations;
}

sub loadGenome {
    my ($file) = @_;
    open INFILE, $file || die "Fail to open $file\n$!\n";
    my $seq = '';
    while (<INFILE>) {
        chomp();
        if (/^>/) {
            # igore id
        } else {
            $seq .= $_;
        }
    }
    close INFILE;

    return $seq;
}
