#!/usr/bin/perl

# Author: Tomasz Gawęda
# Date:   2009-01-17
#
# Description: 
#    The script sorts polypeptides according to the number 
#    of polypeptide occurrences for the defined amino acid sequence or totally at random. 
#    The counting is performed by single amino acids or by a defined polypeptide length.
# 
# Paper: 
#    Statistical distribution of amino acid sequences: a proof of Darwinian evolution
#    Krystian Eitner, Uwe Koch, Tomasz Gawęda and Jędrzej Marciniak
#    Bioinformatics (2010) 26 (23): 2933-2935.
#    DOI: 10.1093/bioinformatics/btq571
#    PMID: 21030460 (http://www.ncbi.nlm.nih.gov/pubmed/21030460) ;
# 
# License: 
#    Apache License Version 2.0, January 2004 (https://tldrlegal.com/ ; http://choosealicense.com/)

use warnings;
use strict;
use File::Basename;
use File::Spec;
use Getopt::Long;
use Fatal qw{ open close mkdir opendir closedir }; # do not install here readdir
Getopt::Long::Configure(qw{no_ignore_case});

##########################################################################
# You may want to change this sort string
my $u_sortString = 'FWLIYCVMDAPTGKHESNQR';
##########################################################################

##########################################################################
##########################################################################
############ Change only if you know what you're doing ###################
##########################################################################
##########################################################################

my $progname = File::Basename::basename($0);
my $version  = q { $Revision: 0.1 $ };
$version =~ s/^[^0-9]+([0-9\.]+).*$/$1/g;

# command line parsing
my ( $help, $verbose );

# prefix u_ - user input
my ( $u_iFName,      $u_oFName, $u_seqLen,     $u_search, $u_sort );
my ( $u_onlyCorrect, $u_norm,   $u_genMissing, $u_perc,   $u_abs );

GetOptions(
	'i|in=s'             => \$u_iFName,
	'o|out=s'            => \$u_oFName,
	'l|sequenceLenght=i' => \$u_seqLen,

	'c|onlyCorrectAA' => \$u_onlyCorrect,
	'a|absolute'      => \$u_abs,
	'p|percentage'    => \$u_perc,
	'n|normalize'     => \$u_norm,

	'g|generateMissing' => \$u_genMissing,
	't|typeSearch=s'    => \$u_search,
	's|sortType=s'      => \$u_sort,

	'v|verbose' => \$verbose,    # produce some information durning run
	'h|?|help'  => \$help,
  )
  or usage(1);
usage(0) if $help;

##########################################################################
### Global variables:
###
##########################################################################

my %uniqAA;                      # count of uniq aa
my %correctAA = (
	'A' => 'ALA',	'R' => 'ARG',
	'N' => 'ASN',	'D' => 'ASP',
	'C' => 'CYS',	'Q' => 'GLN',
	'E' => 'GLU',	'G' => 'GLY',
	'H' => 'HIS',	'I' => 'ILE',
	'L' => 'LEU',	'K' => 'LYS',
	'M' => 'MET',	'F' => 'PHE',
	'P' => 'PRO',	'S' => 'SER',
	'T' => 'THR',	'W' => 'TRP',
	'Y' => 'TYR',	'V' => 'VAL',
	#'B' => 'ASX',	#'Z' => 'GLX',
);

## this is used when changing type of search
my %search = (
	'offset' => \&offsetUniqCount,
	'full'   => \&fullUniqCount,
);

## this is used when selecting type of sort
# sas - sort and save
my %sortAndSaveNoGen = (
	'key-asc' => \&sasKeyAsc,
	'key-usr' => \&sasKeyUsr,
	'val-cnt' => \&sasValCnt,
);

## this is used when selecting type of sort and generating missing AA
# sasm - sort and save with missing
my %sortAndSaveGen = (
	'key-asc' => \&sasmKeyAsc,
	'key-usr' => \&sasmKeyUsr,
	'val-cnt' => \&sasmValCnt,
);

my $sortAndSave;    # holds alternative from %sortAndSaveGen/sortAndSaveNoGen
my $g_max    = -1;
my $g_min    = 99999;
my $g_allCnt = 0;
my $g_zeroCntStr ; # this is used when printing not existing aa into file

##########################################################################
### Functions:
###
##########################################################################

## Display usage information
sub usage {
	my ($ret) = @_;
	print "For usage doc please type: perldoc " . $progname . "\n";
	exit($ret);
}
## Display verbose informations
sub info { print join( ' ', shift ) . "\n" if $verbose; }

sub validateUserParams {

	# check program internals
	if ( length($u_sortString) != scalar( keys %correctAA ) ) {
		die "Sort string is incorrect:\t"
		  . $u_sortString . "\n"
		  . "Sorted string:\t"
		  . join( "", sort split( //, $u_sortString ) )
		  . "All valid AA:\t"
		  . join( "", sort keys %correctAA ) . "\n"
		  . "You're probably missing some AA\n\n";
	}

	die "No input file defined (-i)\n"       if not defined $u_iFName;
	die "Input file doesn't exists (-i)\n"   if not -r $u_iFName;
	die "Not defined sequence lenght (-l)\n" if not defined $u_seqLen;
	die "Seq lenght is not positive number (-l)\n"
	  if not( $u_seqLen =~ m/^\d+$/g );

	die "Normalization requres: ignoring unknown AA (-c)\n"
	  if defined $u_norm
	  and not defined $u_onlyCorrect;

	# validate output format
	$u_norm = 0 if not defined $u_norm;
	$u_perc = 0 if not defined $u_perc;
	$u_abs = 0 if not defined $u_abs;
	$u_abs = 1 if not $u_norm and not $u_perc and not $u_abs;
	
	$u_genMissing = 0 if not defined $u_genMissing;
	$u_onlyCorrect = 0 if not defined $u_onlyCorrect;

	## validate options fallback to defaults
	# validate search type
	$u_search = ( defined $u_search ) ? lc $u_search : 'full';
	if ( not exists $search{$u_search} ) {
		die "Wrong search type - "
		  . $u_search . "\n"
		  . "Possible options for switch -t are: "
		  . join( ", ", keys %search ) . "\n";
	}

	# validate sort type
	$sortAndSave =
	  ( $u_genMissing ) ? \%sortAndSaveGen : \%sortAndSaveNoGen;
	$u_sort = ( defined $u_sort ) ? lc $u_sort : 'val-cnt';
	if ( not exists $sortAndSave->{$u_sort} ) {
		die "Wrong sort type - " . $u_sort . "\n"
		  . "Possible options to switch -s are: \n"
		  . join( ", ", keys %{$sortAndSave} ) . "\n";
	}

	if ( $u_sort eq 'key-str' and not $u_onlyCorrect ) {
		die "Sort type " . $u_sort . " requires option -c\n";
	}

	if ( $u_genMissing and not $u_onlyCorrect ) {
		die "Generating missing AA requires option -c\n";
	}

	if ( not defined $u_oFName ) {
		my $iFBName = File::Basename::basename($u_iFName);
		$iFBName =~ s/\.[a-zA-Z]+$//g;

		my $aInfo = $u_search . '_';
		$aInfo .= ( $u_onlyCorrect ) ? 'only-correct_' : 'full-data_';
		$aInfo .= ( $u_abs )         ? 'absolute_'     : '_';
		$aInfo .= ( $u_norm )        ? 'normalized_'   : '_';
		$aInfo .= ( $u_perc )        ? 'percentage_'   : '_';
		$aInfo .= ( $u_genMissing ) ? 'generated-missing_' : '_';
		$aInfo .= 'sort_' . $u_sort;
		$u_oFName = $iFBName . '_' . $u_seqLen . '_' . $aInfo . '.stat';
	}
	return;
}

sub deleteInvalidAA {
	my ( $validAA, $allAA ) = @_;

	my $regexp = '[^' . join( '|', keys %{$validAA} ) . ']+';
	for my $aa ( keys %{$allAA} ) {
		if ( $aa =~ $regexp ) {
			print STDERR "Invalid AA seq $aa - deleting\n";
			delete $allAA->{$aa};
		}
	}
}

##########################################################################
##### SearchFunctions
##########################################################################
sub fullUniqCount {
	my ( $title, $s, $uSeqLen ) = @_;
	$uSeqLen--;    # [0..2] = 3 elems but we need 2
	               # iterate through whole sequence ($# stands for max index)
	my @seq = split( //, $$s );
	for ( my ( $a, $z ) = ( 0, $uSeqLen ) ; $z < $#seq ; $a++, $z++ ) {
		my $useq = join( '', @seq[ $a .. $z ] );
		$uniqAA{$useq} = 0 if not exists $uniqAA{$useq};
		$uniqAA{$useq}++;
	}
	return 0;
}

sub offsetUniqCount {
	my ( $title, $s, $uSeqLen ) = @_;
	$uSeqLen--;
	my @seq = split( //, $$s );
	for ( my $i = $uSeqLen ; $i < $#seq ; $i += $uSeqLen ) {
		my $useq = join( '', @seq[ $i - $uSeqLen .. $i ] );
		$uniqAA{$useq} = 0 if not exists $uniqAA{$useq};
		$uniqAA{$useq}++;
	}
	return 0;
}
##########################################################################
##### /SearchFunctions
##########################################################################

##########################################################################
##### SortFunctions - No Generating missing
##########################################################################
sub sasValCnt {
	my ( $fileHandle, $allAA ) = @_;
	foreach my $key ( sort { $uniqAA{$b} <=> $uniqAA{$a} } keys %uniqAA ) {
		printAA( $fileHandle, $key, $allAA );
	}
	return;
}

sub sortAndSaveByRandom {
	my ( $fileHandle, $allAA ) = @_;
	while ( my $key = keys %uniqAA ) {
		printAA( $fileHandle, $key, $allAA );
	}
}

sub sasKeyAsc {
	my ( $fileHandle, $allAA ) = @_;
	foreach my $key ( sort { $a cmp $b } keys %uniqAA ) {
		printAA( $fileHandle, $key, $allAA );
	}
	return;
}

sub sasKeyUsr {

	# NOTE: tr/// must have strings in compile time :( at least in perl 5.8
	my ( $fileHandle, $allAA ) = @_;

	my %fwtrTab;
	my %revtrTab;
	my $i = ord('A');

	for my $c ( split( //, $u_sortString ) ) {
		$fwtrTab{$c} = chr($i);
		$revtrTab{ chr($i) } = $c;
		$i++;
	}

	# transliterate keys
	my @trkeys;
	for my $aa ( keys %uniqAA ) {
		my $newAA;
		for my $char ( split( //, $aa ) ) {
			$newAA .= $fwtrTab{$char};
		}
		push @trkeys, $newAA;

		#$printMe{$newAA} = $uniqAA{$aa};
	}

	# sort, transliterate back and print value
	foreach my $key ( sort { $a cmp $b } @trkeys ) {
		my $oldAA;
		for my $char ( split( //, $key ) ) {
			$oldAA .= $revtrTab{$char};
		}
		printAA( $fileHandle, $oldAA, $allAA );
	}
}
##########################################################################
##### /SortFunctions - No Generating missing
##########################################################################

##########################################################################
##### SortFunctions - Generating missing (recursive)
##########################################################################
sub sasmValCnt {
	my ( $fileHandle, $allAA ) = @_;
	sasValCnt( $fileHandle, $allAA );
	my @alph = sort keys %correctAA;
	sasmValCntGen( $fileHandle, $allAA, \@alph, $u_seqLen, my @tmp );
	return;
}

sub sasmValCntGen {
	my ( $fileHandle, $allAA, $alphabet, $genLen, $currGen ) = @_;

	# if we have desired lenght
	if ( $genLen <= 0 ) {
		my $seq = join( '', @{$currGen} );
		printAA( $fileHandle, $seq, $allAA) if not exists $allAA->{$seq};
		return;
	}

	# generating - recursion starts
	for ( my $i = 0 ; $i < scalar(@$alphabet) ; $i++ ) {
		$currGen->[ $u_seqLen - $genLen ] = $alphabet->[$i];
		sasmValCntGen( $fileHandle, $allAA, $alphabet, $genLen - 1, $currGen );
	}
}

sub sasmKeyAsc {
	my ( $fileHandle, $allAA ) = @_;
	my @alph = sort keys %correctAA;
	sasmKeyGen( $fileHandle, $allAA, \@alph, $u_seqLen, my @tmp );
}

sub sasmKeyUsr {
	my ( $fileHandle, $allAA ) = @_;

	# generate reverse tables
	my @alph = split( //, $u_sortString );
	sasmKeyGen( $fileHandle, $allAA, \@alph, $u_seqLen, my @tmp );
}

sub sasmKeyGen {
	my ( $fileHandle, $allAA, $alphabet, $genLen, $currGen ) = @_;

	# if we have desired lenght
	if ( $genLen <= 0 ) {
		printAA( $fileHandle, join( '', @{$currGen} ), $allAA );
		return;
	}

	# generating - recursion starts
	for ( my $i = 0 ; $i < scalar(@$alphabet) ; $i++ ) {
		$currGen->[ $u_seqLen - $genLen ] = $alphabet->[$i];
		sasmKeyGen( $fileHandle, $allAA, $alphabet, $genLen - 1, $currGen );
	}
}
##########################################################################
##### /SortFunctions - Generating missing (recursive)
##########################################################################

sub printAA {
	my ( $fileHandle, $name, $allAA ) = @_;
	my $outStr = '';

	my $count = ( exists $allAA->{$name} ) ? $allAA->{$name} : 0;

	# print absolute values
	if ( exists $allAA->{$name} ) {
		if ($u_abs) { # absolute output
			$outStr .= $count;
		}

		# print normalized values
		if ($u_norm) { # normalized output
			my $norm = $count;
			$norm -= $g_min;
			$norm /= $g_max - $g_min;
			$outStr .= "\t" . $norm;
			#$count = sprintf('%.5d', $count);
		}

		# print percentage values
		if ($u_perc) { # percentage output
			my $perc = (100 * $count)/$g_allCnt;
			$outStr .= "\t".$perc;
		}
	}
	else {    # if not exists $allAA->{$name}
		$outStr =  $g_zeroCntStr;
	}

	print $fileHandle $name . "\t" . $outStr . "\n";
	return;
}

sub getNormalizationFactors {
	my ($allAA) = @_;

	for my $tmp ( values %{$allAA} ) {
		$g_min = $tmp if $tmp < $g_min;
		$g_max = $tmp if $tmp > $g_max;
	}

	# search for missing AA if found set 0 as minimum
	if ( $u_genMissing ) {
		my @tmp = keys %correctAA;
		$g_min = 0 if anyMissingAASeq( \@tmp, $u_seqLen, my @s );
		$g_max = 0 if 0 > $g_max;
	}
	return;
}

sub anyMissingAASeq {
	my ( $alphabet, $genLen, $currGen ) = @_;
	if ( $genLen <= 0 ) {
		my $seq = join( '', @{$currGen} );
		return 1 if ( not exists $uniqAA{$seq} );
		return 0;
	}
	for ( my $i = 0 ; $i < scalar(@$alphabet) ; $i++ ) {
		$currGen->[ $u_seqLen - $genLen ] = $alphabet->[$i];
		return 1 if anyMissingAASeq( $alphabet, $genLen - 1, $currGen );
	}
}

##########################################################################
##### MAIN:
##########################################################################
sub main {

	my $currentSeq;
	my $title;
	my @z = ();

	#genSeq( [ 'a', 'b', 'c' ], $u_seqLen, \@z );

	validateUserParams();

	# minus one len seq 3 is from 0-3 (4 elements)
	#$u_seqLen--;

	# Go throught fasta file
	open( FASTA, '<', $u_iFName );    # not able to open file module Fatal!!
	print <FASTA> if 0;
	while ( my $entry = <FASTA> ) {
		map { s/^\s+//g; s/\s+$//g; } $entry;    # clean input from white chars

		if ( $entry =~ m/^(>|;)/ ) {             # fasta header starts by > or ;
			if ( defined $currentSeq && defined $title ) {

				#print $title." ".$seq."\n";
				# for each entry
				$search{$u_search}->( $title, \$currentSeq, $u_seqLen );
				undef $currentSeq;
			}
			$title = $entry;
			next;
		}

		# validate FASTA
		$entry =~ tr/[a-z]/[A-Z]/;    # upcase all
		$currentSeq .= $entry;
	}
	close(FASTA);
	$search{$u_search}->( $title, \$currentSeq, $u_seqLen );

	# cleanup - invalid aa
	deleteInvalidAA( \%correctAA, \%uniqAA ) if ( $u_onlyCorrect );

	# start normalization
	getNormalizationFactors( \%uniqAA ) if ( $u_norm );
	
	# count all found AA
	$g_allCnt = 0;
	map { $g_allCnt += $_ } values %uniqAA if $u_perc;
	

	$g_zeroCntStr = 0        if $u_abs;
	$g_zeroCntStr .= "\t" . 0 if $u_norm;
	$g_zeroCntStr .= "\t" . 0 if $u_perc;

	# SAVE RESULTS
	open( my $fileHandle, '>', $u_oFName );

	# print file header:
	print $fileHandle "AminoAcid group\t";
	print $fileHandle "Absolute count\t"      if $u_abs;
	print $fileHandle "After normalisation\t" if $u_norm;
	print $fileHandle "Percentage count\t"    if $u_perc;
	print $fileHandle "\n";
	

	$sortAndSave->{$u_sort}->( $fileHandle, \%uniqAA );

	close($fileHandle);

	return 0;
}

main();

__END__

=head1 NAME

fasteUniqAASeq.pl - Searches and list all uniq sequence in protein file


=head1 EXAMPLES

The simplest run (full search, result sorted by abs values)

fasteUniqAASeq.pl -i B<1A49.fasta> -l 4

Exact search with normalization and generating missing AA:

fastaUniqAASeq.pl -i 001.fasta -l 3 -t full -s key-usr -n -g -c

Offset search with normalization

fastaUniqAASeq.pl -i 001.fasta -l 3 -t offset -s key-usr -n -c

=head1 REQUIRED ARGUMENTS

=over 4

=item inputFile

=item uniqSequenceLenght

=back

=head1 OPTIONS

=over 4

=item B<-i, --in>

Input file name and path.

=item B<-l, --sequenceLength>

Lenght of uniq sequence - if higher value you use here the more memmory you have to got.
B<WARNING:> The best way to search for uniq FASTA (without ending memmory is):

MAX_INPUT_FASTA_FILE_SIZE = (AVAIL_PHYS_MEMMORY)/(LENGTH_OF_UNIQ_PROT_SEQ - 1.50)

=item B<-o, --out>

Output file name.

=item B<-f, --offsetSearch>

Enables offset mode search. Offset mode is turned off in default.

=item B<-s, --sortType>

How to sort results?

=over 6

=item key-asc

Sort alphabetically by key (protein short name).

=item key-str

Sort by sort string in program (user specific order).

=item val-count

Sort by popularity of sequence. This is default.

=back

=item B<-c, --onlyCorrectAA>

Without unknown amino acids acids with 'X' default this is not set.

=item B<-a, --absolute>

Prints absolute numbers of found AA - can be used with B<-n> and B<-p>.

=item B<-n, --normalize>

Perform normalization on all values - can be used with B<-a> and B<-p>..

=item B<-p, --percentage>

Prints percentage count of each AA - can be used with B<-a> and B<-n>.

=item B<-g, --generateMissing>

Generate missing AA, it's good for statistical check.

=item B<-v, --verbose>

Be verbose

=item B<-h, --help>

Display short usage help

=back

=head1 DESCRIPTION

This software generates statistics of triplets quartets and quintets across a class of proteines passed as input.

=head1 DEPENDENCIES

=over 4

=item B<perl 5.6> or newer with standard modules

=back
