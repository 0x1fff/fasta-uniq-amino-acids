#!/usr/bin/perl

# Author: Tomasz Gawęda
# Date:   2008-04-25
#
# Description: 
#    This script can convert PDB file to FASTA using 
#      - SEQ Marker data 
#      - using only data which is in PDB file
# 
# 
# License: 
#    Apache License Version 2.0, January 2004 (https://tldrlegal.com/ ; http://choosealicense.com/)


use strict;
use warnings;
use Getopt::Long;
use Fatal qw(open close);

Getopt::Long::Configure ( qw{no_ignore_case} );

my $progname = $0;
$progname =~ s/(\S+)\..*/$1/g;

# command line parsing
my ( @u_input_files, $u_out_file_name, $debug, $verbose, $help, $u_use_seq_marker, $u_no_seq_marker );
$u_use_seq_marker = 1;    # default try to use seq marker
GetOptions(
	"i|in=s{1,}"       => \@u_input_files,
	"o|out=s"          => \$u_out_file_name,
	"s|use_seq_marker" => \$u_use_seq_marker,
	"n|no_seq_marker"  => \$u_no_seq_marker,
	"v|verbose"        => \$verbose,
	"h|?|help"         => \$help,
) or usage();
usage() if ( not scalar @u_input_files ) or defined($help);
$u_use_seq_marker = 0 if ($u_no_seq_marker);

##########################################################################
### Functions:
###
##########################################################################

sub usage {
	print STDERR $progname
	  . ' -i file.pdb [file2.pdb file3.pdb ...] -o [output file]' . "\n\n";
	exit(1);
}

sub getResidue {
	my ($pdbline) = @_ or &err("No data in pdbSplit __CODE_ERROR__");

## Using substr instead of regexp - easier to maintain and probably faster

#	ATOM    297  CG2 ILE A  50       4.799  57.014  22.113  1.00 18.67           C
#	HETATM31834  C1  OXL A 533      10.168  68.452  11.859  1.00 53.37           C
#
# http://www.wwpdb.org/documentation/format23/sect9.html
#
	my @ret;

	#push( @ret, substr($pdbline,  6, 6) );     # SERIAL NUM- columns 7-11
	push( @ret, substr( $pdbline, 17, 4 ) );    # RES NAME  - column 18 - 20
	push( @ret, substr( $pdbline, 22, 4 ) );    # RES SEQ   - column 23 - 26

	# strip whitespace
	map { s/^\s+//; s/\s+$//; $_ = uc $_; } @ret;

	#print "m: ".join(' ', @ret)."\n";
	return @ret;
}

#from http://www.rcsb.org/pdb/docs/format/pdbguide2.2/part_79.html
my %aamap = (
	'ALA' => 'A', 	'ARG' => 'R',
	'ASN' => 'N',	'ASP' => 'D',
	'ASX' => 'B',	'CYS' => 'C',
	'GLN' => 'Q',	'GLU' => 'E',
	'GLX' => 'Z',	'GLY' => 'G',
	'HIS' => 'H',	'ILE' => 'I',
	'LEU' => 'L',	'LYS' => 'K',
	'MET' => 'M',	'PHE' => 'F',
	'PRO' => 'P',	'SER' => 'S',
	'THR' => 'T',	'TRP' => 'W',
	'TYR' => 'Y',	'VAL' => 'V',
	'PCA' => 'X'
);
$aamap{'UNK'} = ' ';

sub main {
	open( OUT, '>', $u_out_file_name ) if defined $u_out_file_name;

	for my $file (@u_input_files) {
		my $lastres   = -1;
		my $titleSeq  = q{};
		my $markerSeq = q{};
		my $atomSeq   = q{};

		if ( not -s $file ) {
			warn( $file . ' is not an file is empty or does not exist' . "\n" );
			next;
		}
		#print 'Parsing ' . $file . "\n";

		if ( not defined $u_out_file_name ) {
			my $outFileName = $file;
			$outFileName =~ s/\.\w{0,3}$//g;
			$outFileName .= '.fasta';

			## cleanup
			open( OUT, '>', $outFileName );
		}

		open( IN, '<', $file );
		print <IN> if 0;    # dummy entry
		while (<IN>) {
			s/^s+//g;
			s/\s+$//g;
			if ( substr( $_, 0, 5 ) eq 'TITLE' ) {    #/^TITLE.*/) {
				$_ = substr( $_, 10 );
				$titleSeq .= $_;
				next;
			}

			if ( $u_use_seq_marker and substr( $_, 0, 6 ) eq 'SEQRES' ) {
				$_ = substr( $_, 19 );
				my $curr;
				for $curr ( split( /\s+/, $_ ) ) {
					$curr = uc $curr;
					#warn("Invalid residue $curr in $file")
					next  if not exists $aamap{$curr};
					$markerSeq .= $aamap{$curr};
				}
				next;
			}    #SEQRES

			if ( substr( $_, 0, 4 ) eq 'ATOM' ) {
				my @res = getResidue($_);
				if ( $lastres != $res[1] ) {
					my $curr = uc $res[0];
					#warn("Invalid residue $curr in $file")
					next  if not exists $aamap{$curr};
					$lastres = $res[1];
					$atomSeq .= $aamap{$curr};
				}
				next;
			}    # ATOM
		}    # <IN>
		close(IN);

		print OUT ">" . $titleSeq . "\n";
		if ( $u_use_seq_marker and ( length $markerSeq ) ) {
			print OUT $markerSeq . "\n";
		}
		else {
			if ( ( not length($markerSeq) ) and $u_use_seq_marker ) {
				warn( 'No SEQ markers using SEQ from atoms' . "\n" );
			}
			print OUT $atomSeq . "\n";
		}
		close(OUT) if ( not defined $u_out_file_name );
	}
}

main() unless caller;


__END__

=head1 NAME

pdb2fasta.pl - convert files in PDB format to FASTA format.

=head1 EXAMPLES

=over 4

=item B<pdb2fasta.pl -s -i 1AQ2.pdb>

Convert one PDB file to fasta

=item B<pdb2fasta.pl -s -i 1AQ2.pdb 1AYL.pdb -o 01.fasta>

Convert many PDB files to one FASTA file (with many entries)

=back

=head1 OPTIONS

=over 4

=item B<-i, --in>

Input files (at least one)

=item B<-o, --out>

Output file name

=item B<-s, --use_seq_marker>

Try to convert PDB using data from SEQ section

=item B<-n, --no_seq_marker>

Ignore SEQ seqction in PDB file - just do conversion

=item B<-v, --verbose>

Be verbose

=item B<-h, --?, --help>

Print help and exit

=back
