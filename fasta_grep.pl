#!/usr/bin/perl
#
#
# Author: Tomasz Gawęda
# Date:   2008-12-19
#
# Description: 
#    The script can grep fasta files using organism description field
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
#

##########################################################################
############ Required third party modules and it's config ################
##########################################################################
use strict;
use warnings;
use Fatal qw{ open close mkdir opendir closedir }; 
use Getopt::Long;
Getopt::Long::Configure(qw{no_ignore_case});

# use Data::Dumper; # debugging

##########################################################################
##########################################################################
######## The code - change only if you know what you're doing ############
##########################################################################
##########################################################################

# command line parsing
my ( $help, $verbose, $debug ) = ( 0, 1, 1 );
usage(0) if $help;

##########################################################################
### Global variables:
##########################################################################

##########################################################################
### Functions:
##########################################################################
## UI
sub debug { print '[d] ' . join( ' ', @_ ) . "\n" if $debug; }
sub info  { print '[i] ' . join( ' ', @_ ) . "\n" if $verbose; }
sub warning { print '[w] ' . join( ' ', @_ ) . "\n"; }

## Display usage information
sub usage {
	my ($exitCode) = @_;
	$exitCode |= 1;
	print "For usage doc please type: perldoc " . $0 . "\n";
	exit($exitCode);
}

##########################################################################
### Main:
##########################################################################
sub main {
	
	
	usage(1) if not scalar(@ARGV) == 2;
	if ( $ARGV[0] =~ m/^-/) {
		usage(1);
	} 
	
	my $pattern = $ARGV[0];
	#$pattern =~ s/[^\w\s=]//ig;
	
	my $inputFileName = $ARGV[1];
	die "File ".$inputFileName." doesn't exists.\n" if not -r $inputFileName;
	
	
	
	my $seq;
	my $title;
	my $i         = 0;

	open( FASTA, '<', $inputFileName );
	print <FASTA> if 0;
	while (<FASTA>) {
		s/^\s+//g;
		s/\s+$//g;                        # clean input from white chars
		if (m/^(>|;)/) {                  # fasta can start by > or ;
			$i++;
			if ( defined $seq and defined $title ) {
				if ( $title =~ m/$pattern/i) {
					print $title . "\n" . $seq . "\n";
				}
				undef $seq;
			}
			$title = $_;
			next;
		}

		# validate FASTA? hell no! :P
		#tr/[a-z]/[A-Z]/;    # upcase all
		#s/[^A-Z]//g;        # ignore all which is not A-Z char
		$seq .= $_;
	}
	close(FASTA);
	
	# code for last item in fasta file
	if ( defined $seq and defined $title ) {
		if ( $title =~ m/$pattern/i) {
			print $title . "\n" . $seq . "\n";
		}
		undef $seq;
	}
	return 0;
}

main();


__END__

=head1 NAME

fasta_grep.pl - grep fasta files using organism description field

=head1 EXAMPLES

fasta_grep.pl 'pattern' file.fast

=head1 REQUIRED ARGUMENTS

=over 4

=item B<pattern>

Pattern should be valid regexp.

=item B<filename>

File name must exist and be readable

=back
