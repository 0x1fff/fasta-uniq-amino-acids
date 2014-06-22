fasta-uniq-amino-acids
======================

![DOI FROM zenodo.org](https://zenodo.org/badge/4195/0x1fff/fasta-uniq-amino-acids.png)


The repository contains set of scripts which were used to analyze the quantity of amino acids, dipeptides and tripeptides for all proteins available in the UNIPROT–TREMBL database and the listing for selected species and enzymes (Fasta Unique Sequences Amino Acids Search). UNIPROT–TREMBL contains protein sequences associated with computationally generated annotations and large-scale functional characterization.

Result of this research was published in Bioinformatics.

Paper describing usage of applications
---------------------------------------------

Version 1.0 of this applications were used durning research in paper

> Statistical distribution of amino acid sequences: a proof of Darwinian evolution
> Krystian Eitner, Uwe Koch, Tomasz Gawęda and Jędrzej Marciniak
>
> Bioinformatics (2010) 26 (23): 2933-2935.
> DOI: 10.1093/bioinformatics/btq571
> PMID: 21030460 (http://www.ncbi.nlm.nih.gov/pubmed/21030460) ;

Online version of paper: http://bioinformatics.oxfordjournals.org/content/26/23/2933.short (Full version available for free)

In the article we have proved that the distribution of amino acids, dipeptides and tripeptides is statistical which confirms that the evolutionary biodiversity development model is subject to the theory of independent events. It seems interesting that certain short peptide combinations occur relatively rarely or even not at all. First, it confirms the Darwinian theory of evolution and second, it opens up opportunities for designing pharmaceuticals among rarely represented short peptide combinations. Calculation results prove that the number of di- and tripeptide sequences is consistent with their theoretical number calculated based on the amino acid content in the TREMBL database. This results from the fact that amino acids have similar or even identical biosynthesis pathways in different species. Therefore, the number of amino acids and their combinations (polypeptides) is generally similar. Our calculations provide a new strategy for inhibitor design based on short natural amino acid sequences.


Fasta Unique Sequences Amino Acids Search Script description
-------------------------------------------------------------------

A script fasta_uniq_seq.pl counts the number of polypeptides with a defined length. Due to software limitations, only the listing of all the available di- (400) and tripeptides (8000) could be currently visualized. The script sorts polypeptides according to the number of polypeptide occurrences for the defined amino acid sequence or totally at random. The counting is performed by single amino acids or by a defined polypeptide length.

FASTA is a linear format for the description of proteins and DNA/RNA. Files in the FASTA format contain information about the types of amino acids in a protein. The format is very useful for searching for similar proteins in databases (BLAST, CLUSTALW projects). The largest FASTA databases are AS follows: UNIPROT, Sprot and Protein Data Bank (PDB). Owing to the software developed (fasta_uniq_seq.pl), we can search for recurring sequences of any length in FASTA files in two modes:

<dl>
  <dt>Normal search (precise search, searching sequences one by one amino acid).</dt>
  <dd>The precise search involves the listing of all amino acid sequences with a defined length (shift by one amino acid). </dd>

  <dt>Offset search (search with jumps by a defined sequence length).</dt>
  <dd>The offset search, in turn, involves shifting by a defined length (offset) with fragments at a distance of the offset being included. This is best explained by an example in Table 1 in publication</dd>
</dl>

Furthermore, the script can filter results to avoid situations when unknown amino acids (X) occur in the listing. One of the options is to generate sequences which do not occur in a group of sequences. This is useful for statistical analysis. 

The program can sort results of searching for unique amino acid sequences in various ways:

 * key-asc—alphabetic sorting, ascending by the amino acid name.

 * key-random—no sorting.

 * key-usr—sorting by a user-defined sequence (changed in the programm source code).

 * val-count—sorting by the frequency of an amino acid, descending.

The program can print results using three methods. The most important is to count the physical number of sequences. In addition, it provides normalized quantity and percentage of a sequence. 


When analyzing the whole UNIPROT–TREMBL file, there may be need to split  FASTA files while retaining their structure, according to the number of entries in the resulting FASTA file or size of the resulting file (enable the analysis of large protein sets using computers without sufficient memory).


Fasta Grep Amino Acids Script description
-------------------------------------------------------------------

fasta_grep.pl is a script used to retrieve respective protein or enzyme types or any combinations of characters in the description of a sequence (following the ‘>’ character) from large FASTA files. 


Bonus script: Convert PDB File to FASTA format
------------------------------------------------------

Script pdb2fasta.pl can easily convert files in PDB format to FASTA format. 
It can work in two modes:

<dl>
  <dt>Using SEQRES Marker Section from PDB file</dt>
  <dd>Each PDB formatted file includes "SEQRES records" which list the primary sequence of the polymeric molecules present in the entry. This sequence information is also available as a FASTA download. This listing includes the sequence of each chain of linear, covalently-linked standard or modified amino acids or nucleotides. It may also include other residues that are linked to the standard backbone in the polymer. Chemical components or groups covalently linked to side-chains (in peptides) or sugars and/or bases (in nucleic acid polymers).</dd>

  <dt>Ignoring SEQRES Marker Section from PDB file</dt>
  <dd>This mode ignores SEQRES Marker and tries to convert PDB file using only informations which are defined in other sections of PDB file. This mode is fallback when there is no SEQRES section available.</dd>
</dl>


Requirements for all scripts
----------------------------------

* Perl 5.6 or newer
* perl-doc


Additional documentation:
-------------------------------------------------------------------

You can read additional documentation using perl-doc package.

	$ perldoc fasta_uniq_seq.pl
    $ perldoc fasta_grep.pl
    $ perldoc pdb2fasta.pl


License:
-------------

Apache License Version 2.0, January 2004 (https://tldrlegal.com/ ; http://choosealicense.com/)

