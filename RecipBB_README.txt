Author: Marcus Gallagher-Jones
Institution: UCLA Department of chemistry and biochemistry
Email: marcusgj13@gmail.com

			  ###### Introduction ######

Recip_BB is a set of functions for performing reciprocal best blast homology searches
for a library of genes identified by either genomic or proteomic experiments. 
The code is implemented in python3 and requires installation of the biopython suite 
of programmes. 

There are two main scripts for running reciprocal best blast searches.
The first script enables you to perform searches on any machine by sending queries to
the online BLAST servers and parsing those results. This works best with a small list
of genes as when the online servers are repeatedly accessed their can be a large ammount
of overhead. The second script allows blast to be performed locally provided that the
BLAST+ suite of programmes have been installed.

			###### Installation guide ######

To use these programmes it is reccomended that you have python3 installed (although they
are backwards compatible with python2.7)

Install the biopython library

python -m pip install biopython

If you would like to run a local version of the script you will need to install the BLAST+
suite of programs. Instructions for thise can be found here:

https://www.ncbi.nlm.nih.gov/books/NBK52637/ (Windows)
https://www.ncbi.nlm.nih.gov/books/NBK52640/ (Linux)
https://www.ncbi.nlm.nih.gov/books/NBK279671/ (Mac)

In the same directory it is recommended to make a directory containing blast databases
for each genome/proteome you wishto perform BLASTs against. This can be done as follows:

makeblastdb -in some_proteome.fasta -parse_seqids -dbtype prot.


		###### Running the programs #######

There are two different versions of the program, one requires that you have BLAST+ 
installed on your machine (local) the other requires an internet connection so that
the BLAST servers can be accessed. 

Running Recip_BB_homologs_local:

python Recip_BB_homologs_local.py sequences_to_check.fasta ‘output_filename’ ‘query_organism_name’ ‘subject_organism_name’ 
query_database.fasta subject_database.fasta position_in_fasta_file eValue  

e.g.

python Recip_BB_homologs_local.py test.fasta ‘test’ ‘Trypanosoma_brucei_brucei’ ‘Homo_Sapiens’ 
Trypanosoma_proteome.fasta Human_proteome.fasta 1 0.1

Running Recip_BB_homologs.py requires you to give organism names in the form of an Entrez query:

python Recip_BB_homologs.py sequences_to_check.fasta 'output_filename' 'query_organism_name' 'subject_organism_name'
position_in_fasta_file eValue 

e.g.

python Recip_BB_homologs.py test.fasta 'test' 'Trypanosoma brucei brucei [organism]' 'Homo sapiens [organism]'
1 0.1
