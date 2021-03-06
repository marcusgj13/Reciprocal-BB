from sys import argv
import time
import signal
from Bio import SeqIO
from Bio.Seq import *
from Bio.Alphabet import *
from Bio.SeqUtils import *
from Bio.Blast import NCBIXML, NCBIWWW


# This set of functions are for parsing through a fasta file to extract gene sequences
# performing blast on those sequences are returning the reciprocal best blasts for each
# gene using the Biopython module.
# Program can be run by: python Recip_BB_tryp_homologs.py filename organism_a organism_b
# Marcus Gallagher-Jones    Department of Chemistry and Biochemistry    2017

# These lines are here to handle the fact that BLAST sometimes randomly hangs
# when connection to the server is lost.
class TimeoutException(Exception):   # Custom exception class
    time.sleep(60)
    pass


def timeout_handler(signum, frame):   # Custom signal handler
    raise TimeoutException


# Change the behavior of SIGALRM
signal.signal(signal.SIGALRM, timeout_handler)


def do_blast(seq, organism, eVal):

    while True:
        signal.alarm(120)
        try:
            result = NCBIWWW.qblast("blastp",
                                    "nr",
                                    seq,
                                    entrez_query=organism,
                                    expect=eVal)
            break
        except TimeoutException:
            print("Server timeout, trying again")
            continue

    signal.alarm(0)
    return result


def get_recip_BB(top_three, organism_a, accession_id):
    # This function performs a reciprocal blast on the top three blast hits
    # to determine the most likely orthologous hit.

    test_scores = []
    print("Performing reciprocal BLASTs on top three results\n")
    print("Top result\n")
    top_blast = do_blast(top_three[0].hsps[0].sbjct, organism_a, eVal)
    top_record = NCBIXML.read(top_blast)
    if top_record.alignments[0].accession == accession_id:
        test_scores.append(top_record.alignments[0].hsps[0].score)
    else:
        test_scores.append(float(0))



    print("Second result\n")
    try:
        second_blast = do_blast(top_three[1].hsps[0].sbjct, organism_a, eVal)
        second_record = NCBIXML.read(second_blast)
        if second_record.alignments[0].accession == accession_id:
            test_scores.append(second_record.alignments[0].hsps[0].score)
        else:
            test_scores.append(float(0))
    except IndexError:
        print("No second result\n")
        test_scores.append(float(0))

    print("Third result\n")
    try:
        third_blast = do_blast(top_three[2].hsps[0].sbjct, organism_a, eVal)
        third_record = NCBIXML.read(third_blast)
        if third_record.alignments[0].accession == accession_id:
            test_scores.append(third_record.alignments[0].hsps[0].score)
        else:
            test_scores.append(float(0))
    except IndexError:
        print("No third result\n")
        test_scores.append(float(0))

    print("BLAST complete searching for result with best score\n")


    best_score = max(test_scores)

    if best_score != 0:
        best_index = test_scores.index(best_score)
        recip_BB = top_three[best_index]
    else:
        recip_BB = 'No result found'

    return recip_BB


def write_gene_ids(title2, title1, t, y, output_file, position):
    # This function just writes out the gene ids of the Reciprocal_BB as a text file
    print("=========Starting to write gene ID=========")
    if position == 1:
        with open(output_file, 'w') as f2:
            if type(y) != str:
                print("writing " + y.title + "'s gene ID to file")
                f2.write("> Gene position" + "\t" )
                f2.write(title1 + " geneID \t")
                f2.write(title2 + " geneID \t")
                f2.write("eValue \n")
                title_string = str(y.title)
                title_string2 = str(t.title)
                id = title_string.split('|', 4)
                id2 = title_string2.split('|', 4)
                f2.write(str(position) + "\t")
                f2.write(id[3] + "\t")
                f2.write(id2[3] + "\t")
                f2.write(str(y.hsps[0].expect) + "\n")

            else:
                f2.write("> Gene position" + "\t" )
                f2.write(title1 + " geneID \t")
                f2.write(title2 + " geneID \t")
                f2.write("eValue \n")

    else:
        with open(output_file, 'a') as f2:

            if type(y) != str:
                print("writing " + y.title + "'s gene ID to file")
                title_string = str(y.title)
                title_string2 = str(t.title)
                id = title_string.split('|', 4)
                id2 = title_string2.split('|', 4)
                f2.write(str(position) + "\t")
                f2.write(id[3] + "\t")
                f2.write(id2[3] + "\t")
                f2.write(str(y.hsps[0].expect) + "\n")

    print("=========done writing to file=========")


def write_to_file(y, output_filename, position):
    # This function writes out the results of the Reciprocal_BB as a text file
    print("=========Starting to write output=========")
    if position == 1:
        with open(output_filename, 'w') as f:
            #for y in homologs:
            if type(y) != str:
                print("writing " + y.title + " to file")
                f.write(str(position) + "\n")
                f.write("sequence:" + y.title + "\n")
                f.write("length:" + str(y.length) + "\n")
                f.write("e value:" + str(y.hsps[0].expect) + "\n")
                f.write("score:" + str(y.hsps[0].score) + "\n")
                f.write("Query\n")
                f.write(y.hsps[0].query + "\n")
                f.write(y.hsps[0].match + "\n")
                f.write(y.hsps[0].sbjct + "\n")
                f.write("Subject\n")
                f.write("\n")
                f.write("\n")
            else:
                f.write(str(position) + "\n")
                f.write(y + "\n")
                f.write("\n")
                f.write("\n")
    else:
        with open(output_filename, 'a') as f:
            # for y in homologs:
            if type(y) != str:
                print("writing " + y.title + " to file")
                f.write(str(position) + "\n") 
                f.write("sequence:" + y.title + "\n")
                f.write("length:" + str(y.length) + "\n")
                f.write("e value:" + str(y.hsps[0].expect) + "\n")
                f.write("score:" + str(y.hsps[0].score) + "\n")
                f.write("Query\n")
                f.write(y.hsps[0].query + "\n")
                f.write(y.hsps[0].match + "\n")
                f.write(y.hsps[0].sbjct + "\n")
                f.write("Subject\n")
                f.write("\n")
                f.write("\n")
            else:
                f.write(str(position) + "\n")
                f.write(y)
                f.write("\n")
                f.write("\n")

    print("=========done writing to file=========")


def parse_sequences(filename, gene_source, organism_a, organism_b,
                    start_point, eVal):
    # This function takes a fasta and extracts gene sequences which
    # are subsequently run through BLAST against the genome of a specific organism.
    # The top hits are then used for a reciprocal BLAST and the reciprocal
    # best BLAST hits are output to a text file. Organism_a is the origin of
    # the query sequence and organism_b is the genome to compare to. These
    # should be specified as entrez queries e.g
    # organism_a ='Homo sapiens [organism]'

    # First read in the fasta file
    Sequence_list = [] # These will eventually store the results.
    Homolog_list  = []
    Scores = []
    Gene_ids = []
    E_values = []
    position = 1
    print("You are comparing " 
            + organism_a + " genes against " 
            + organism_b + " genes\n")

    for entry in SeqIO.parse(filename,"fasta",IUPAC.unambiguous_dna):
        if position >= start_point:
            print("Working on entry " + str(entry.description) + 
                    '\n' + "at position " + str(position) + "\n")
            Sequence_list.append(entry)
            # Perform blast and extract the top three hits
            # Have to blast sequence against its own genome to get correct accession ID

            test = do_blast(entry.seq.translate(), organism_a, eVal)

            record = NCBIXML.read(test)
            try:
                top = record.alignments[0]
                accession_id = top.accession

                print("accession ID identified now BLASTing against " + organism_b + " genome\n")

                # Now perform BLAST against search organism

                result = do_blast(entry.seq.translate(), organism_b, eVal)

                blast_record = NCBIXML.read(result)
                top_three = blast_record.alignments[0:3]
                if top_three != []:
                    recip_BB = get_recip_BB(top_three, organism_a, accession_id)
                    Homolog_list.append(recip_BB)
                else:
                    recip_BB = "No result found"
                    Homolog_list.append(recip_BB)
            except IndexError:
                recip_BB = "No result found"
                top = 'No Gene found'

            if type(recip_BB) != str:
                Scores.append(recip_BB.hsps[0].score)
            else:
                Scores.append(0)

            write_to_file(recip_BB, 'Reciprocal_BB_results_' + gene_source +
                          '.txt',position)
            write_gene_ids(organism_a.replace(" ", "_"), organism_b.replace(" ", "_"), top, recip_BB,
                           'Reciprocal_BB_results_' + gene_source + 'geneIDs.txt', position)

            position += 1
        else:
            position += 1


# Get variables from input

script, filename, gene_source, organism_a, organism_b, start_point, eVal = argv
start_point = int(start_point)
eVal = float(eVal)

parse_sequences(filename, gene_source, organism_a,
                organism_b, start_point, eVal)












