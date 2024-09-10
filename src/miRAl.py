#!/usr/bin/env python3

#####################################################################################################################
#                                                                                                                    #
#                                             miRNA couples alignment tools                                            #
#                                                                                                                    #
# Needleman-Wunsch alignment given a list of couples of miRNAs in .tsv format structured as follows:                #
#                 mature miRNA1     |    mature miRNA2                                                                     #
# and a fasta file containing the sequences of the miRNAs perform the alignmetn between the sequences                 #
# of the two miRNAs in the couple                                                                                    #
#                                                                                                                    #
# output is a .tsv file structured as follows:                                                                        #
#                 mature miRNA1     |    mature miRNA2     |    normalized_score                                            #
#                                                                                                                    #
#  -f -a options add a second output file in a simil-fasta format,                                                     #
# where the header contains the names of the two aligned sequences                                                     #
# while the "body" contains the alignment score and the alignment itself                                            #
#                                                                                                                    #
# Author: Leonardo Agasso, 2023                                                                                        #
#                                                                                                                    #
# idea from:                                                                                                        #
# https://bmcgenomics.biomedcentral.com/articles/10.1186/1471-2164-13-218                                            #
#                                                                                                                    #
#####################################################################################################################




#____________________________________________Libraries & Methods_____________________________________________________
import sys
import errno
import warnings
import numpy as np

                                            #BioPython required
from optparse import OptionParser
from Bio import SeqIO
from Bio import pairwise2
from Bio.pairwise2 import format_alignment
#____________________________________________________________________________________________________________________



#______________________________________Global Variables & Parameters_________________________________________________
# Ordered nucleotides to be substituted in seed
nonseed_nt = "ACGU"    # standard nucleotides of a RNA sequence
seed_nt = "ZVRB"        # A→Z, C→V, G→R, U→B : translation for the nucleotides in the seed
all_nt = nonseed_nt+seed_nt

seed_nt_map = {
    "A":"Z",
    "C":"V",
    "G":"R",
    "U":"B"
}

seed_nt_map_reverse = {v: k for k, v in seed_nt_map.items()}
#____________________________________________________________________________________________________________________



#____________________________________________General Functions_______________________________________________________
def ignore_broken_pipe(func):    # Avoid broken pipe error 
    try:
        func()
    except IOError as e:
        if e.errno == errno.EPIPE:
            sys.exit(0)
        else:
            raise


def parse_args():    # Parse command line arguments
    
    parser = OptionParser(usage=format_usage('''
        %prog [OPTIONS] fasta.fa [fasta2.fa] < couples.tsv >output.tsv

        Perform an all-vs-all alignment (using the Needleman-Wunsch algorithm) between the sequences in FASTA,
        returns a qFASTA file where the header contains the names of the two aligned sequences
        while the "body" contains the alignment score and the alignment itself.\033[0m
        Available seed types are (\033[34mblue\033[0m identifies the seed, \033[31mred\033[0m is the last non-seed nucleotide at the 5'):

          •        \033[7m'8merA'\033[0m:                            \033[7m'8merB'\033[0m:
                           ...321                              ...321
                      3'-...N\033[34mNNNNNNNN\033[0m\033[31mN\033[0m-5'                 3'-...N\033[34mNNNNNNNN\033[0mN\033[31mN\033[0m-5'                  
                       ||||||||                            |||||||| 
                ORF...N\033[33mNNNNNNNN\033[0mN...                 ORF...N\033[33mNNNNNNNN\033[0mN...

          •        \033[7m'7merA'\033[0m (default):                  \033[7m'7merB'\033[0m:
                         ...321                               ...321
                      3'-...N\033[34mNNNNNNN\033[0m\033[31mN\033[0m-5'              3'-...N\033[34mNNNNNNN\033[0mN\033[31mN\033[0m-5'
                       |||||||                             |||||||
                ORF...N\033[33mNNNNNNN\033[0mN...                  ORF...N\033[33mNNNNNNN\033[0mN...        
         
          •        \033[7m'6merA'\033[0m:                            \033[7m'6merB'\033[0m:
                        ...321                               ...321 
                      3'-...N\033[34mNNNNNN\033[0m\033[31mN\033[0m-5'                  3'-...N\033[34mNNNNNN\033[0mN\033[31mN\033[0m-5'
                       ||||||                              ||||||
                ORF...N\033[33mNNNNNN\033[0mN...                   ORF...N\033[33mNNNNNN\033[0mN...
        

        In-seed nucleotides are converted in the the following way: A→Z, C→V, G→R, U→B    
    '''    

        
    ))

    parser.add_option('-t', '--seed-type', default='6mer', metavar='STRING',
                        help='Specify the format of the seed you want to consider (default: %default). Options are: 6mer, 7mer, 8mer or null')

    parser.add_option('-g', '--gap-open', default=-5, metavar='INT',
                        help='Specify the gap opening penalty (default: %default)')

    parser.add_option('-e', '--gap-extend', default=-4, metavar='INT',
                        help='Specify the gap extension penalty (default: %default)')

    parser.add_option('-x', '--mismatch-score', default=-1, metavar='INT',
                        help='Specify the mismatch score (default: %default)')
    
    parser.add_option('-m', '--match-score', default=1, metavar='INT',
                        help='Specify the match score (default: %default)')

    parser.add_option('-y', '--seed-mismatch-score', default=-3, metavar='INT',
                        help='Specify the seed mismatch score (default: %default)')

    parser.add_option('-n', '--seed-match-score', default=3, metavar='INT',
                        help='Specify the seed match score (default: %default)')

    parser.add_option('-p', '--penalize-end-gaps', default=False, action='store_true',
                        help='Penalize end gaps (default: %default)')

    parser.add_option('-o', '--one-alignment', default=False, action='store_true',
                        help='Return only one alignment per pair (default: %default)')

    parser.add_option('-s', '--score-only', default=False, action='store_true',
                        help='Return only the alignment score (default: %default)')

    parser.add_option('-f', '--fasta-out', default=False, action="store_true",
                        help="""Change the output format. Produce a simil-fasta format, where the header contains the names
of the two aligned sequences while the "body" contains the alignment
score (absolute and normalized) and the alignment itself.  (default: %default)""")

    parser.add_option('-b', '--not-in-seed-nucleotides', default=1, metavar='INT',
                        help='Specify the number of not-in-seed nucleotides at the 5\' end \
                        (default: %default)')


    return parser.parse_args()


def format_usage(usage):    # Format the usage string
    
    def prefix_length(line):
        length = 0
        while length < len(line) and line[length] in (' ', '\t'):
            length += 1
        return length

    lines = usage.split('\n')
    while len(lines) and len(lines[0].strip()) == 0:
        del lines[0]
    while len(lines) and len(lines[-1].strip()) == 0:
        del lines[-1]

    plen = min(prefix_length(l) for l in lines if len(l.strip()) > 0)
    return '\n'.join(l[plen:] for l in lines)


def read_couples(): # Read the miRNA couples from the file and store them in a list of couples
    couples = []
    
    for line in sys.stdin:
        columns = line.strip().split("\t")
        couples.append((columns[0], columns[1]))
    return couples


def check_sequences(seq1, row):    # Check if the sequences contain characters that are not in the alphabet
    for i in seq1:
        if i not in all_nt:
            raise ValueError("sequence {} (row {}, first sequence) contains the unknown character {}".format(seq1, row, i))

#_____________________________________________Seed Management________________________________________________________
def seed_modify(sequence):        # Convert the seed nucleotides of a miRNA to a new alphabet (AUCG--->ZVRB)
    l1 = ConvertToList(sequence)
    for i in myRange(final_notseed_bps, final_notseed_bps+seed_length-1, 1):
        for j in range(len(nonseed_nt)):
            if(l1[i] == nonseed_nt[j]):
                l1[i] = seed_nt[j]
    
    mod_sequence = ConvertToString(l1)
    return mod_sequence

def seed_modify_all(mirnas):    # Convert the seed nucleotides of all the miRNAs in a dictionary to a new alphabet (AUCG--->ZVRB)
    for key in mirnas.keys():
        mirnas[key] = seed_modify(mirnas[key])
    return mirnas

def seed_restore(sequence):        # Restore the seed nucleotide to the standard alphabet (ZVRB--->AUCG) to display them
    l1 = ConvertToList(sequence)
    
    for i in myRange(len(sequence)-seed_length-final_notseed_bps, len(sequence)-final_notseed_bps-1, 1):
        for j in range(len(nonseed_nt)):
            if(l1[i] == seed_nt[j]):
                l1[i] = nonseed_nt[j]
    
    rest_sequence = ConvertToString(l1)
    return rest_sequence


def ConvertToList(string):        # Convert a string to a list
    list1 = []
    list1[:0] = string
    return list1


def ConvertToString(list):        # Convert a list to a string
    str1 = ""
    for elem in list:
        str1 += elem
    return str1


def myRange(start,end,step):    # Define a range for the for-loop
    i = start
    while i < end:
        yield i
        i += step
        yield end


def seed_definition(options):    # Define the seed length from the command line argument
    global seed_length
    global final_notseed_bps
    
    final_notseed_bps = options.not_in_seed_nucleotides

    # Switch structure to define the seed length
    if options.seed_type=='6mer':
        seed_length = 6
    elif options.seed_type=='7mer':
        seed_length = 7
    elif options.seed_type=='8mer':
        seed_length = 8
    elif options.seed_type=='null':
        seed_length = 0
#____________________________________________________________________________________________________________________



#_____________________________________________Matrices_______________________________________________________________
def subs_matrix(options):                                    # Define the substitution matrix
    
    V = options.match_score                        # Out-of-seed match score
    X = options.mismatch_score                    # Out-of-seed mismatch score
    Y = X                                        # Out-of-seed TRANSITION mismatch score (As for now we are not considering transitions and transversions differently)

    sV = options.seed_match_score                # In-seed match score
    sX = options.seed_mismatch_score            # In-seed mismatch score
    sY = sX                                     # In-seed TRANSITION mismatch score (As for now we are not considering transitions and transversions differently)
    
    # Visualize to better handle the substitution matrix in an intuitive way (since the matrix is symmetric just the superior triangular matrix has to be filled)
    # One improvement could be to differently handle transversion and transition mismatches:
    #    Transition:             A ↔ G, C ↔ T                        (Purine into purine or pyrimidine into pyrimidine), 
    #    Transversion:             A ↔ C, A ↔ T, G ↔ C, G ↔ T            (Purine into pyrimidine or pyrimidine into purine)

    subs_matrix = np.array([
        #    A    |            C    |            G    |            U    |            Z    |            V    |            R    |            B    |
        [    int(V),            int(X),            int(Y),            int(X),            int(sV),        int(X),            int(Y),            int(X)    ],    # A

        [    int(X),            int(V),            int(X),            int(X),            int(X),            int(sV),        int(X),            int(Y)    ],    # C

        [    int(Y),            int(X),            int(V),            int(X),            int(Y),            int(X),            int(sV),        int(X)    ],    # G                
                                                                                                                                                        #-----------------------------IN-SEED PORTION-------------------------------    
        [    int(X),            int(X),            int(X),            int(V),            int(X),            int(X),            int(X),            int(sV)    ],    # U                Z    |                V    |                R    |                B                                                                                                                
        
        [    int(sV),        int(X),            int(Y),            int(X),                                                                                        int(sV),            int(sX),            int(sY),            int(sX)    ],    # Z
                                                                                                                
        [    int(X),            int(sV),        int(X),            int(Y),                                                                                        int(sX),            int(sV),            int(sX),            int(sY)    ],    # V

        [    int(Y),            int(X),            int(sV),        int(X),                                                                                        int(sY),            int(sX),            int(sV),            int(sX)    ],    # R

        [    int(X),            int(Y),            int(X),            int(sV),                                                                                    int(sX),            int(sY),            int(sX),            int(sV)    ]    # B
    ])

    return subs_matrix


def Sim_matr_to_dic(string1, string2, options):                # Convert the substitution matrix into a dictionary (necessary for the alignment function)

    #convert the substitution matrix into a dictionary
    sm = subs_matrix(options)
    list1 = ConvertToList(string1)
    list2 = ConvertToList(string2)
    list3 = list1+list2
    string3 = string1+string2

    dic_subs_matrix = {}
    for i in range(0,len(string3)):
        for j in range(i,len(string3)):
            dic_subs_matrix[list3[i],list3[j]] = sm[i,j]

    return dic_subs_matrix
#____________________________________________________________________________________________________________________



#___________________________________________Output/Sequence Format___________________________________________________
def print_fasta_header(name1, name2, seq1, seq2, file):        # Print the header of the fasta file

    #print every alignment in a quasi FASTA format
    print(">{}: {}\t | \t{}: {}".format(name1, seed_restore(seq1), name2, seed_restore(seq2)), file=file)


def decorate_alignment_output_single_row(row):                # Decorate the single row of the output to highlight the seed
    decorated=""
    for i in row:
        if i in seed_nt_map_reverse.keys():
            decorated += seed_nt_map_reverse[i]
        else:
            decorated += i.lower()
    return decorated
            

def decorate_alignemnt_output(alignment_string):            # Decorate the alignment output to highlight the seed (for every row)
    rows = alignment_string.split("\n")

    seq1 = rows[0]
    seq2 = rows[2]
    seq_symbols = rows[1]

    seq1, seq2, seq_symbols = fix_alignments_match(seq1, seq2, seq_symbols)

    rows[0] = decorate_alignment_output_single_row(rows[0])
    rows[1] = seq_symbols
    rows[2] = decorate_alignment_output_single_row(rows[2])
    return "\n".join(rows)


def fix_alignments_match(seq1, seq2, seq_symbols):

    # convert the sequences into lists
    list1 = ConvertToList(seq1)
    list2 = ConvertToList(seq2)
    list3 = ConvertToList(seq_symbols)

    # for every element of the list
    for i in range(0,len(list3)):
        # if the element is a match
        if list3[i] == ".":
            # if the mismatch involve a seed nucleotide and a non-seed nucleotide
            if list1[i] in seed_nt and list2[i] in nonseed_nt:
                # check whether they are the same nucleotide, just translated
                if seed_nt_map_reverse[list1[i]] == list2[i]:
                    # if they are, change the symbol from "." to "|"
                    list3[i] = "|"
                else:
                    # otherwise, keep the mismatch
                    list3[i] = "."

            if list2[i] in seed_nt and list1[i] in nonseed_nt:
                if seed_nt_map_reverse[list2[i]] == list1[i]:
                    list3[i] = "|"
                else:
                    list3[i] = "."

    # convert the lists back to strings
    seq1 = ConvertToString(list1)
    seq2 = ConvertToString(list2)
    seq_symbols = ConvertToString(list3)

    return seq1, seq2, seq_symbols



#____________________________________________________________________________________________________________________


def read_fasta_file(fasta_file, allow_na=False):    # Read the fasta file and store the sequences in a dictionary
    seq_dict = {}
    with open(fasta_file,"r") as handle:
        record_number = 0
        for record in SeqIO.parse(handle, "fasta"):
            record_number += 1
            # Check if the sequences contain characters that are not in the alphabet
            record.seq=record.seq.replace("T", "U")
            check_sequences(record.seq, record_number)
            if record.id == "N.A." and allow_na:
                continue
            else:
                seq_dict[record.id] = record.seq
    return seq_dict

def main():

    # Ignore the deprecationwarning from Bio.pairwise2
    warnings.filterwarnings("ignore", category=PendingDeprecationWarning)

    # Parse command line arguments
    options, args = parse_args()
    
    # Seed properties
    seed_definition(options)

    # Similarity matrix
    dic_sssmatr = Sim_matr_to_dic(nonseed_nt, seed_nt, options)
    
    # Dictionary to store the names of the mirnas associated with the sequence
    mirseq = {}

    fasta_files_number = len(args)
    if len(args) == 0:
        print("No fasta file provided", file=sys.stderr)
        exit(1)
    if len(args) == 1:
        mirseq = read_fasta_file(args[0], allow_na=True)
        # Couples from stdin 
        couples = read_couples()
    if len(args) == 2:
        mirseq = seed_modify_all(read_fasta_file(args[0]))
        mirseq2 = seed_modify_all(read_fasta_file(args[1]))
        #populate couples with all the possible couples between the two fasta files
        couples = []
        for i in mirseq.keys():
            for j in mirseq2.keys():
                couples.append((i, j))
    else:
        print("Too many arguments", file=sys.stderr)
        exit(1)

    # For every couple of miRNAs
    row = 0
    for i,j in couples:
        if fasta_files_number == 1:
            seq_i = seed_modify(mirseq[i])
            seq_j = seed_modify(mirseq[j])
        else: 
            seq_i = seed_modify(mirseq[i])
            seq_j = seed_modify(mirseq2[j])
            #already modified

#
#    All alignment functions have the following arguments:
#        - Two sequences:                     strings, Biopython sequence objects or lists. Lists are useful for supplying sequences which contain residues that are encoded by more than one letter.
#        - penalize_extend_when_opening:     boolean (default: False). Whether to count an extension penalty when opening a gap. If false, a gap of 1 is only penalized an "open" penalty, otherwise it is penalized "open+extend".
#        - penalize_end_gaps:                 boolean. Whether to count the gaps at the ends of an alignment. By default, they are counted for global alignments but not for local ones. Setting penalize_end_gaps to (boolean, boolean) 
#                                             allows you to specify for the two sequences separately whether gaps at the end of the alignment should be counted.
#        - gap_char:                         string (default: '-'). Which character to use as a gap character in the alignment returned. If your input sequences are lists, you must change this to ['-'].
#        - force_generic:                     boolean (default: False). Always use the generic, non-cached, dynamic programming function (slow!). For debugging.
#        - score_only:                         boolean (default: False). Only get the best score, don't recover any alignments. The return value of the function is the score. Faster and uses less memory.
#        - one_alignment_only:                 boolean (default: False). Only recover one alignment.
#    SOURCE:    https://biopython.org/docs/1.75/api/Bio.pairwise2.html    


    ####################################################################################################
        a = pairwise2.align.globalds(seq_i,                                                            #                
                                    seq_j,                                                                #
                                    dic_sssmatr,                                                        #
                                    options.gap_open,                                                    #
                                    options.gap_extend,                                                #
                                                                                                        #
                                    one_alignment_only = options.one_alignment,                        #
                                    score_only = options.score_only,                                    #
                                    penalize_end_gaps = options.penalize_end_gaps,                        #
                                    penalize_extend_when_opening = False)                               #
    ####################################################################################################

        #normalized_score = mathematically_accurate_normalized_score(a[0][2], options, seq_i, seq_j, dic_sssmatr)
        
        

        #Standard output with the first mirna, the second mirna and the score as a .tsv file
        formatted_alignment = decorate_alignemnt_output(format_alignment(*a[0]))
        if options.fasta_out:
            print("{}\t{}\t{}".format(i, j, a[0][2]))
        else:
            #generate another output file in quasi-FASTA format
            print(">{}\t{}\t{}".format(i, j, a[0][2]))
            print(formatted_alignment, file=sys.stdout)
    
    row += 1



if __name__=='__main__':

    # Ignore the deprecationwarning from Bio.pairwise2
    warnings.filterwarnings("ignore", category=PendingDeprecationWarning) 
    
    ignore_broken_pipe(main)



    