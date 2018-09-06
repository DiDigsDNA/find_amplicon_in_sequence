#!/usr/bin/env python
'''This script is the first in a set of custom scripts to facilitate in-silico screening of 
sequences for exact PCR primers and probe hybridization sites to screen for strains of target
organism that may not be detected by lab-developed PCR assays. 
This script searches a nucleotide fasta file for an exact match to the PCR amplicon, using information
contained in an assay-specific oligonucleotide (csv) definition file. The amplicon is created
from primers and probe in the correct order, with intervening regions of 0-100 N's.
Each query fasta sequence is read as a Biopython SeqRecord and searched for an exact match to the amplicon
(or its reverse complement) - allowing degenerate bases - and sequences lacking an exact match are output
to fasta for downstream analysis.
A summary containing record IDs of non-matching sequences is printed to a text file.'''

'''cmd line example: python find_amplicon_in_sequence.py ncbi_sequences.fasta oligos.csv my_output'''
'''Author: Diane Eisler, Molecular Microbiology & Genomics, BCCDC Public Health Laboratory,Feb 2018'''

import sys,string,os, time, Bio, re
from Bio import Seq, SeqIO, SeqUtils, Alphabet, SeqRecord
from Bio.Alphabet import IUPAC
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import Bio.Data.IUPACData

fasta_to_parse = sys.argv[1] #fasta file of nucleotide sequences to search in
oligos_to_find = sys.argv[2] # csv primer definition file of oligo sequences to search for
gb_accessions_handle = sys.argv[3] + ".txt" #user-specified output file name
gb_accessions= open(gb_accessions_handle,'w') # text output file of seqences without exact match

def output_results_file(seqRecord_list):
    '''Prints results report, with sequence identifiers of non-matching sequences, for the end-user.'''
    localtime = time.asctime(time.localtime(time.time())) #date and time of analysis
    gb_accessions.write("---------------------------------------------------------------------------\n")
    gb_accessions.write("RESULTS OF SEARCH FOR EXACT MATCHES TO OLIGONUCLEOTIDES IN QUERY SEQUENCES:\n")
    gb_accessions.write("---------------------------------------------------------------------------\n\n")
    gb_accessions.write("Sequences Searched: %s \n\n" % (fasta_to_parse))
    gb_accessions.write("Oligonucleotides Searched for: %s \n\n" % (oligos_to_find))
    gb_accessions.write("Analysis as of: " + localtime + "\n\n")
    gb_accessions.write("The following Genbank accession numbers require further analysis in Geneious:\n")
    gb_accessions.write("\n\n---------------------------------------------------------------------------\n")
    if len(seqRecord_list) > 0: #if there are sequences to investigate, print id to list
        print("Exact match to amplicon NOT found in %i sequences:" % (len(seqRecord_list))) #message to console
        output_fasta_handle = sys.argv[3] + ".fasta" #from user=specified output name
        output_fasta = open(output_fasta_handle, 'w') #output fasta containing sequences requiring further analysis
        for record in seqRecord_list:
            print(record.id) #print record id to console
            print_recordID_to_text_file(record) #write id's to txt file
        SeqIO.write(seqs_mismatched_to_amplicon, output_fasta_handle, "fasta")
    else: #if no sequences require further analysis, print a message to outfile
        gb_accessions.write("No sequences requiring further investigation.")
        print("No sequences requiring further investigation.")
    print("\nEND OF RESULTS\n")
    gb_accessions.write("\n---------------------------------------------------------------------------")
    gb_accessions.write("\n\nEND OF RESULTS")
    gb_accessions.write("\nOutput file: " + gb_accessions_handle)
    print("Search Results Output file: " + gb_accessions_handle) #print search results txt file
    if len(seqRecord_list) > 0: #direct users to unmatching sequences fasta file and print to console
        gb_accessions.write("\nFasta sequences for analysis: " + output_fasta_handle + "\n")
        print ("Fasta sequences for further analysis: " + output_fasta_handle + "\n")
    else: #no fasta file produced
        gb_accessions.write("\nFasta sequences for analysis: N/A\n")
        print("No fasta sequences require further analysis.")
    return

def print_recordID_to_text_file(record):
    '''Print the SeqRecord id to text information file.'''
    gb_accessions.write('\n' + record.id)
    return

def create_amplicon(primer_list):
    '''Creates a theoretical "amplicon" from a list of primer/probe SeqRecords, with a variable
    length spacer of any N's between oligos and returns an amplicon Sequence object.'''
    amplicon_search_string = ''
    spacer = 'N{0,100}' #intervening sequence between oligos (1-100 N's)
    for i in range(len(primer_list)):
        oligo = primer_list[i]
        if i > 0: #precede oligo sequence with spacer for all but first oligo
            amplicon_search_string += spacer + str(oligo.seq)
        else: #just add oligo sequence
            amplicon_search_string += str(oligo.seq)
    #print("amplicon_search_string = ", amplicon_search_string)
    amplicon = Seq(amplicon_search_string, IUPAC.ambiguous_dna)
    return amplicon

def nt_search_regex(seq, subseq): 
    """Search for a DNA subseq in sequence (forward strand only), allowing using degenerate nucleotides
    (e.g. N = A or T or C or G, R = A or G etc.)""" 
    pattern = '' 
    for nt in subseq:
        if nt in Bio.Data.IUPACData.ambiguous_dna_values.keys():
            value = Bio.Data.IUPACData.ambiguous_dna_values[nt]
        else:
            value = nt
        if len(value) == 1: 
            pattern += value 
        else: 
            pattern += '[%s]' % value 

    pos = -1 
    result = [pattern] 
    l = len(seq) 
    while True: 
        pos += 1 
        s = seq[pos:]
        m = re.search(pattern, s) 
        if not m: 
            break 
        pos += int(m.start(0)) 
        result.append(pos) 
    return result

def parse_seqRecord_for_amplicon(record, amplicon):
    '''Parse record(SeqRecord) for amplicon (or reverse complement), returning True if found and False if not.'''
    results = nt_search_regex(str(record.seq),str(amplicon)) #search in SeqRecord sequence for amplicon
    if (len(results) > 1): #if the results contain match and position, the amplicon has been found
        print("Amplicon FOUND in %s at %i position(s):" % (record.id, (len(results)-1)))
        return True
    else:
        return False

#read in oligo properties from csv and create SeqRecords to represent each oligo
with open(oligos_to_find, 'r') as oligo_file:
    oligo_file.readline() #skip the header line
    lines = oligo_file.readlines()
    f_primer = '' #placeholders for oligos
    f_primer = ''
    probe = ''
    for line in lines: #split remaining lines at commas
        oligo_values= line.rstrip().split(',') #assign oligo properties
        [oligo_name,sequence,oligo_type,orientation] = oligo_values[0:4]
        new_seq = Seq(sequence, IUPAC.ambiguous_dna) #create new Sequence object from oligo
        if (oligo_type.lower() == 'primer'): #check if oligo is primer and which strand
            if (orientation.upper() == 'F'): #create forward primer SeqRecord
                f_primer = SeqRecord(new_seq, id = oligo_name)
            elif (orientation.upper() == 'R'): #rev-comp Sequence -->create R primer SeqRecord
                r_primer = SeqRecord(new_seq.reverse_complement(), id = oligo_name)
            else: #print mssg if primer orientation not defined in oligo file
                print("ERROR: primer orientation should be 'F' or 'R'!")
        elif (oligo_type.lower() == 'probe'): #check if oligo is probe and on which strand
            if (orientation.upper() == 'F'): #create probe SeqRecord
                probe = SeqRecord(new_seq, id = oligo_name)
            elif (orientation.upper() == 'R'): #rev-comp Sequence --> create probe SeqRecord
                probe = SeqRecord(new_seq.reverse_complement(), id = oligo_name)
            else: #print mssg if probe orientation not defined in oligo file
                print("ERROR: probe orientation should be 'F' or 'R'!")
        else: #print mssg to console if oligo not defined as primer or probe in file
            print("ERROR in oligonucleotide definition file!")
    #create list of primers along with a reversed list and create amplicon and revcomp to search
    primer_list = [f_primer, probe, r_primer]
    amplicon = create_amplicon(primer_list)
    #reverse complement each member of primer_list and reverse the list using list comprehension
    rc_primer_list = [p.reverse_complement() for p in primer_list[::-1]]
    revcomp_amplicon = create_amplicon(rc_primer_list)
    print("Amplicon        : " + str(amplicon))
    print("Amplicon RevComp: " + str(revcomp_amplicon))

with open(fasta_to_parse,'r') as inFile:
    #read nucleotide sequences from fasta file into SeqRecords, uppercase and add to a list
    seq_list = [rec.upper() for rec in list(SeqIO.parse(inFile, "fasta", alphabet=IUPAC.ambiguous_dna))]
    seqs_mismatched_to_amplicon = [] #empty list to hold SeqRecords without exact match to amplicon
    print("\nSearching for oligonucleotides in %i sequences..." % len(seq_list))
    for record in seq_list: #parse SeqRecord for exact match to amplicon or reverse complement
        found_on_forward_strand = parse_seqRecord_for_amplicon(record, amplicon)
        found_on_reverse_strand = parse_seqRecord_for_amplicon(record, revcomp_amplicon)
        #if exact match to amplicon or its reverse complement not found
        if (found_on_forward_strand == False) and (found_on_reverse_strand == False):
            seqs_mismatched_to_amplicon.append(record) #add SeqRecord to list for further investigation
    #write detailed results report for follow-up by end-user      
    output_results_file(seqs_mismatched_to_amplicon)  
    #output list of sequences needing further investigation to fasta file
    #SeqIO.write(seqs_mismatched_to_amplicon, output_fasta_handle, "fasta")           

inFile.close()
gb_accessions.close()
