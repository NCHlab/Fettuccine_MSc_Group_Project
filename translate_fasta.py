
# coding: utf-8

# In[31]:


##from https://github.com/mtarbit/Rosalind-Problems/blob/master/e017-orf.py
##Translates sequences from fasta file to retrieve all possible peptides.
##Specification of minimum peptide length possible.
##User specifies input and output file name in command line using sys.argv.

DNA_CODON_TABLE = {
    'TTT': 'F',     'CTT': 'L',     'ATT': 'I',     'GTT': 'V',
    'TTC': 'F',     'CTC': 'L',     'ATC': 'I',     'GTC': 'V',
    'TTA': 'L',     'CTA': 'L',     'ATA': 'I',     'GTA': 'V',
    'TTG': 'L',     'CTG': 'L',     'ATG': 'M',     'GTG': 'V',
    'TCT': 'S',     'CCT': 'P',     'ACT': 'T',     'GCT': 'A',
    'TCC': 'S',     'CCC': 'P',     'ACC': 'T',     'GCC': 'A',
    'TCA': 'S',     'CCA': 'P',     'ACA': 'T',     'GCA': 'A',
    'TCG': 'S',     'CCG': 'P',     'ACG': 'T',     'GCG': 'A',
    'TAT': 'Y',     'CAT': 'H',     'AAT': 'N',     'GAT': 'D',
    'TAC': 'Y',     'CAC': 'H',     'AAC': 'N',     'GAC': 'D',
    'TAA': 'Stop',  'CAA': 'Q',     'AAA': 'K',     'GAA': 'E',
    'TAG': 'Stop',  'CAG': 'Q',     'AAG': 'K',     'GAG': 'E',
    'TGT': 'C',     'CGT': 'R',     'AGT': 'S',     'GGT': 'G',
    'TGC': 'C',     'CGC': 'R',     'AGC': 'S',     'GGC': 'G',
    'TGA': 'Stop',  'CGA': 'R',     'AGA': 'R',     'GGA': 'G',
    'TGG': 'W',     'CGG': 'R',     'AGG': 'R',     'GGG': 'G'
}


def translate_codon(codon):
    protein = None
    if len(codon) == 3 and codon in DNA_CODON_TABLE:
        protein = DNA_CODON_TABLE[codon]
    return protein


def reverse_complement(dna):
    lookup = {'A':'T', 'T':'A', 'G':'C', 'C':'G'}
    rev_str= ''
    for c in reversed(dna):
        if c in "ACTG":
            rev_str+= lookup[c]
        else: rev_str+= "N"
    return rev_str


def possible_protein_strings(sequence,minimum_size):
    '''Find all possible translated products in a DNA sequence. 
       Save DNA seq of translated product
    '''
    
    results = []
    start_indices = []
    
    #capture indices of Methionine (start) codon
    seq_len = len(sequence)
    for i in range(seq_len):
        protein = translate_codon(sequence[i:i+3])
        if protein and protein == 'M':
            start_indices.append(i)

    #Build protein sequence from start codon index to "stop"...
    #...Discard peptides under constructuing when ambiuity base is hit.
    for start in start_indices:
        found_stop = False
        protein_string = ''

        for j in range(start, seq_len, 3):
            protein = translate_codon(sequence[j:j+3]) #capture amino
            

            if protein == 'Stop':  #stop building
                found_stop = True
                break
            
            if not protein:
                break
            
            protein_string += protein
        
        #Store peptides longer than user-set threshold
        if found_stop and len(protein_string) > minimum_size:
            results.append(protein_string)

    return results

import sys
from Bio import SeqIO
record_list  = list(SeqIO.parse(sys.argv[1], "fasta"))




# In[35]:


#Store peptide sequences for each DNA sequence as individual lists.
#...Concatinate list of these lists.

protein_list = []
for i in range(0, len(record_list)):
    forward = possible_protein_strings(record_list[i].seq,  50)                    #SET MINIMUM SIZE = 50
    reverse = possible_protein_strings(reverse_complement(record_list[i].seq),50)       #SET MINIMUM SIZE =50)
    as_list = [record_list[i].id] + forward + reverse
    protein_list.append(as_list)



# In[ ]:



#Write protein list to file with index identical to translated...
#...DNA index in first column. Peptides from same index stored...
#...as adjacent rows.
import csv
with open(sys.argv[2], 'w') as f:
    writer = csv.writer(f)
    for row in protein_list:
        writer.writerow(row)

