#!/usr/bin/env python3

# Get full alpha and beta TCR sequences when providing variable gene, cdr3, 
# and j gene to be used in structural models and TCR cloning
# Run by calling ./chain_seqs.py chain_example.csv

from Bio import SeqIO
from pprint import pprint
import difflib
import pandas as pd
import sys
import os

def get_seqs(fasta_filename):
    with open(fasta_filename, 'r') as f:
        sequences = SeqIO.parse(f, 'fasta')
        fasta_with_seqs = {}
        for sequence in sequences: 
            # get data from each record by index
            vals = sequence.description.split('|')
            # do not add any rows with pseudogenes
            if 'P' in vals[3]:   
                continue 
            
            aa_sequence = str(sequence.seq) 
            fasta_with_seqs[vals[1]] = aa_sequence
            
    return fasta_with_seqs

def assemble(str_list, min=2, max=15):
    if len(str_list) < 2:
        return set(str_list)
    output = set()
    string = str_list.pop()
    for i, candidate in enumerate(str_list):
        matches = set()
        if candidate in string:
            matches.add(string)
        elif string in candidate:
            matches.add(candidate)
        for n in range(min, max + 1):
            if candidate[:n] == string[-n:]:
                matches.add(string + candidate[n:])
            if candidate[-n:] == string[:n]:
                matches.add(candidate[:-n] + string)
        for match in matches:
            output.update(assemble(str_list[:i] + str_list[i + 1:] + [match]))
    return output

def get_overlap(s1, s2):
    s = difflib.SequenceMatcher(None, s1, s2, autojunk=False)
    pos_a, pos_b, size = s.find_longest_match(0, len(s1), 0, len(s2))
    print(list(s.get_matching_blocks())) 
    print(s.find_longest_match(0, len(s1), 0, len(s2)))
    return s1[pos_a:pos_a+size]
    # want match between beginning of one seq and end of another 
    
def main():

    df = pd.read_csv(sys.argv[1])
    filename = sys.argv[1]

    # all fasta files created from IMGT
    # get variable gene amino acid sequence dictionary
    trvs = get_seqs('IMGT_human_trv_ref.txt')
    pprint(trvs)
    
    # full lpart1+lpart2 gene amino acid sequence dictionary
    l12 = get_seqs('IMGT_human_lpart1+lpart2.txt')

    # remove allele suffix from key and keep unique entries
    clean_l12 = {}
    for key,value in l12.items():
        if value not in clean_l12.values():
            newkey = key[:-3]
            print(newkey)
            clean_l12[newkey] = value
    pprint(clean_l12) 
    
    # get j gene amino acid sequence dictionary
    trjs = get_seqs('IMGT_human_trj_ref.txt') 
    pprint(trjs)
    
    # add column to df for final sequence
    df['chain'] = ''
    for row in df.itertuples():
        print("running", row.id)
        cdr3 = row.cdr3
        # look up sequences based on input gene names
        v_noallele = row.v[:-3]
        l = clean_l12[v_noallele]
        v = trvs[row.v]
        j = trjs[row.j]
       
        print("signal sequence:", l)
        print("v gene sequence:", v)
        print("cdr3 sequence:", cdr3)
        print("j gene sequence:", j)
 
        # add signal sequence and V gene sequence
        # no overlap here, added side by side
        lv = l+v
        print('L-part1+Lpart2+V:', lv)
        
        try: 
            # get latest overlap of vgene end at beginning of cdr3
            # this is hard coded for the last 6 aas, may need adjusted             
            overlap = get_overlap(lv[-6:], cdr3) 
            print('overlap bw V and CDR3:', overlap)
            v_index = lv.rindex(overlap) + len(overlap)
            print('vindex:', v_index)
            # get end of V gene from CDR3
            vend = lv[:v_index]
            print('v end of chain:', vend)
            # merge V+CDR3 sequence (head to tail) with clean overlap
            vc = list(assemble([vend, cdr3]))  
            print('combined chain:', vc[0])

            # get earliest overlap of trj at end of cdr3            
            overlap = get_overlap(cdr3, j)
            print('overlap bw CDR3 and J:', overlap)
            j_index = j.rindex(overlap)
            # get end of J gene from CDR3
            end = j[j_index:]
            print('J end of chain:', end)
            # merge V+CDR3 sequence (head to tail) with clean J overlap
            vcj = list(assemble([vc[0], end])) 
            print('combined chain:', vcj[0])
            df.at[row.Index, 'chain'] = vcj[0]
            print("----------------------------------")
 
        except IndexError as e: 
            print('No overlap detected between segments, skipping TCR')
            print("----------------------------------") 
            continue # continue to next TCR if this error is encountered

        #except KeyError as e: 
        #    print('Missing reference, skipping TCR')
        #    continue # continue to next TCR if this error is encountered

    outfile = os.path.splitext(filename)[0] + "_out.csv"    
    df.to_csv(outfile)

main()
