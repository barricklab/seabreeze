#!/usr/bin/env python3

''' imports '''

import os
from collections import deque
import pandas as pd
import numpy as np
import argparse

''' fetch arguments from the command line '''

parser = argparse.ArgumentParser(description='replichore_arms_analyse.py, a script to compute lenght of replichore arms and degree of imbalance for genome assemblies')
parser.add_argument('--genomes', help='path to folder containing assemblies')
parser.add_argument('--data', help='path to the data.csv file')
parser.add_argument('--sequences', help='path to the csv file with the sequences for the ori and dif loci')
parser.add_argument('--output', help='output filename')
parser.add_argument('--verbose', '-v', action='store_true', help='Enable verbose mode')
parser.add_argument('--noarms', action='store_true', help='Enable replichore arm size analysis. If false, only ori and dif locations reported')


def reverse_complement(string):

    ''' Return the reverse complement of a string '''

    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N', # upper case typical
                  'a': 't', 'c': 'g', 'g': 'c', 't': 'a', 'n': "n", # lower case typical
                  'W': 'W', 'S': 'S', 'M': 'K', 'K': 'M', 'Y': 'R', 'R': 'Y', # upper case rare pt1
                  'w': 'w', 's': 's', 'm': 'k', 'k': 'm', 'y': 'r', 'r': 'y', # lower case rare pt1
                  'B': 'V', 'D': 'H', 'H': 'D', 'V': 'B', # upper case rare pt2,
                  'b': 'v', 'd': 'h', 'h': 'd', 'v': 'b', # lower case rare pt2
                  '-': '-',}

    if "u" in string or "U" in string:
        complement["u"] = "a"
        complement["U"] = "A"
        complement["A"] = "U"
        complement["a"] = "u"
        del complement["T"]
        del complement["t"]

    return ''.join([complement.get(base, base) for base in reversed(string)])

def coords(fasta_name, oric, dif):

    ''' return start coords of ori and dif sequence specified '''

    fasta = deque() #stores the lines of the fasta files
    contig_count=0 #stores the number of contigs in the fasta file
    with open(fasta_name, "r") as file:
        for line in file:
            line = line.strip()
            if line[0] == ">":
                contig_count+=1
            else:
                fasta.append(line)
    if contig_count > 1: #skipping analysis if more than one contig
        print("more that 1 contig, skipping ", fasta_name)
        return 0
    sequence=''.join(fasta) #stores the entire assembly as a string
    length=len(sequence)
    oric_start_coords=0
    dif_start_coords=0

    coord_dict={"oric_start":oric_start_coords, "dif_start":dif_start_coords, "length":len(sequence)} #stores the start coords

    #try forward sequence

    oric_start_coords = sequence.find(oric)
    if oric_start_coords != -1: #this is false if the oric is not found in the sequence
        coord_dict["oric_start"]=oric_start_coords


    dif_start_coords = sequence.find(dif)
    if dif_start_coords != -1: #this is false if the dif is not found in the sequence
        coord_dict["dif_start"]=dif_start_coords

    #try reverse complement

    sequence_rc=reverse_complement(sequence)

    oric_start_coords = sequence_rc.find(oric)
    if oric_start_coords != -1: #this is false if the oric is not found in the sequence
        coord_dict["oric_start"]=length-oric_start_coords

    dif_start_coords = sequence_rc.find(dif)
    if dif_start_coords != -1: #this is false if the dif is not found in the sequence
        coord_dict["dif_start"]=length-dif_start_coords

    return coord_dict


def get_sequence(sequences_df, name,seq_type):

    ''' Accept the sequences_df dataframe, and return the ori/dif sequences (depending on the value supplied on seq_type) for a given ancestor '''

    if name in sequences_df['ancestor'].values:
        sequence = sequences_df.loc[sequences_df['ancestor'] == name, seq_type].iloc[0]
        return sequence
    else:
        return None  # Or raise an exception, depending on the use case


def main(genomes, data, sequences, output, noarms):

    from fasta_stats import read_data # import statement is here because otherwise is causes trouble with pytest not finding this module


    data_df=read_data(data) # table of genome assemblies and their ancestor
    sequences_df=read_data(sequences) # table of ancestors and their ori and dif sequences
    # fasta_files=get_fasta_names(genomes) # list all of the fasta file in the genomes directory
    assembly_to_ancestor_dict = dict(zip(data_df['assembly'], data_df['ancestor']))


    assemblies=data_df['assembly'].tolist() # a list of all of the assemblies specified in the data.csv file

    if noarms:

        ''' only report the ori and dif locations '''

        df=pd.DataFrame(np.nan, index=range(len(assemblies)), columns=["clone", "ori","dif","length"]) #stores final table
        row_idx=0
        for name in assemblies:
            print(name)
            file_name=f"{name}.fasta"
            file_path=os.path.join(genomes,file_name)
            oric=get_sequence(sequences_df, assembly_to_ancestor_dict[name],'ori')
            dif=get_sequence(sequences_df, assembly_to_ancestor_dict[name],'dif')
            coord_dict=coords(file_path,oric,dif)
            df.loc[row_idx,"clone"]=name
            df.loc[row_idx,'ori']=coord_dict["oric_start"]
            df.loc[row_idx,'dif']=coord_dict["dif_start"]
            df.loc[row_idx,'length']=coord_dict["length"]
            row_idx+=1
        print(df)
        df.to_csv(output, index=False)

    if not(noarms):

        ''' perform replichore arm size analysis as well '''

        df=pd.DataFrame(np.nan, index=range(len(assemblies)), columns=["clone", "ori","dif","length","arm_1", "arm_2", "ratio", "percent"]) #stores final table
        row_idx=0
        for name in assemblies:

            file_name=f"{name}.fasta"
            file_path=os.path.join(genomes,file_name)
            oric=get_sequence(sequences_df, assembly_to_ancestor_dict[name],'ori')
            dif=get_sequence(sequences_df, assembly_to_ancestor_dict[name],'dif')
            coord_dict=coords(file_path,oric,dif)

            df.loc[row_idx,"clone"]=name
            df.loc[row_idx,'ori']=coord_dict["oric_start"]
            df.loc[row_idx,'dif']=coord_dict["dif_start"]
            df.loc[row_idx,'length']=coord_dict["length"]
            df.loc[row_idx,"arm_1"]=int(coord_dict["dif_start"])-int(coord_dict["oric_start"])
            df.loc[row_idx,"arm_2"]=int(coord_dict["length"])-int(coord_dict["dif_start"])

            df.loc[row_idx,"ratio"]= (max(df.loc[row_idx,"arm_1"],df.loc[row_idx,"arm_2"]))/(min(df.loc[row_idx,"arm_1"],df.loc[row_idx,"arm_2"])) # ratio of long arm to short arm
            df.loc[row_idx,"percent"]= ((max(df.loc[row_idx,"arm_1"],df.loc[row_idx,"arm_2"]))/df.loc[row_idx,"length"])*100 # percent of the total length that the long arm is
            row_idx+=1

        print(df)

        df.to_csv(output, index=False)

if __name__ == '__main__':
    args = parser.parse_args()
    main(args.genomes, args.data, args.sequences, args.output, args.noarms)
