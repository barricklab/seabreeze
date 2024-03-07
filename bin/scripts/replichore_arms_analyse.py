#!/usr/bin/env python3
#This script accepts an ori sequence, dif sequence and a folder with fasta files of genome assemblies
#returns a csv file that has the length of replichore arms  and ratio of long arm:short arm

''' imports '''

import os
from collections import deque
import pandas as pd
import numpy as np
import argparse

''' fetch arguments from the command line '''

parser = argparse.ArgumentParser(description='replichore_arms_analyse.py, a script to compute lenght of replichore arms and degree of imbalance for genome assemblies')
parser.add_argument('--assemblies', help=' path to folder containing assemblies reindexed to the oric sequence')
parser.add_argument('--oric', help='string representing oric sequence')
parser.add_argument('--ancestor', help='fasta file of the ancestor')
parser.add_argument('--dif', help='string representing dif sequence')
parser.add_argument('--output', help='output filename')
parser.add_argument('--verbose', '-v', action='store_true', help='Enable verbose mode')
parser.add_argument('--noarms', action='store_true', help='Enable replichore arm size analysis. If false, only ori and dif locations reported')


def get_fasta_names(assemblies): 
    
    ''' returns a list of all the fasta files in the assemblies folder '''
    
    fasta_names = []
    os.chdir(assemblies)
    for entry in os.scandir():
        if entry.is_file() and entry.name.endswith(".fasta") and len(entry.name) > len(".fasta"): #check it is a fasta file which has a file name 
            fasta_names.append(entry.name)
    return fasta_names

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

def main(assemblies, oric, dif, ancestor, output, noarms):

    fasta_files=get_fasta_names(assemblies) #list with names of all fasta files

    if noarms:

        ''' only report the ori and dif locations '''

        df=pd.DataFrame(np.nan, index=range(len(fasta_files)), columns=["clone", "ori","dif","length"]) #stores final table
        row_idx=0
        for name in fasta_files:
            print(name)
            coord_dict=coords(name,oric,dif)
            df.loc[row_idx,"clone"]=name
            df['clone'] = df['clone'].str.replace('.fasta', '') #remove.fasta extension
            df.loc[row_idx,'ori']=coord_dict["oric_start"]
            df.loc[row_idx,'dif']=coord_dict["dif_start"]
            df.loc[row_idx,'length']=coord_dict["length"]
            row_idx+=1
        print(df)
        df.to_csv(output, sep='\t', index=False)

    if not(noarms):

        ''' perform replichore arm size analysis as well '''

        df=pd.DataFrame(np.nan, index=range(len(fasta_files)), columns=["clone", "ori","dif","length","arm_1", "arm_2", "ratio", "percent"]) #stores final table
        row_idx=0
        for name in fasta_files:
            print(name)
            coord_dict=coords(name,oric,dif)
            df.loc[row_idx,"clone"]=name
            df['clone'] = df['clone'].str.replace('.fasta', '') #remove.fasta extension
            df.loc[row_idx,'ori']=coord_dict["oric_start"]
            df.loc[row_idx,'dif']=coord_dict["dif_start"]
            df.loc[row_idx,'length']=coord_dict["length"]
            df.loc[row_idx,"arm_1"]=int(coord_dict["dif_start"])-int(coord_dict["oric_start"])
            df.loc[row_idx,"arm_2"]=int(coord_dict["length"])-int(coord_dict["dif_start"])
            # what the hell is this verbose ass code STOP
            # if df.loc[row_idx,"arm_1"] >= df.loc[row_idx,"arm_2"]:
            #     df.loc[row_idx,"ratio"]= df.loc[row_idx,"arm_1"]/df.loc[row_idx,"arm_2"]
            # if df.loc[row_idx,"arm_1"] < df.loc[row_idx,"arm_2"]:
            #     df.loc[row_idx,"ratio"]= df.loc[row_idx,"arm_2"]/df.loc[row_idx,"arm_1"]
            # df.loc[row_idx,'fraction']=abs(0.5-(coord_dict["dif_start"]/coord_dict["length"]))/0.5
            df.loc[row_idx,"ratio"]= (max(df.loc[row_idx,"arm_1"],df.loc[row_idx,"arm_2"]))/(min(df.loc[row_idx,"arm_1"],df.loc[row_idx,"arm_2"])) # ratio of long arm to short arm
            df.loc[row_idx,"percent"]= ((max(df.loc[row_idx,"arm_1"],df.loc[row_idx,"arm_2"]))/df.loc[row_idx,"length"])*100 # percent of the total length that the long arm is
            row_idx+=1
        print(df)
        # index_of_anc = df[df['clone'] == ancestor].index[0] #row index of ancestor
        # for row_idx in range(len(df)):
        #     df.loc[row_idx,"change"]= ((df.loc[row_idx,"ratio"]-df.loc[index_of_anc,"ratio"])/df.loc[index_of_anc,"ratio"])*100
        # print(df)
        df.to_csv(output, sep='\t', index=False)

if __name__ == '__main__':
    args = parser.parse_args()
    main(args.assemblies, args.oric, args.dif, args.ancestor, args.output, args.noarms)