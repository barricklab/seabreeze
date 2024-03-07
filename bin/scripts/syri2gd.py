#!/usr/bin/env python3
#Author: ira-zibbu
#run from 06_Syri_Output
#This script converts all the deletions and inversions in a syri.out file to lines in genome diff format used by breseq

''' imports '''

import os
import pandas as pd
import numpy as np
import argparse

''' fetch arguments '''

parser = argparse.ArgumentParser(description='syri2gd.py, a script to convert a syri.out file to genome diff format used by breseq')
parser.add_argument('--syri', help='syri.out file')
parser.add_argument('--output', help='name of output gd file')
parser.add_argument('--inversion', action='store_true', help='include inversions in gd file')
parser.add_argument('--deletion', action='store_true', help='include deletions in gd file')
parser.add_argument('--verbose', '-v', action='store_true', help='Enable verbose mode')


def syri2gd(syri,output,inversion,deletion):
    
    df_syri = pd.read_table(syri,delimiter='\t', header=None) 
    chr_id=df_syri.iloc[0,0] #reference chromosome ID/ contig ID
    df_syri_del = df_syri[df_syri.iloc[:,8].str.contains('DEL')]
    df_syri_inv = df_syri[df_syri.iloc[:,8].str.contains('INV')]
    df_gd_del=pd.DataFrame(np.nan, index=range(len(df_syri_del)), columns=range(6))
    df_gd_inv=pd.DataFrame(np.nan, index=range(len(df_syri_inv)), columns=range(6))

    header_line="#=GENOME_DIFF\t1.0"
    with open(output, 'w') as file:
        file.write(header_line + '\n')  # Write the comment header line
    row_idx = 0

    if deletion: # add genome diff lines for deletions

        for row_idx in range(len(df_syri_del)):
            df_gd_del.iloc[row_idx,0]="DEL"
            df_gd_del.iloc[row_idx,1]=str(row_idx+1)
            df_gd_del.iloc[row_idx,2]="."
            df_gd_del.iloc[row_idx,3]=chr_id
            start_coords=df_syri_del.iloc[row_idx,1]
            end_coords=df_syri_del.iloc[row_idx,2]
            del_length=int(end_coords)-int(start_coords)
            df_gd_del.iloc[row_idx,4]=str(start_coords)
            df_gd_del.iloc[row_idx,5]=str(del_length)

        with open(output, 'a') as file:
            df_gd_del.to_csv(file, sep='\t', index=False, header=False)
    
    if inversion: # add genome diff lines for inversions

        index_list=range(row_idx+2,(len(df_syri_inv)+row_idx+2),1) # to continue the indices from the gd file
        for row_idx in range(len(df_syri_inv)):
            df_gd_inv.iloc[row_idx,0]="INV"
            df_gd_inv.iloc[row_idx,1]=str(index_list[row_idx])
            df_gd_inv.iloc[row_idx,2]="."
            df_gd_inv.iloc[row_idx,3]=chr_id
            start_coords=df_syri_inv.iloc[row_idx,1]
            end_coords=df_syri_inv.iloc[row_idx,2]
            inv_length=int(end_coords)-int(start_coords)
            df_gd_inv.iloc[row_idx,4]=str(start_coords)
            df_gd_inv.iloc[row_idx,5]=str(inv_length)

        with open(output, 'a') as file:
            df_gd_inv.to_csv(file, sep='\t', index=False, header=False)
        

def main(syri,output,inversion,deletion):

    syri2gd(syri,output,inversion,deletion)

if __name__ == '__main__':
    args = parser.parse_args()
    main(args.syri,args.output,args.inversion,args.deletion)

