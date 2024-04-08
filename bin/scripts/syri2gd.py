#!/usr/bin/env python3
#Author: ira-zibbu
#run from 06_Syri_Output
#This script converts all the deletions, inversions and duplications in a syri.out file to lines in genome diff format used by breseq

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
parser.add_argument('--amplification', action='store_true', help='include amplifications in gd file')
parser.add_argument('--verbose', '-v', action='store_true', help='Enable verbose mode')


def syri2gd(syri,output,inversion,deletion,amplification):
    
    df_syri = pd.read_table(syri,delimiter='\t', header=None) 
    chr_id=df_syri.iloc[0,0] #reference chromosome ID/ contig ID
    df_syri_del = df_syri[df_syri.iloc[:,10].str.contains('DEL')]
    # df_syri_inv = df_syri[df_syri.iloc[:,10.str.contains('INV')]
    # df_syri_amp = df_syri[df_syri.iloc[:,10].str.contains('DUP|INVDP')]
    df_syri_inv = df_syri[df_syri.iloc[:, 10] == 'INV']
    df_syri_amp = df_syri[df_syri.iloc[:, 10].isin(['DUP', 'INVDP','CPG'])]


    df_syri_amp = df_syri_amp.reset_index(drop=True)

    df_gd_del=pd.DataFrame(np.nan, index=range(len(df_syri_del)), columns=range(6))
    df_gd_inv=pd.DataFrame(np.nan, index=range(len(df_syri_inv)), columns=range(6))
    df_gd_amp=pd.DataFrame(np.nan, index=range(len(df_syri_amp)), columns=range(7))


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
        
    if amplification: # add genome diff lines for amplfications

        flag = any(df_syri.iloc[:, 9].str.contains('DUP|INVDP|CPG', na=False)) # check to see if there are any amplifications. Skip if none

        flag = flag if flag else False #true is there an an amplification
        
        if flag:

            index_list=range(index_list[-1]+1,(len(df_syri_amp)+index_list[-1])+1,1) # to continue the indices from the gd file
            print(index_list)
            #index_list=range(row_idx+2,(len(df_syri_amp)+row_idx+2),1) # to continue the indices from the gd file
            for row_idx in range(len(df_syri_amp)):
                df_gd_amp.iloc[row_idx,0]="AMP"
                df_gd_amp.iloc[row_idx,1]=str(index_list[row_idx])
                df_gd_amp.iloc[row_idx,2]="."
                df_gd_amp.iloc[row_idx,3]=chr_id
                start_coords=min(df_syri_amp.iloc[row_idx,1],df_syri_amp.iloc[row_idx,2]) # ensuring that the start coords are always smaller
                end_coords=max(df_syri_amp.iloc[row_idx,1],df_syri_amp.iloc[row_idx,2])
                amp_length=int(end_coords)-int(start_coords)
                df_gd_amp.iloc[row_idx,4]=str(start_coords)
                df_gd_amp.iloc[row_idx,5]=str(amp_length)
                df_gd_amp.iloc[row_idx,6]=str(2) # SyRI only reports duplications so this copy number will be 2 always


            with open(output, 'a') as file:
                df_gd_amp.to_csv(file, sep='\t', index=False, header=False)

def main(syri,output,inversion,deletion,amplification):

    syri2gd(syri,output,inversion,deletion,amplification)

if __name__ == '__main__':
    args = parser.parse_args()
    main(args.syri,args.output,args.inversion,args.deletion,args.amplification)

