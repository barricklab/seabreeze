#!/usr/bin/env python3
# this script accepts a folder with the clone_boundaries.tsv files
# a returns  single csv file that has all the clones, and the count of how many deletions
# per type of mechanism

''' imports '''

import os
import pandas as pd
pd.options.mode.chained_assignment = None  # default='warn'
import numpy as np
import argparse
from collections import Counter
import subprocess
from Bio import SeqIO

''' fetch arguments from the command line '''

parser = argparse.ArgumentParser(description='classify_deletions.py, a script to determine the mechanism of rearrangements')
parser.add_argument('--folder', help='name of the folder with the boundaries_tsv files')
# parser.add_argument('--folder_fasta', help='name of the folder with the assembly fasta files')
parser.add_argument('--output', help='name of output folder')
parser.add_argument('--verbose', '-v', action='store_true', help='Enable verbose mode')
parser.add_argument('--deletion',action='store_true', help='Only report deletions')
parser.add_argument('--inversion',action='store_true', help='Only report inversions')

def get_csv_names(folder):

    ''' return a list of all the tsv files in the folder '''

    os.chdir(folder)
    csv_names = []
    for entry in os.scandir():
        if entry.is_file() and entry.name.endswith("boundaries.csv") and len(entry.name) > len(".csv"):
            csv_names.append(entry.name)
    return csv_names


# def extract_sequence(start,stop,fasta):

#     with open(fasta, "r") as fasta_file:
#         records = list(SeqIO.parse(fasta_file, "fasta"))

#         if len(records) == 0:
#             print("Error: No sequences found in the FASTA file.")
#             return None

#         # Get the first sequence. Assuming only one sequence
#         sequence = str(records[0].seq)

#         # Extract the substring from start to end
#         subseq = sequence[start - 1:stop]

#         return subseq



def classify_deletions(file,inversion,deletion):

    ''' returns a dict with the name of clone, and count of each type of IS'''

    df=pd.read_csv(file)
    print(df)
    clone_name=file.replace('_boundaries.csv','')

    if deletion:
        print(f"Counting deletions")
        summary_dict={'clone':'', 'total':0,'between_IS':0, 'IS_mediated':0, 'other':0}
        summary_dict['clone']=file.replace('_boundaries.tsv','')
        df_del = df[df.loc[:,'tag_3'].str.contains('DEL')]
        df_del['Mechanism']=['']*(len(df_del)) # stores the mechanism of that SV
        df_del['Evidence']=['']*(len(df_del)) # stores the mechanism of that SV
        df_del.reset_index(drop=True, inplace=True)
        summary_dict['total']=len(df_del)
        for row_idx in range(len(df_del)):

            if (df_del.loc[row_idx,'L_ref'] == df_del.loc[row_idx,'R_ref'] == df_del.loc[row_idx,'L_query']) and pd.notna(df_del.loc[row_idx, 'L_ref']):
                summary_dict['between_IS']+=1
                df_del.loc[row_idx,"Mechanism"]='betweeen_IS'
                df_del.loc[row_idx,"Evidence"]='full'
                continue

            if df_del.loc[row_idx,'L_ref'] == df_del.loc[row_idx,'R_ref'] and pd.notna(df_del.loc[row_idx,'L_ref']):
                summary_dict['between_IS']+=1
                df_del.loc[row_idx,"Mechanism"]='betweeen_IS'
                df_del.loc[row_idx,"Evidence"]='incomplete'
                continue

            if (df_del.loc[row_idx,'L_ref'] == df_del.loc[row_idx,'L_query'] or df_del.loc[row_idx,'R_ref'] == df_del.loc[row_idx,'L_query']) and pd.notna(df_del.loc[row_idx,'L_query']):
                summary_dict['IS_mediated']+=1
                df_del.loc[row_idx,"Mechanism"]='IS_mediated'
                df_del.loc[row_idx,"Evidence"]='full'
                continue

            if pd.notna(df_del.loc[row_idx,'L_ref']) or pd.notna(df_del.loc[row_idx,'R_ref']):
                summary_dict['IS_mediated']+=1
                df_del.loc[row_idx,"Mechanism"]='IS_mediated'
                df_del.loc[row_idx,"Evidence"]='incomplete'
                continue

            if pd.notna(df_del.loc[row_idx,'L_query']):
                summary_dict['IS_mediated']+=1
                df_del.loc[row_idx,"Mechanism"]='IS_mediated'
                df_del.loc[row_idx,"Evidence"]='evolved'
                continue

            else:
                summary_dict['other']+=1
                df_del.loc[row_idx,"Mechanism"]='other'
                df_del.loc[row_idx,"Evidence"]='NA'

               # path_to_clone =
        print(summary_dict)
        out_filename= clone_name + "_deletion.csv"
        #columns_to_drop=["tag_1","tag_2","tag_3","tag_4"]

        #df_del=df_del.drop(columns=columns_to_drop)

        df_del.to_csv(out_filename, index=False,float_format='%.0f')

    if inversion:
        print (f"Counting inversions")
        summary_dict={'clone':'', 'total':0,'between_IS':0, 'IS_mediated':0, 'other':0}
        summary_dict['clone']=file.replace('_boundaries.tsv','')
        df_inv = df[df.loc[:,'tag_3']=='INV']
        df_inv['Mechanism']=['']*(len(df_inv)) # stores the mechanism of that SV
        df_inv['Evidence']=['']*(len(df_inv)) # stores the mechanism of that SV
        df_inv.reset_index(drop=True, inplace=True)
        summary_dict['total']=len(df_inv)
        for row_idx in range(len(df_inv)):

            if (df_inv.loc[row_idx,'L_ref'] == df_inv.loc[row_idx,'R_ref'] == df_inv.loc[row_idx,'L_query'] == df_inv.loc[row_idx,'R_query']) and pd.notna(df_inv.loc[row_idx,'L_ref']):
                summary_dict['between_IS']+=1
                df_inv.loc[row_idx,"Mechanism"]='betweeen_IS'
                df_inv.loc[row_idx,"Evidence"]='full'
                continue

            if ((df_inv.loc[row_idx,'L_ref'] == df_inv.loc[row_idx,'L_query'] == df_inv.loc[row_idx,'R_query']) or  (df_inv.loc[row_idx,'R_ref'] == df_inv.loc[row_idx,'L_query'] == df_inv.loc[row_idx,'R_query'])) and pd.notna(df_inv.loc[row_idx,'L_query']):
                summary_dict['IS_mediated']+=1
                df_inv.loc[row_idx,"Mechanism"]='IS_mediated'
                df_inv.loc[row_idx,"Evidence"]='full'
                continue

            if df_inv.loc[row_idx,'L_query'] == df_inv.loc[row_idx,'R_query'] and pd.notna(df_inv.loc[row_idx,'L_query']):
                summary_dict['IS_mediated']+=1
                df_inv.loc[row_idx,"Mechanism"]='IS_mediated'
                df_inv.loc[row_idx,"Evidence"]='evolved'
                continue

            if df_inv.loc[row_idx,'L_query'] == df_inv.loc[row_idx,'R_query'] and pd.notna(df_inv.loc[row_idx,'L_query']):
                summary_dict['IS_mediated']+=1
                df_inv.loc[row_idx,"Mechanism"]='IS_mediated'
                df_inv.loc[row_idx,"Evidence"]='evolved'
                continue

            if (df_inv.loc[row_idx,'L_ref'] == df_inv.loc[row_idx,'R_ref']):
                summary_dict['between_IS']+=1
                df_inv.loc[row_idx,"Mechanism"]='between_IS'
                df_inv.loc[row_idx,"Evidence"]='incomplete'
                continue

            else:
                summary_dict['other']+=1
                df_inv.loc[row_idx,"Mechanism"]='other'
                df_inv.loc[row_idx,"Evidence"]='NA'

        out_filename= clone_name + "_inversion.csv"
        print(summary_dict)
        #columns_to_drop=["tag_1","tag_2","tag_3","tag_4"]

        #df_inv=df_inv.drop(columns=columns_to_drop)

        df_inv.to_csv(out_filename, index=False,float_format='%.0f')

    return summary_dict


def main(folder,output,inversion,deletion):
    csv_names=get_csv_names(folder)
    print(csv_names)
    print("test")
    #df_summary=pd.DataFrame(np.nan, index=range(len(csv_names)), columns=['clone','between_IS', 'IS_mediated','other'])
    summary_list=[]
    for file in csv_names:
        print (file)
        summary_dict=classify_deletions(file,inversion,deletion)
        summary_list.append(summary_dict)
    df=pd.DataFrame(summary_list)
    df.to_csv(output, index=False,float_format='%.0f')
    #print(f"output name is {output}")
    # if os.path.isfile(output):
    #     print(f"The file '{output}' exists in the current working directory.")
    # else:
    #     print(f"The file '{output}' does not exist in the current working directory.")
    # print(df)


if __name__ == '__main__':
    args = parser.parse_args()
    main(args.folder,args.output,args.inversion,args.deletion)
