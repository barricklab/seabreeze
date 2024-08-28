#!/usr/bin/env python3
#this script accepts a folder with the csv files produced by ISescan, a returns a single 
# summary csv file which contains the count of each type of IS element. 
# This also switches the names from the ISesca nomenclature to that of REL606

''' imports '''

import os
import pandas as pd
import numpy as np
import argparse
from collections import Counter

''' fetch arguments from the command line '''

parser = argparse.ArgumentParser(description='ISescan_summary.py, a script to return a summary of all of the IS elements')
parser.add_argument('--isescan', help=' path to folder containing ISescan csv files')
parser.add_argument('--output', help='output filename')
parser.add_argument('--ancestor', help='name of the ancestor relative to which copy number changes are computed')
parser.add_argument('--verbose', '-v', action='store_true', help='Enable verbose mode')

def get_csv_names(isescan):

    ''' return a list of all the csv files in the folder '''
    os.chdir(isescan)
    csv_names = []
    for entry in os.scandir():
        if entry.is_file() and entry.name.endswith(".csv") and len(entry.name) > len(".csv"):
            csv_names.append(entry.name)
    return csv_names


def fetch_isescan(filename):

    ''' returns a dataframe of ISescan table '''
    
    df = pd.read_csv(filename)
    column_to_list='cluster'
    cluster_names = df[column_to_list].tolist()
    return cluster_names

def IS_summary(csv_names):

    '''return a dataframe that has a list of clones and the number of each IS element '''

    IS_dict={"IS1_316":"IS1","IS3_168":"IS3","IS3_61":"IS150","IS4_107":"IS4","IS4_169":"IS186"}
    IS_list = list(IS_dict.values())
    IS_list.insert(0, "clone")
    df = pd.DataFrame(np.nan, index=range(len(csv_names)), columns=IS_list) #stores the summary dataframe
    
    row_idx=0
    for file in csv_names:
        cluster_list = fetch_isescan(file) #list of IS elements
        renamed_IS_list = [IS_dict[key] for key in cluster_list if key in IS_dict]
        IS_count = dict(Counter(renamed_IS_list))
        df.loc[row_idx,'clone']=file.replace('.csv','')

        for key, value in IS_count.items():
            if key in df.columns:
                df.loc[row_idx,key] = value
        row_idx+=1

    return df

def IS_summary_relative(df_summary,ancestor):

    ''' return a new dataframe which contains the copy number change relative to the ancestor'''

    reference_row = df_summary[df_summary['clone'] == ancestor]
    df_relative = df_summary.copy()  # Make a copy of the original DataFrame
    df_relative.iloc[:, 1:] -= reference_row.iloc[:, 1:].values
    return df_relative


def main(isescan,output, ancestor):
    csv_names = get_csv_names(isescan)
    df_summary = IS_summary(csv_names)
    df_relative=IS_summary_relative(df_summary,ancestor)

    df_summary.to_csv(output+".csv", sep=',', index=False,  float_format='%.0f')
    df_relative.to_csv(output+"_copy_change.csv", sep=',', index=False, float_format='%.0f')



if __name__ == '__main__':
    args = parser.parse_args()
    main(args.isescan, args.output,args.ancestor)