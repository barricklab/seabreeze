#!/usr/bin/env python3
# This script takes a syri.out file and returns a tsv file
# that describes which SVs were flanked by IS elements, distance to nearest Is element.

''' imports '''
import pandas as pd
import numpy as np
import os
import argparse


''' fetch arguments '''
parser = argparse.ArgumentParser(description='IS_SV_border.py, a script to determine which SVs were mediated by IS elements')
parser.add_argument('--ancestor', help='path to isescan csv file of the ancestor')
parser.add_argument('--evolved', help='path to isescan csv file of the evolved clone')
parser.add_argument('--syri', help='path to syri.out file of the evolved clone')
parser.add_argument('--output', help='filename for .tsv output')
parser.add_argument('--verbose', '-v', action='store_true', help='Enable verbose mode')


def fetch_isescan(filename):

    ''' returns a dataframe of ISescan table '''

    columns_to_read=["cluster", "isBegin", "isEnd"]
    df = pd.read_table(filename, delimiter=',', usecols=columns_to_read)
    return df

def fetch_syri(filename):

    ''' return syri.out file as a dataframe, and add column headings '''

    df = pd.read_table(filename, delimiter='\t', header=None, dtype=str)
    df.columns=["ref_ID", "ref_start", "ref_stop", "seq_deleted", "seq_inserted", "query_ID", "query_start", "query_stop", "tag_1", "tag_2", "tag_3", "tag_4"] #adding columns for easy reference

    #additional columns which stores the name of IS elements and their distance from the boundary of the SV

    df['L_ref']=['NA']*(len(df))
    df['L_ref_distance']=[0]*(len(df))
    df['R_ref']=['NA']*(len(df))
    df['R_ref_distance']=[0]*(len(df))
    df['L_query']=['NA']*(len(df))
    df['L_query_distance']=[0]*(len(df))
    df['R_query']=['NA']*(len(df))
    df['R_query_distance']=[0]*(len(df))

    return df

def annotate_reference(df_syri,df_ancestor):

    '''annotate L-ref and R-ref and their distances of SVs based on coords of IS elements in ISescn '''

    threshold = 20 #the gap permissible between an IS element and SV
    for row_idx_syri in range(len(df_syri)):

        start_SV=df_syri.loc[row_idx_syri,'ref_start']
        stop_SV=df_syri.loc[row_idx_syri,'ref_stop']

        if start_SV == '-' or stop_SV == '-': #skipping SV that do not have reference coords listed
            continue

        for row_idx_anc in range(len(df_ancestor)):
            start_IS=df_ancestor.loc[row_idx_anc,'isBegin']
            stop_IS=df_ancestor.loc[row_idx_anc,'isEnd']
            name_IS=df_ancestor.loc[row_idx_anc,'cluster']
            #name_IS=convert_is_name(name_IS)

            left_difference_in = float(start_IS) - float(start_SV) # the IS element is inside the SV, near the start of SV
            right_difference_in = float(stop_SV) - float(stop_IS) # the IS element is inside the SV, near the end of SV
            left_difference_out = float(stop_IS) - float(start_SV) # the IS element is outside the SV, near the start of the SV
            right_difference_out = float(stop_SV) - float(start_IS) # the IS element is outside the SV, near the end of the SV

            if abs(left_difference_in)<threshold:
                df_syri.loc[row_idx_syri,"L_ref"] = name_IS
                df_syri.loc[row_idx_syri,"L_ref_distance"]=str(int(left_difference_in))

            if abs(right_difference_in)<threshold:
                df_syri.loc[row_idx_syri,"R_ref"] = name_IS
                df_syri.loc[row_idx_syri,"R_ref_distance"]=str(int(right_difference_in))

            if abs(left_difference_out)<threshold:
                df_syri.loc[row_idx_syri,"L_ref"] = name_IS
                df_syri.loc[row_idx_syri,"L_ref_distance"]=str(int(left_difference_out))

            if abs(right_difference_out)<threshold:
                df_syri.loc[row_idx_syri,"R_ref"] = name_IS
                df_syri.loc[row_idx_syri,"R_ref_distance"]=str(int(right_difference_out))

    return df_syri

def annotate_query(df_syri,df_evolved):

    '''annotate L-query and R-query of SVs based on coords of IS elements in ISescn of the evolved clone '''

    threshold = 20 #the gap permissible between an IS element and SV
    for row_idx_syri in range(len(df_syri)):

        start_SV=df_syri.loc[row_idx_syri,'query_start']
        stop_SV=df_syri.loc[row_idx_syri,'query_stop']

        if start_SV == '-' or stop_SV == '-': #skipping SV that do not have reference coords listed
            continue

        for row_idx_evol in range(len(df_evolved)):
            start_IS=df_evolved.loc[row_idx_evol,'isBegin']
            stop_IS=df_evolved.loc[row_idx_evol,'isEnd']
            name_IS=df_evolved.loc[row_idx_evol,'cluster']
            #name_IS=convert_is_name(name_IS)

            left_difference_in = float(start_IS) - float(start_SV) # the IS element is inside the SV, near the start of SV
            right_difference_in = float(stop_SV) - float(stop_IS) # the IS element is inside the SV, near the end of SV
            left_difference_out = float(stop_IS) - float(start_SV) # the IS element is outside the SV, near the start of the SV
            right_difference_out = float(stop_SV) - float(start_IS) # the IS element is outside the SV, near the end of the SV

            if abs(left_difference_in)<threshold:
                df_syri.loc[row_idx_syri,"L_query"] = name_IS
                df_syri.loc[row_idx_syri,"L_query_distance"]=str(int(left_difference_in))

            if abs(right_difference_in)<threshold:
                df_syri.loc[row_idx_syri,"R_query"] = name_IS
                df_syri.loc[row_idx_syri,"R_query_distance"]=str(int(right_difference_in))

            if abs(left_difference_out)<threshold:
                df_syri.loc[row_idx_syri,"L_query"] = name_IS
                df_syri.loc[row_idx_syri,"L_query_distance"]=str(int(left_difference_out))

            if abs(right_difference_out)<threshold:
                df_syri.loc[row_idx_syri,"R_query"] = name_IS
                df_syri.loc[row_idx_syri,"R_query_distance"]=str(int(right_difference_out))

    return df_syri

def write_syri(df_syri,output_filename):

    '''accepts a syri dataframe and writes out a tab separated file'''

    columns_to_write=["ref_start", "ref_stop","query_start", "query_stop", "tag_3", "L_ref", "L_ref_distance","R_ref", "R_ref_distance","L_query","L_query_distance" ,"R_query", "R_query_distance"]
    df_syri=df_syri[columns_to_write]
    df_syri.to_csv(output_filename, index=False,  float_format='%.0f')

def convert_is_name(isescan_name):
    # This function is not in use in this version of SyRI
    ''' Accept a name of an IS element with ISescan nomenclature, and return the corresponding name from REL606.gff3 annotation as from this table: https://app.box.com/file/1350490004900 '''

    IS_dict={"IS1_316":"IS1","IS3_168":"IS3","IS3_61":"IS150","IS4_107":"IS_4","IS4_169":"IS186"}
    if isescan_name in IS_dict:
        return IS_dict[isescan_name]
    else:
        return isescan_name


def main(ancestor, evolved, syri,output):
    df_ancestor=fetch_isescan(ancestor)
    df_evolved=fetch_isescan(evolved)
    df_syri=fetch_syri(syri)
    df_syri=annotate_reference(df_syri,df_ancestor)
    df_syri=annotate_query(df_syri,df_evolved)
    write_syri(df_syri,output)

if __name__ == '__main__':
    args = parser.parse_args()
    main(args.ancestor, args.evolved, args.syri, args.output)
