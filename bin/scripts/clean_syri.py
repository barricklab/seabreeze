#!/usr/bin/env python3
# Last update: 2024-02-29 [Now it uses ISEScan tables for both the the subject and query. This does not generate a different output vs using the gff3 file for the subject]
# Author: Ira Zibbu
# This script cleans up the syri.out files 

''' imports '''

import pandas as pd
import fileinput
import re
import numpy as np
import argparse
import os

''' fetch arguments '''

parser = argparse.ArgumentParser(description='clean_syri.py: clean up a syri.out file')
#parser.add_argument('--reference', help='Path to the folder with the reference gff3 file')
#parser.add_argument('--prefix', help='prefix for the reference gff3 file')
parser.add_argument('--syri', help='File name of syri.out file')
parser.add_argument('--isescan_subject', help='File name of isescan csv file of the subject (ancestor)')
parser.add_argument('--isescan_query', help='File name of isescan csv file of the subject (evolved clone)')
parser.add_argument('--verbose', '-v', action='store_true', help='Enable verbose mode')

# def fetch_repeat_regions(ref_filepath, prefix): 

# 	""" accepts reference .gff3 file and returns a gff3 file with only the lines corresponding to IS elements and no header """

# 	wd=os.getcwd()
# 	os.chdir(ref_filepath) #go to folder with reference genome assembly
# 	input_filename=prefix+".gff3"
# 	output_file = prefix+"_repeat.gff3"
# 	with fileinput.input(files=(input_filename,), inplace=False) as f:
# 		with open(output_file, "w") as output:
# 			for line in f:
# 				if re.search(r"\trepeat_region\t", line):
# 					output.write(line)
# 	#print("Got repeat lines")
# 	os.chdir(wd)

# def read_gff3(ref_filepath, prefix): 
	
# 	""" read gff3 file into a dataframe """

# 	wd=os.getcwd()
# 	os.chdir(ref_filepath) #go to folder with reference genome
# 	input_filename=prefix+"_repeat.gff3"
# 	df_gff3=pd.read_table(input_filename, delimiter='\t', header=None)
# 	gff3_annotation_col=df_gff3.iloc[:,8] #column with IS element names
# 	gff3_IS_list=gff3_annotation_col.tolist() 
# 	idx=0
# 	for name in gff3_IS_list: #only retaining the IS element names, and not other notes from the annotation column
# 		name=name.replace('Name=','')
# 		name=name.replace(";Note=repeat region",'')
# 		gff3_IS_list[idx]=name
# 		idx+=1
# 	df2 = pd.DataFrame({'start': df_gff3.iloc[:,3], 'stop': df_gff3.loc[:,4], 'IS_name': gff3_IS_list}) #dataframe that stores the start/ stop coords of Is elements and 
# 	#print("Read gff3 file")
# 	os.chdir(wd)
# 	return df2


def syri_remove_duplicate_lines(syri_filename): #

	""" remove 	"SYNAL" "INVAL" "TRANSAL" "INVTRAL" "DUPAL" "INVDPAL" tagged lines as they don't affect the synteny plot or further analysis """

	df_syri = pd.read_table(syri_filename,delimiter='\t', header=None, dtype=str)
	rows_to_drop=[] # stores row indices of lines to be removed
	remove_tags=["SYNAL", "INVAL" ,"TRANSAL", "INVTRAL" ,"DUPAL", "INVDPAL", "HDR", "TDM"] #tags to be removed
	rows_to_drop = df_syri[df_syri[8].str.contains('|'.join(remove_tags), na=False)].index
	df_syri = df_syri.drop(rows_to_drop).reset_index(drop=True)

	#print("Removed unecessary lines from syri.out file")
	return (df_syri)

# def syri_annotate_iselements_gff3(df_syri, df_gff3):

# 	""" use the gff3 annotations of IS elements to annotate which SVs contain IS elements by checking if the reference coords of the SV overlap with those of IS elements in the gff3 """

# 	df_syri.columns=["ref_ID", "ref_start", "ref_stop", "seq_deleted", "seq_inserted", "query_ID", "query_start", "query_stop", "tag_1", "tag_2", "tag_3", "tag_4"] #adding columns for easy reference
# 	df_syri_annotated=df_syri #store the syri dataframe with an additional column with binary values. 0 - the alignment is not an IS element and 1 - the alignment is an IS element
# 	df_syri_annotated['annotation_ref']=[0]*(len(df_syri)) #additional column which stores the annotation. Initialized to 0

# 	threshold=500 #minimum number of bases in an SV that should not overlap with an annotated IS element
# 	row_idx_syri=0 #row number counter for df_syri_annotated

# 	while row_idx_syri <= len(df_syri_annotated):

# 		if row_idx_syri == (len(df_syri_annotated)): #this code block is executed when the last line is encountered. It returns a df
# 				row_idx_syri+=1
# 				print("PROGRESS: Annotated IS elements of reference")
# 				return df_syri_annotated

# 		ref_start=df_syri_annotated.loc[row_idx_syri,'ref_start']
# 		ref_stop=df_syri_annotated.loc[row_idx_syri,'ref_stop']

# 		if ref_start != "-" and ref_stop != "-": #non-empty reference coords
# 			ref_start=int(ref_start)
# 			ref_stop=int(ref_stop)
# 			ref_length=ref_stop - ref_start
# 			row_idx_gff3=0

# 			is_count=0 # stores a count of the number of IS elements that overlap with the current SV
# 			is_len=0 # stores the sum of lengths of the IS element that overlaps with the current SV

# 			while row_idx_gff3 < len(df_gff3): #loop over lines of the gff3 Is elements
# 				IS_start=int(df_gff3.loc[row_idx_gff3,'start']) # start coords
# 				IS_stop=int(df_gff3.loc[row_idx_gff3,'stop']) # stop coords
# 				if ref_start >= IS_start and ref_start <= IS_stop and IS_stop <= ref_stop: #this case is when the IS element begins before the SV and extends into the SV
# 					is_count+=1
# 					is_len+=IS_stop-ref_start
# 					row_idx_gff3+=1
# 					continue
# 				if ref_start <= IS_start and IS_start <= ref_stop and IS_stop <= ref_stop: #this case is when the IS element is entirely contained within the SV
# 					is_count+=1
# 					is_len+=IS_stop-IS_start
# 					row_idx_gff3+=1
# 					continue			
# 				if ref_start <= IS_start and IS_start <= ref_stop and ref_stop <= IS_stop: #this case is when the IS element begins inside the SV and then extends beyond it
# 					is_count+=1
# 					is_len+=ref_stop-IS_start
# 					row_idx_gff3+=1
# 					continue
# 				if IS_start <= ref_start and IS_stop >= ref_stop:#this is when the IS element begins before and ends after the SV (i.e SV is contained within an IS element)
# 					is_count+=1
# 					is_len+=ref_stop-ref_start
# 					row_idx_gff3+=1
# 					continue
# 				else:
# 					row_idx_gff3+=1

# 			difference=ref_length-is_len # length of the SV that is not an IS element
			
# 			if ref_length > threshold and difference <= threshold:
# 				df_syri_annotated.loc[row_idx_syri,'annotation_ref']=1

# 		row_idx_syri+=1


def syri_annotate_query_iselements_isescan(df_syri_annotated, isescan_filename): 

	""" use the ISescan results ofthe evolved clone to annotate which SVs contain IS elements by checking if the query coords of the SV overlap with those of IS elements reported by ISescanf3 """

	columns_to_read=["cluster", "isBegin", "isEnd"]
	df_isescan = pd.read_table(isescan_filename, delimiter=',', usecols=columns_to_read)

	df_syri_annotated['annotation_query']=[0]*(len(df_syri_annotated)) #additional column which stores the annotation. Initialized to 0. It is 1 if the coords overlap with Is element
	threshold=500 #minimum number of bases in an SV that should not overlap with an annotated IS element
	row_idx_syri=0 #row number counter for df_syri_annotated


	while row_idx_syri <= len(df_syri_annotated):  

		if row_idx_syri == (len(df_syri_annotated)): #this code block is executed when the last line is encountered. It write out a tsv file
				row_idx_syri+=1
				print("PROGRESS: Annotated IS elements of query")

				return df_syri_annotated

		query_start=df_syri_annotated.loc[row_idx_syri,'query_start']
		query_stop=df_syri_annotated.loc[row_idx_syri,'query_stop']

		if query_start > query_stop: #reversing the coords if this is an inversion
			temp = query_start
			query_start=query_stop
			query_stop=temp

		if query_start != "-" and query_stop != "-": #non-empty query coords
			query_start=int(query_start)
			query_stop=int(query_stop)
			query_length=query_stop - query_start
			row_idx_isescan=0

			is_count=0 # stores a count of the number of IS elements that overlap with the current SV
			is_len=0 # stores the sum of lengths of the IS element that overlaps with the current SV

			while row_idx_isescan < len(df_isescan): #loop over lines of the isescan Is elements
				IS_start=df_isescan.loc[row_idx_isescan,'isBegin'] # start coords
				IS_stop=df_isescan.loc[row_idx_isescan,'isEnd'] # stop coords

				if query_start >= IS_start and query_start <= IS_stop and IS_stop <= query_stop:#this case is when the IS element begins before the SV and extends into the SV
					is_count+=1
					is_len+=IS_stop-query_start
					row_idx_isescan+=1
					continue

				if query_start <= IS_start and IS_start <= query_stop and IS_stop <= query_stop: #this case is when the IS element is entirely contained within the SV
					is_count+=1
					is_len+=IS_stop-IS_start
					row_idx_isescan+=1
					continue			

				if query_start <= IS_start and IS_start <= query_stop and query_stop <= IS_stop: #this case is when the IS element begins inside the SV and then extends beyond it
					is_count+=1
					is_len+=query_stop-IS_start
					row_idx_isescan+=1
					continue

				if IS_start <= query_start and IS_stop >= query_stop:#this is when the IS element begins before and ends after the SV (i.e SV is contained within an IS element)
					is_count+=1
					is_len+=query_stop-query_start
					row_idx_isescan+=1
					continue

				else:
					row_idx_isescan+=1

			difference=query_length-is_len # length of the SV that is not an IS element
			
			if query_length >= threshold and difference <= threshold:
				df_syri_annotated.loc[row_idx_syri,'annotation_query']=1

		row_idx_syri+=1

def syri_annotate_subject_iselements_isescan(df_syri, isescan_filename): 

	""" use the ISescan results ofthe evolved clone to annotate which SVs contain IS elements by checking if the ref coords of the SV overlap with those of IS elements reported by ISescanf3 """

	columns_to_read=["cluster", "isBegin", "isEnd"]
	df_isescan = pd.read_table(isescan_filename, delimiter=',', usecols=columns_to_read)

	df_syri.columns=["ref_ID", "ref_start", "ref_stop", "seq_deleted", "seq_inserted", "query_ID", "query_start", "query_stop", "tag_1", "tag_2", "tag_3", "tag_4"] #adding columns for easy reference
	df_syri_annotated=df_syri #store the syri dataframe with an additional column with binary values. 0 - the alignment is not an IS element and 1 - the alignment is an IS element

	df_syri_annotated['annotation_ref']=[0]*(len(df_syri_annotated)) #additional column which stores the annotation. Initialized to 0. It is 1 if the coords overlap with Is element
	threshold=500 #minimum number of bases in an SV that should not overlap with an annotated IS element
	row_idx_syri=0 #row number counter for df_syri_annotated


	while row_idx_syri <= len(df_syri_annotated):  

		if row_idx_syri == (len(df_syri_annotated)): #this code block is executed when the last line is encountered. It write out a tsv file
				row_idx_syri+=1
				print("PROGRESS: Annotated IS elements of ref")

				return df_syri_annotated

		ref_start=df_syri_annotated.loc[row_idx_syri,'ref_start']
		ref_stop=df_syri_annotated.loc[row_idx_syri,'ref_stop']

		if ref_start > ref_stop: #reversing the coords if this is an inversion
			temp = ref_start
			ref_start=ref_stop
			ref_stop=temp

		if ref_start != "-" and ref_stop != "-": #non-empty ref coords
			ref_start=int(ref_start)
			ref_stop=int(ref_stop)
			ref_length=ref_stop - ref_start
			row_idx_isescan=0

			is_count=0 # stores a count of the number of IS elements that overlap with the current SV
			is_len=0 # stores the sum of lengths of the IS element that overlaps with the current SV

			while row_idx_isescan < len(df_isescan): #loop over lines of the isescan Is elements
				IS_start=df_isescan.loc[row_idx_isescan,'isBegin'] # start coords
				IS_stop=df_isescan.loc[row_idx_isescan,'isEnd'] # stop coords

				if ref_start >= IS_start and ref_start <= IS_stop and IS_stop <= ref_stop:#this case is when the IS element begins before the SV and extends into the SV
					is_count+=1
					is_len+=IS_stop-ref_start
					row_idx_isescan+=1
					continue

				if ref_start <= IS_start and IS_start <= ref_stop and IS_stop <= ref_stop: #this case is when the IS element is entirely contained within the SV
					is_count+=1
					is_len+=IS_stop-IS_start
					row_idx_isescan+=1
					continue			

				if ref_start <= IS_start and IS_start <= ref_stop and ref_stop <= IS_stop: #this case is when the IS element begins inside the SV and then extends beyond it
					is_count+=1
					is_len+=ref_stop-IS_start
					row_idx_isescan+=1
					continue

				if IS_start <= ref_start and IS_stop >= ref_stop:#this is when the IS element begins before and ends after the SV (i.e SV is contained within an IS element)
					is_count+=1
					is_len+=ref_stop-ref_start
					row_idx_isescan+=1
					continue

				else:
					row_idx_isescan+=1

			difference=ref_length-is_len # length of the SV that is not an IS element
			
			if ref_length >= threshold and difference <= threshold:
				df_syri_annotated.loc[row_idx_syri,'annotation_ref']=1

		row_idx_syri+=1


def rename_deletions(df_syri): 

	''' Renaming NOTAL and copyloss SVs to DEL, as they are effectively deletion '''

	row_idx_syri=0
	for row_idx_syri in range(len(df_syri)):
		if df_syri.loc[row_idx_syri,'tag_4'] == "copyloss": #for lines that DUP or INVDUP but are noted as copyloss, effectively a deletion
			df_syri.loc[row_idx_syri, 'tag_1']="DEL"
			df_syri.loc[row_idx_syri,'tag_3']="DEL"
			df_syri.loc[row_idx_syri,'tag_4']="-"
			df_syri.loc[row_idx_syri,'query_start']="-"
			df_syri.loc[row_idx_syri,'query_stop']="-"

		#if "CPL" in df_syri.loc[row_idx_syri,'tag_1']: #for lines that DUP or INVDUP but are noted as copyloss, effectively a deletion
		#	df_syri.loc[row_idx_syri, 'tag_1']="DEL"
		#	df_syri.loc[row_idx_syri,'tag_3']="DEL"

		if "NOTAL" in df_syri.loc[row_idx_syri,'tag_1']:
			if "-" == df_syri.loc[row_idx_syri, 'query_start'] and "-" == df_syri.loc[row_idx_syri,'query_stop']:
				df_syri.loc[row_idx_syri, 'tag_1']="DEL"
				df_syri.loc[row_idx_syri,'tag_3']="DEL"

	#print("renamed tags syri")

	return df_syri

def remove_IS_matches(df_syri): 
	
	''' removing lines from the syri dataframe based on annotations provided by syri_annotate_iselements_isescan() and syri_annotate_iselements_gff3() '''

	df_syri=convert_column_int(df_syri,'ref_start')
	df_syri=convert_column_int(df_syri,'ref_stop')

	threshold=25 #a grace region allowed to when checking if two SVs are adjacent to each other
	row_idx_1=0
	rows_to_drop=[] # stores row indices of rows to be removed from the df_syri dataframe 
	remove_tags=["CPG", "INS", "INVDP", "DUP"] #tags to be removed if they correspond to an IS element
	for row_idx_1 in range(len(df_syri)): #loop over lines of df_syri

		if df_syri.loc[row_idx_1,'ref_start'] == '-' or df_syri.loc[row_idx_1,'ref_stop'] == '-': #skip any SV that does not have its reference coords listed
			continue

		if "CPL" in df_syri.loc[row_idx_1,'tag_1']:
			if df_syri.loc[row_idx_1, "annotation_ref"] == 0.0 and df_syri.loc[row_idx_1, "annotation_query"] == 1.0:

				#when a deletion occurs between two IS elements, and leaves behind one IS element

				df_syri.loc[row_idx_1, 'tag_1']="DEL"
				df_syri.loc[row_idx_1,'tag_3']="DEL"
				df_syri.loc[row_idx_1,'tag_4']="-"
				df_syri.loc[row_idx_1,'query_start']="-"
				df_syri.loc[row_idx_1,'query_stop']="-"
				df_syri.loc[row_idx_1, "annotation_query"] = 0.0

		if df_syri.loc[row_idx_1, "annotation_ref"] == 1.0 and df_syri.loc[row_idx_1, "annotation_query"] == 1.0:

			#this is either an IS element transposition event, or a random match between IS elements that got misclassified.

			ref_start_1=int(df_syri.loc[row_idx_1,"ref_start"])
			ref_stop_1=int(df_syri.loc[row_idx_1,"ref_stop"])
			tag=df_syri.loc[row_idx_1, "tag_1"]
			row_idx_2=0

			for row_idx_2 in range(len(df_syri)): #loop to find which SV the IS element is nested in. If it is nested inside another SV, we can drop it from the dataframe

				if df_syri.loc[row_idx_2,'ref_start'] == '-' or df_syri.loc[row_idx_1,'ref_stop'] == '-': #skip any SV that does not have its reference coords listed
					continue

				if row_idx_2 != row_idx_1: #prevent self matches
					ref_start_2=int(df_syri.loc[row_idx_2,"ref_start"])
					ref_stop_2=int(df_syri.loc[row_idx_2,"ref_stop"])
					# checking to see if the IS element is contained inside another SV, with a grace region of 
					if ref_start_2-ref_start_1 < threshold and ref_stop_1-ref_stop_2 < threshold:
						
						rows_to_drop.append(row_idx_1)

		if any(tag in df_syri.loc[row_idx_1, "tag_1"] for tag in remove_tags): #current line has a tag from the above list
			if df_syri.loc[row_idx_1, "annotation_ref"] == 0.0 and df_syri.loc[row_idx_1, "annotation_query"] == 1.0: 

				# this corresponds to a new IS element insertion in the query (evolved clone)
				ref_start_1=int(df_syri.loc[row_idx_1,"ref_start"])
				ref_stop_1=int(df_syri.loc[row_idx_1,"ref_stop"])
				tag=df_syri.loc[row_idx_1, "tag_1"]
				row_idx_2=0
				for row_idx_2 in range(len(df_syri)): #loop to find which SV the IS element is nested in
					if df_syri.loc[row_idx_2,'ref_start'] == '-' or df_syri.loc[row_idx_1,'ref_stop'] == '-': #skip any SV that does not have its reference coords listed
						continue
					if row_idx_2 != row_idx_1: #prevent self matches
						ref_start_2=int(df_syri.loc[row_idx_2,"ref_start"])
						ref_stop_2=int(df_syri.loc[row_idx_2,"ref_stop"])
						# checking to see if the IS element is contained inside another SV, with a grace region of 
						if ref_start_2-ref_start_1 < threshold and ref_stop_1-ref_stop_2 < threshold:
							rows_to_drop.append(row_idx_1)
	
	df_syri=df_syri.drop(rows_to_drop)
	df_syri.reset_index(drop=True, inplace=True)
	#print("Dropped cpg ins IS elements")

	print("PROGRESS: Dropped IS elements")
	return df_syri



def deleted_IS_reassign(df_syri): 
	
	''' when an IS element was part of a larger deletion '''

	row_idx_1=0
	for row_idx_1 in range(len(df_syri)): #loop over lines of df_syri
		if df_syri.loc[row_idx_1,'ref_start'] == '-' or df_syri.loc[row_idx_1,'ref_stop'] == '-': #skip any SV that does not have its reference coords listed
			continue
		if df_syri.loc[row_idx_1, "annotation_ref"] == 1 and df_syri.loc[row_idx_1, "annotation_query"] == 1: #annotated as 1 i.e this SV is an IS element
			if df_syri.loc[((row_idx_1)-1),"tag_1"] == "DEL" and df_syri.loc[((row_idx_1)+1),"tag_1"] == "DEL": #is the SV flanked by deletions
				df_syri.loc[row_idx_1, 'tag_1']="DEL"
				df_syri.loc[row_idx_1,'tag_3']="DEL"
				df_syri.loc[row_idx_1,'query_start']="-"
				df_syri.loc[row_idx_1,'query_stop']="-"

	#print("renamed deleted IS elements")
	return df_syri


def SV_IS_reassign(df_syri): 

	''' When an IS element spans the border of two adjacent SVs, which usually have an IS element '''

	rows_to_drop=[] #stores indices of rows that need to be dropped
	row_idx_1=1 #this index should go from the second line to the penultimate line
	for row_idx_1 in range(len(df_syri)-1): #loop over lines of df_syri
		if df_syri.loc[row_idx_1,'ref_start'] == '-' or df_syri.loc[row_idx_1,'ref_stop'] == '-': #skip any SV that does not have its reference coords listed
			continue
		if df_syri.loc[row_idx_1, "annotation_ref"] == 1.0 and df_syri.loc[row_idx_1, "annotation_query"] == 1.0: 
			stop_coord_upstream=int(df_syri.loc[(row_idx_1-1), 'ref_stop']) #stop coords of the SV directly upstream
			start_coord_downstream=int(df_syri.loc[(row_idx_1+1),'ref_start']) #start coords of the SV directly downstream
			if start_coord_downstream < stop_coord_upstream: #the upstream and downstream SV overlap
				rows_to_drop.append(row_idx_1)
	df_syri=df_syri.drop(rows_to_drop)
	df_syri.reset_index(drop=True, inplace=True)

	return df_syri


def rename_cpl(df_syri): 

	''' Renaming 'CPL' as they are effectively deletion '''

	for row_idx_1 in range(len(df_syri)): #loop over lines of df_syri

		if "CPL" in df_syri.loc[row_idx_1,'tag_1']:
			df_syri.loc[row_idx_1, 'tag_1']="DEL"
			df_syri.loc[row_idx_1,'tag_3']="DEL"
			df_syri.loc[row_idx_1,'query_start']="-"
			df_syri.loc[row_idx_1,'query_stop']="-"

	#print("renamed cpl")
	return df_syri

def syn_IS_reassign(df_syri):

	'''An IS element is deleted/ IS element transposition that is flanked by syntenic regions '''

	threshold = 25
	rows_to_drop=[]
	row_idx_1=0
	for row_idx_1 in range(len(df_syri)): #loop over lines of df_syri
		if df_syri.loc[row_idx_1,'ref_start'] == '-' or df_syri.loc[row_idx_1,'ref_stop'] == '-': #skip any SV that does not have its reference coords listed
			continue
		if df_syri.loc[row_idx_1, "annotation_ref"] == 1.0 and df_syri.loc[row_idx_1, "annotation_query"] == 1.0: 
			if "SYN" in df_syri.loc[((row_idx_1)-1),"tag_1"] and "SYN" in df_syri.loc[((row_idx_1)+1),"tag_1"]: #is the SV flanked by deletions?
				stop_coord_upstream=int(df_syri.loc[(row_idx_1-1), 'query_stop']) #stop coords of the SV directly upstream
				start_coord_downstream=int(df_syri.loc[(row_idx_1+1),'query_start']) #start coords of the SV directly downstream
				if start_coord_downstream - stop_coord_upstream < threshold:
					rows_to_drop.append(row_idx_1)
					df_syri.loc[(row_idx_1-1),"ref_stop"]=df_syri.loc[row_idx_1,"ref_stop"]
	df_syri=df_syri.drop(rows_to_drop)
	df_syri.reset_index(drop=True, inplace=True)
	#print("merged Is adjacent to syn")

	return df_syri

def merge_syn_lines(df_syri): 

	'''merge adjacent syntenic lines'''

	#making a new empty data frame of the same shape as df_syri that will store the new syri lines after merging
	df_syri_filtered = pd.DataFrame(np.nan, index=range(len(df_syri)), columns=df_syri.columns)
	row_idx_1=0
	row_idx_2=0
	while row_idx_1 < len(df_syri): #iterate over all licdnes
		if row_idx_1 == (len(df_syri)-1): #this code block is executed when the last line is encountered. It write out a tsv file
			df_syri_filtered.loc[row_idx_2,:]=df_syri.iloc[row_idx_1,:]
			row_idx_2+=1
			df_syri_filtered=df_syri_filtered.dropna(how='all')

			df_syri=convert_column_int(df_syri,'ref_start')
			df_syri=convert_column_int(df_syri,'ref_stop')
			print("PROGRESS: Merged adjacent synteni regions")

			return df_syri_filtered
		if "SYN" in df_syri.loc[row_idx_1,"tag_1"]: #checks to see if the ith row is a 'SYN' line
			# print (f"this is the line causing an issue:{row_idx_1}")
			# print (f"date type of df_syri.loc[row_idx_1,'ref_stop'] {type(df_syri.loc[row_idx_1,'ref_stop'])}")
			df_syri_filtered.loc[row_idx_2,:]=df_syri.iloc[row_idx_1,:]
			row_range=range((row_idx_1+1),(len(df_syri)))
			for row_idx_3 in row_range: #iterates from current line to the end of the data frame looking for adjacent 'SYN' lines
				if "SYN" in df_syri.loc[row_idx_3,"tag_1"]: #found an adjacent 'SYN' line. Next few lines update coords
					stop_coords_ref=df_syri.loc[row_idx_3,"ref_stop"]
					stop_coords_query=df_syri.loc[row_idx_3,"query_stop"]
					df_syri_filtered.loc[row_idx_2,"ref_stop"] = stop_coords_ref
					df_syri_filtered.loc[row_idx_2,"query_stop"] = stop_coords_query
					continue
				if not("SYN" in df_syri.loc[row_idx_3,"tag_1"]): #found a line that is not 'SYN', no more 'SYN' lines ahead left to merge, break out of the loop
					break
			row_idx_2+=1
			row_idx_1=row_idx_3 #row has the row immediately after the last of the continuous 'SYN' lines, skipping to that line
			continue 
		if not("SYN" in df_syri.loc[row_idx_1,'tag_1']): #if not a SYN line, copy the line as is, and move on
			df_syri_filtered.loc[row_idx_2,:]=df_syri.iloc[row_idx_1,:]
			row_idx_2+=1
			row_idx_1+=1
			continue

def merge_del_lines(df_syri): 

	'''merge adjacet deletions'''

	row_idx_1=0
	row_idx_2=0
	#making a new empty data frame of the same shape as df_syri that will store the new syri lines after merging
	df_syri_new = pd.DataFrame(np.nan, index=range(len(df_syri)), columns=df_syri.columns)
	while row_idx_1 < len(df_syri): #iterate over all lines
		if row_idx_1 == (len(df_syri)-1): #this code block is executed when the last line is encountered. It write out a tsv file
			df_syri_new.iloc[row_idx_2,:]=df_syri.iloc[row_idx_1,:]
			row_idx_2+=1
			row_idx_1+=1

			df_syri=convert_column_int(df_syri,'ref_start')
			df_syri=convert_column_int(df_syri,'ref_stop')

			print("PROGRESS: Merged adjacent deletions")

			df_syri_new=df_syri_new.dropna(how="all")
			return df_syri_new
		if "DEL" in df_syri.loc[row_idx_1,"tag_1"]: #checks to see if the ith row is a 'DEL' line
			df_syri_new.loc[row_idx_2,:]=df_syri.loc[row_idx_1,:]
			stop_coords_ref_1=int(df_syri.iloc[row_idx_1,2])
			for row_idx_3 in range((row_idx_1+1),(len(df_syri))): #iterates from current line to the end of the data frame looking for adjacent 'DEL' lines
				if df_syri.loc[row_idx_3,'ref_start'] == '-' or df_syri.loc[row_idx_3,'ref_stop'] == '-': #skip any SV that does not have its reference coords listed
					continue
				if "DEL" in df_syri.loc[row_idx_3,"tag_1"]: #found an adjacent 'DEL' line. Next few lines update coords
					start_coords_ref_2=int(df_syri.loc[row_idx_3,"ref_start"])
					stop_coords_ref_2=int(df_syri.loc[row_idx_3,"ref_stop"])
					difference=start_coords_ref_2-stop_coords_ref_1
					if difference <= 50:
						df_syri_new.loc[row_idx_2,"ref_stop"]=stop_coords_ref_2
						stop_coords_ref_1=stop_coords_ref_2
					if difference > 50:
						break
				if not("DEL" in df_syri.loc[row_idx_3,"tag_1"]): #found a line that is not 'DEL', no more 'DEL' lines ahead left to merge, break out of the loop
					break
			row_idx_2+=1
			row_idx_1=row_idx_3 #jth row has the row immediately after the last of the continuous 'DEL' lines, skipping to that line
			continue 
		if not("DEL" in df_syri.loc[row_idx_1,"tag_1"]): #if not a DEL line, copy the line as is, and move on
			df_syri_new.loc[row_idx_2,:]=df_syri.iloc[row_idx_1,:]
			row_idx_2+=1
			row_idx_1+=1

def deleted_IS_near_syn(df_syri): 
	
	''' when an IS element was part of a larger deletion, with a syntenic region upstream and deletion downstream '''

	row_idx_1=0
	for row_idx_1 in range(len(df_syri)): #loop over lines of df_syri
		if df_syri.loc[row_idx_1,'ref_start'] == '-' or df_syri.loc[row_idx_1,'ref_stop'] == '-': #skip any SV that does not have its reference coords listed
			continue
		if df_syri.loc[row_idx_1, "annotation_ref"] == 1.0 and df_syri.loc[row_idx_1, "annotation_query"] == 1.0: #annotated as 1 i.e this SV is an IS element
			if "SYN" in df_syri.loc[((row_idx_1)-1),"tag_1"] and df_syri.loc[((row_idx_1)+1),"tag_1"] == "DEL": #is the SV flanked by a syntenic region and deletion downstream
				df_syri.loc[row_idx_1, 'tag_1']="DEL"
				df_syri.loc[row_idx_1,'tag_3']="DEL"
				df_syri.loc[row_idx_1,'tag_4']="-"
				df_syri.loc[row_idx_1,'query_start']="-"
				df_syri.loc[row_idx_1,'query_stop']="-"

	#print("renamed deleted IS elements")
	return df_syri

def IS_in_IS(df_syri): 
	
	''' REL606 has an IS1 inside an IS150 at coords 588494-590471. If this region is annotated as an SV other than a deletion, it is likely misclassified / misidentified  '''

	IS_start_ref=588494
	IS_stop_ref=590471
	threshold= 50
	row_idx_1=0
	rows_to_drop=[]
	for row_idx_1 in range(len(df_syri)): #loop over lines of df_syri
		if df_syri.loc[row_idx_1,'ref_start'] == '-' or df_syri.loc[row_idx_1,'ref_stop'] == '-': #skip any SV that does not have its reference coords listed
			continue
		if df_syri.loc[row_idx_1, "annotation_ref"] == 0 and df_syri.loc[row_idx_1, "annotation_query"] == 0: #not an IS element, skip			
			continue
		if df_syri.loc[row_idx_1,"ref_start"]-IS_start_ref < threshold and df_syri.loc[row_idx_1,"ref_stop"]-IS_stop_ref < threshold and "copygain" in df_syri.loc[row_idx_1,"tag_4"]: #syri has tagged this IS1 in IS150 as a duplication 
			ref_start_1=int(df_syri.loc[row_idx_1,"ref_start"])
			ref_stop_1=int(df_syri.loc[row_idx_1,"ref_stop"])
			tag=df_syri.loc[row_idx_1, "tag_1"]
			row_idx_2=0
			for row_idx_2 in range(len(df_syri)): #loop to find which SV the IS element is nested in
				if df_syri.loc[row_idx_2,'ref_start'] == '-' or df_syri.loc[row_idx_1,'ref_stop'] == '-': #skip any SV that does not have its reference coords listed
					continue
				if row_idx_2 == row_idx_1: #prevent self matches
					continue
				ref_start_2=int(df_syri.loc[row_idx_2,"ref_start"])
				ref_stop_2=int(df_syri.loc[row_idx_2,"ref_stop"])
				# checking to see if the IS element is contained inside another SV, with a grace region of 
				if ref_start_2-ref_start_1 < threshold and ref_stop_1-ref_stop_2 < threshold:
					rows_to_drop.append(row_idx_1)
	
	df_syri=df_syri.drop(rows_to_drop)
	df_syri.reset_index(drop=True, inplace=True)
	return df_syri

	#print("renamed deleted IS elements")
	return df_syri

def convert_column_int(df_syri, col_name): #accepts a dataframe and a column name (string) of that dataframe in which all float values should be converted to string while leaving
	for row_idx in range(len(df_syri)):
		if df_syri.loc[row_idx,col_name] != '-':
			df_syri.loc[row_idx,col_name]=int(df_syri.loc[row_idx,col_name])
	return df_syri

def rename_missing_chrid(df_syri):

	''' Some rows are missing the chromosome ID of the query. This is necessary for plotsr '''

	chr_ID=df_syri.loc[0,"ref_ID"]
	df_syri["query_ID"]=chr_ID
	return df_syri


def add_deletion_coords(df_syri):

	''' query coords for some deletions are missing and need to be added '''

	for row_idx_1 in range(len(df_syri)):
		if "DEL" in df_syri.loc[row_idx_1,"tag_1"]:
			
			del_start_ref=int(df_syri.loc[row_idx_1,"ref_start"])
			del_stop_ref=int(df_syri.loc[row_idx_1,"ref_stop"])

			for row_idx_2 in range(len(df_syri)): #loop to find which SV the IS element is nested in
				
				if df_syri.loc[row_idx_2,'ref_start'] == '-' or df_syri.loc[row_idx_2,'ref_stop'] == '-': #skip any SV that does not have its reference coords listed
					continue
				
				if df_syri.loc[row_idx_2,'query_start'] == '-' or df_syri.loc[row_idx_2,'query_stop'] == '-': #skip any SV that does not have its query coords listed
					continue
				
				ref_start_SV=int(df_syri.loc[row_idx_2,"ref_start"])
				ref_stop_SV=int(df_syri.loc[row_idx_2,"ref_stop"])
				query_start_SV=int(df_syri.loc[row_idx_2, "query_start"])
				query_stop_SV=int(df_syri.loc[row_idx_2, "query_stop"])
				SV_name=df_syri.loc[row_idx_2,"tag_1"]

				# checking to see if the IS element is contained inside another SV
				if ref_start_SV<del_start_ref and del_stop_ref<ref_stop_SV:
					distance = 0
					del_query_coords = 0 # will hold the coords to be assigned to the deletion in the query sequences

					if "SYN" in SV_name:
						distance=del_start_ref-ref_start_SV
						del_query_coords=query_start_SV+distance
					
					if "INV" in SV_name:
						distance=ref_stop_SV-del_stop_ref
						del_query_coords=query_start_SV+distance
					
					df_syri.loc[row_idx_1,"query_start"]=del_query_coords
					df_syri.loc[row_idx_1,"query_stop"]=del_query_coords
	
	df_syninv = df_syri[df_syri['tag_1'].str.contains("SYN|INV")] # only contains the lines corresponding to inversions and syntenic regions of df_syri
	index_list = df_syninv.index.tolist()

	for row_idx_1 in range(len(df_syri)):
		if "DEL" in df_syri.loc[row_idx_1,"tag_1"]:
			
			if df_syri.loc[row_idx_1,"query_start"] != '-' and df_syri.loc[row_idx_1,"query_stop"] != '-': # if the query coords exist, skip
				continue

			del_start_ref=int(df_syri.loc[row_idx_1,"ref_start"])
			del_stop_ref=int(df_syri.loc[row_idx_1,"ref_stop"])
			query_coords_del=0

			for row_idx_2 in range(len(index_list)-1): # loop over the indices of df_syninv until the penultimate index

				index_1=index_list[row_idx_2]
				index_2=index_list[(row_idx_2)+1]

				if index_1 < row_idx_1 and index_2 > row_idx_1:
					query_coords_del=df_syninv.loc[index_1,"query_stop"]
					
			df_syri.loc[row_idx_1,"query_start"]=query_coords_del
			df_syri.loc[row_idx_1,"query_stop"]=query_coords_del
	
	print("PROGRESS: Added query coords for deletions")
	return df_syri
								


def main(syri, verbose, isescan_subject, isescan_query):

	# fetch_repeat_regions(reference,prefix)

	# df2=read_gff3(reference, prefix)

	df_syri= syri_remove_duplicate_lines(syri)

	print("removed duplicate lines: ", df_syri)

	#df_syri_annotated=syri_annotate_iselements_gff3(df_syri, df2)

	#print("syri_annotate_iselements_gff3: ", df_syri)

	df_syri = syri_annotate_subject_iselements_isescan(df_syri,isescan_subject)

	print("syri_annotate_iselements_isescan subject ", df_syri)

	df_syri = syri_annotate_query_iselements_isescan(df_syri,isescan_query)

	print("syri_annotate_iselements_isescan: ", df_syri)

	df_syri = rename_deletions(df_syri)

	print("rename_deletions", df_syri)

	df_syri=deleted_IS_reassign(df_syri)

	print("deleted IS reassign", df_syri)

	df_syri=SV_IS_reassign(df_syri)

	print("SV_IS_reassign", df_syri)

	df_syri=rename_cpl(df_syri)

	print("rename_cpl", df_syri)

	df_syri=remove_IS_matches(df_syri)

	print("remove_IS_matches", df_syri)

	df_syri=syn_IS_reassign(df_syri)

	print("syn_IS_reassign", df_syri)

	df_syri=deleted_IS_near_syn(df_syri)

	df_syri = IS_in_IS(df_syri)

	df_syri=merge_del_lines(df_syri)

	print("merge_del_linesn", df_syri)

	df_syri=merge_syn_lines(df_syri)

	print("erge_syn_lines", df_syri)

	df_syri=add_deletion_coords(df_syri)

	print("add_deletion_coords", df_syri)

	df_syri=rename_missing_chrid(df_syri)

	print("rename_missing_chrid", df_syri)

	output_filename=syri+"_v2"
	columns_to_write=["ref_ID", "ref_start", "ref_stop", "seq_deleted", "seq_inserted", "query_ID", "query_start", "query_stop", "tag_1", "tag_2", "tag_3", "tag_4"]
	df_syri=df_syri[columns_to_write]

	df_syri=convert_column_int(df_syri, "ref_start")
	df_syri=convert_column_int(df_syri, "ref_stop")


	df_syri.to_csv(output_filename, sep='\t', index=False, header=False,  float_format='%.0f')


if __name__ == '__main__':
    args = parser.parse_args()
    main(args.syri, args.verbose, args.isescan_subject, args.isescan_query)