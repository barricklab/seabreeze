#!/usr/bin/env python3
# this script accepts the paths to a set of fasta files to generate a series or batch of synteny plots

''' imports '''

import argparse
import os
import subprocess

''' fetch arguments '''

parser = argparse.ArgumentParser(description='align_and_visual_timeseries.py: a script to generate syneny plots for time series fasta files with rearrangement intermediates')
parser.add_argument('--verbose', '-v', action='store_true', help='Enable verbose mode')
parser.add_argument('--plotsr', help='Path to plotsr-bin')
parser.add_argument('--batch', '-b', action='store_true', help='Batch mode. If not used, series mode is run by default')
parser.add_argument('fasta', nargs='+', help='List of fasta files, in the order of generating plots')

def check_contig_count(fasta):

    ''' accepts a list of fasta files to check that all have only one entry / one contig '''

    for count, name in enumerate(fasta):

        with open(name, "r") as file:
            contig_count=0
            for line in file:
                line = line.strip()
                if line[0] == ">":
                    contig_count+=1
            if contig_count > 1:
                print("Aborting process. ", name, " has more than one contig")

def aligner(seq1, seq2, n):

    ''' accepts file names of two fasta files and aligns them, and writes out the alignment files, and returns'''

    # generate alignment
    prefix="file"+str(n)
    nucmer_command = [
        'nucmer',
        '--maxmatch',
        '-c', '100',
        '-b', '500',
        '-l', '50',
        '-p', prefix,
        seq1, seq2]

    try:
        # Run the subprocess
        subprocess.run(nucmer_command, check=True)
        print(f"Subprocess nucmer completed successfully.")
    except subprocess.CalledProcessError as e:
        print(f"Error running subprocess nucmer: {e}")

    # filter alignments for seq ID and length

    file_delta = prefix + ".delta"
    file_filtered_delta = prefix + ".filtered.delta"
    filter_command = ['delta-filter', '-i', '95', '-l', '100', file_delta]

    try:
        # Run the subprocess and redirect stdout to a file. DONT USE the PIPE > THIS IS NOT A SHELL COMMAND
        with open(file_filtered_delta, 'w') as output_file:
            result = subprocess.run(filter_command, check=True, stdout=output_file, stderr=subprocess.PIPE, text=True)

        # Check errors
        if result.stderr:
            print(f"Subprocess delta filter completed with errors:\n{result.stderr}")
        else:
            print("Subprocess delta filter completed successfully.")
    except subprocess.CalledProcessError as e:
        print(f"Error running subprocess delta filter: {e}")

    # delta show coords, to change delta file format to coords format

    file_coords=prefix+".filtered.coords"
    show_coords_command=['show-coords','-THrd',file_filtered_delta]

    try:
        # Run the subprocess and redirect stdout to a file. DONT USE the PIPE > THIS IS NOT A SHELL COMMAND
        with open(file_coords, 'w') as output_file:
            result = subprocess.run(show_coords_command, check=True, stdout=output_file, stderr=subprocess.PIPE, text=True)

        # Check errors
        if result.stderr:
            print(f"Subprocess delta show coords completed with errors:\n{result.stderr}")
        else:
            print("Subprocess delta show coords completed successfully.")
    except subprocess.CalledProcessError as e:
        print(f"Error running subprocess delta show coords: {e}")


def call_variants(delta_coords, delta_filtered, seq1, seq2, n):

    ''' accepts specified paths and uses syri to call structural variants. Returns a syri file.'''

    prefix = "file" + str(n)
    syri_command=['syri','--nosnp','-c',delta_coords,'-d',delta_filtered,'-r',seq1,'-q',seq2,'--prefix',prefix]

    try:
        # Run the subprocess
        subprocess.run(syri_command, check=True)
        print(f"Subprocess syri completed successfully.")
    except subprocess.CalledProcessError as e:
        print(f"Error running subprocess syri: {e}")

    files_to_remove = ["file"+str(n)+"syri.log","file"+str(n)+"syri.summary","file"+str(n)+"syri.vcf"] # not needed, removing to reduce clutter
    for file in files_to_remove:
        try:
            # Remove the file
            os.remove(file)
            print(f"{file} removed successfully.")
        except OSError as e:
            print(f"Error removing {file}: {e}")


def make_genomes_file(fasta,names,count):

    ''' accepts a list of fasta file paths and a list of names of the sequences and generates the tsv file needed by plotsr '''

    file_name="file"+str(count)+".genomes.tsv"

    # Writing the header to the file
    header = "#file\tname\ttags\n"
    with open(file_name, 'w') as file:
        file.write(header)

    # Adding the lines to the file
    for n, entry in enumerate(fasta): # add one line per entry in this list of fasta files
        line = f"{entry}\t{names[n]}\tlw:1.5\n"
        with open(file_name, 'a') as file:
            file.write(line)

def makeplots(syri, genomes, count, batch,length,plotsr): #syri can be single file (for batch mode) or a list of file names (for series mode)

    ''' generate the synteny plots '''

    if batch: # generate pairwise syteny plot

        prefix="file"+str(count)
        plotsr_command=[plotsr,'-s','500','--genomes',genomes,'--sr',syri,'-H','5','-W','10','-o',f"{prefix}_plotsr.pdf"]

        try:
            # Run the subprocess
            subprocess.run(plotsr_command, check=True)
            print(f"Subprocess plotsr completed successfully.")
        except subprocess.CalledProcessError as e:
            print(f"Error running plotsr nucmer: {e}")

    if not(batch): #accept an ordered list of syri files to generate one large synteny plots in series

        prefix="file"
        plotsr_command=[plotsr,'-s','500','--genomes',genomes]

        for file in syri: # add the syri option and entries to the command
            plotsr_command.append('--sr')
            plotsr_command.append(file)

        plotsr_command+=['-H',length,'-W','10','-o',f"{prefix}_plotsr.pdf"]

        try:
            # Run the subprocess
            subprocess.run(plotsr_command, check=True)
            print(f"Subprocess plotsr completed successfully.")
        except subprocess.CalledProcessError as e:
            print(f"Error running plotsr nucmer: {e}")

def main(batch, fasta,plotsr):
    check_contig_count(fasta) # abort process if fasta files have more than one entry


    if batch: # generate pair wise genomes files if true

        for count, name in enumerate(fasta): # pick up one sequence at a time to generate clone1 vs ancestor, clone2 vs ancestor and so on

            if count == len(fasta)-1: # once the last element is reached, quit
                break

            aligner(fasta[0], fasta[count+1], count)
            prefix="file"+str(count)
            delta_coords=prefix+".filtered.coords"
            delta_filtered = prefix + ".filtered.delta"
            call_variants(delta_coords, delta_filtered, fasta[count], fasta[count+1], count)

        print("Alignment and structural varaint calling done. Starting plotsr")

        for count, name in enumerate(fasta): # pick up two sequences at  time
            if count == len(fasta)-1: # once the last element is reached, quit
                break
            list_of_files=[fasta[0], fasta[count+1]]
            names=[f"Ancestor",f"Sequence_{count+1}"] # sequences are named sequence_1, sequence_2 etc in the synteny plot
            make_genomes_file(list_of_files,names,count)
            prefix="file"+str(count)
            syri=prefix+"syri.out"
            genomes=prefix+".genomes.tsv"
            makeplots(syri,genomes,count,batch,'',plotsr)

    if not(batch):

        for count, name in enumerate(fasta): # pick up one sequence at a time to generate ancestor vs clone1, clone1 vs clone 2 and so on

            if count == len(fasta)-1: # once the last element is reached, quit
                break

            aligner(fasta[count], fasta[count+1], count)
            prefix="file"+str(count)
            delta_coords=prefix+".filtered.coords"
            delta_filtered = prefix + ".filtered.delta"
            call_variants(delta_coords, delta_filtered, fasta[count], fasta[count+1], count)

        print("Alignment and structural varaint calling done. Starting plotsr")

        names=[]
        syri=[]
        for count,file in enumerate(fasta):
            names.append(f"Sequence_{count}")
            if count == len(fasta)-1:
                break
            syri.append(f"file{count}syri.out")
        make_genomes_file(fasta,names,'') # only one file is generated
        count='' # the files for series mode only have "file" as the prefix with no following numbers
        prefix="file"+str(count)
        genomes=prefix+".genomes.tsv"
        makeplots(syri,genomes,count,batch,str(5*len(fasta)),plotsr)

if __name__ == '__main__':
    args = parser.parse_args()
    main(args.batch, args.fasta,args.plotsr)
