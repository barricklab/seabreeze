#! /usr/bin/python3
# This script reindexes (shifts) the sequence in the input file to a requested subsequence

# ---- Imports ----
import argparse
import os
import importlib.util
from io import StringIO

# Optional dependency
if importlib.util.find_spec("rich") is not None:
    from rich import print


# ---- Functions ----

def reverse_complement(string):
    """Return the reverse complement of a string"""
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

def string_splitter(string, length):
    """Generator. Splits a string into chunks of a given length"""
    string = StringIO(string)
    while True:
        substring = string.read(length)
        if len(substring) == 0:
            break
        yield substring
        if len(substring) == length:
            continue
        break

def reindex_fasta(file_in, file_out, reindex_bases):
    current_contig = ""
    contig_label = ""
    contigs = []
    # Open file for reading
    with open(file_in, "r") as f:
        # Separate contigs
        for line in f:
            readline = line.strip()
            if readline.startswith(">"):
                if current_contig:
                    contigs.append([contig_label, current_contig])
                contig_label = readline
                current_contig = ""
            else:
                current_contig += readline
        if current_contig:
            contigs.append([contig_label, current_contig])

    # Attempt to find reindex bases in each contig
    for i, contig in enumerate(contigs):
        reindex_start = None

        try: # the forward sequence
            reindex_start = contig[1].index(reindex_bases)
            contig[1] = contig[1][reindex_start:]+contig[1][:reindex_start]
            print("Reindexed {} contig {}".format(file_in, i+1))
            continue
        except ValueError:
            pass

        try: # the reverse complement
            contig_rc = reverse_complement(contig[1])
            reindex_start = contig_rc.index(reindex_bases)
            contig[1] = contig_rc[reindex_start:]+contig_rc[:reindex_start]
            print("Reindexed {} contig {} as reverse complement".format(file_in, i+1))
            continue
        except ValueError:
            pass

        # Abandon hope if we can't find the reindex bases
        print("Could not find reindex bases in sequence {} contig {}".format(file_in, i+1))
        continue

    
    # Delete the output file if it already exists
    if os.path.exists(os.path.join(file_out)):
        os.remove(os.path.join(file_out))

    # Write to output file. 50 bases per line
    with open(os.path.join(file_out), "w") as f:
        for contig in contigs:
            f.write(contig[0] + "\n")
            for line in string_splitter(contig[1], 50):
                f.write(line + "\n")


def main():
    # Get arguments
    parser = argparse.ArgumentParser(description="Reindex sequence file to a conserved sequence")
    parser.add_argument("-b", "--reindex_bases", help="Bases to reindex on", action="store", required=True)
    parser.add_argument("-i", "--input", help="Folder containing target files", action="store", required=True)
    parser.add_argument("-o", "--output", help="Folder containing target files", action="store", required=True)
    parser.add_argument("-t", "--filetype", help="Filetype to reindex [fasta]", action="store", default="fasta")
    args = parser.parse_args()

    # Check if filetype is valid
    valid_filetypes = ["fasta"]
    if args.filetype.lower() not in valid_filetypes:
        print("Invalid filetype. Valid filetypes are: {}".format(", ".join(valid_filetypes)))
        exit(1)
    
    
    reindex_fasta(args.input, args.output, args.reindex_bases)
    exit(0)
    
    # This is old code for batching a folder below
    
    # Create output folder if it doesn't exist
    output_folder = os.path.join(args.input_folder, "reindexed")
    if not os.path.exists(output_folder):
        os.mkdir(output_folder)
    
    # Get list of files in input folder
    files = os.listdir(args.input_folder)

    # Filter to input filetype
    files = [file for file in files if file.endswith("." + args.filetype.lower())]
    print(f"Found {len(files)} target files in {args.input_folder}...")

    # Reindex each file
    for file in files:
        output_file = os.path.join(output_folder, file)

        if args.filetype.lower() == "fasta":
            reindex_fasta(os.path.join(args.input_folder, file), output_file, args.reindex_bases)
        else:
            assert args.filetype in valid_filetypes
            print("Something went wrong. Please report this bug.")
            exit(1)



if __name__ == '__main__':
    main()
