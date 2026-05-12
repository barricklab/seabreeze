#!/usr/bin/env python3

__author__ = "Ira Zibbu and ChatGPT"
__version__ = "0.1.0"

import argparse

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def verify_unique_reindex_bases(sequence, reindex_bases, record_id):
    """
    Verify that reindex_bases occurs exactly once in the sequence.

    Raises:
        ValueError: If the sequence is absent or occurs multiple times.
    """
    seq_str = str(sequence).upper()
    query = reindex_bases.upper()

    count = seq_str.count(query)

    if count == 0:
        raise ValueError(
            f"Sequence '{record_id}' does not contain reindex bases: {query}"
        )

    if count > 1:
        raise ValueError(
            f"Sequence '{record_id}' contains reindex bases multiple times: {query}"
        )

    return True


def reindex_circular_sequence(record, reindex_bases, input):
    """
    Rotate a circular sequence so it starts at reindex_bases.

    Returns:
        SeqRecord: Rotated sequence record.
    """
    seq_str = str(record.seq).upper()
    query = reindex_bases.upper()

    # Verify uniqueness before reindexing
    verify_unique_reindex_bases(seq_str, query, record.id)

    start_pos = seq_str.find(query)

    print(f"Sequence {query} found in {input} at location {start_pos}")

    # skip reindexing if the sequence starts at 0

    if start_pos==0:
        print(f"Sequence already starts at specified bases. Skipping. No change required")
        return record

    # Circular rotation
    rotated_seq = seq_str[start_pos:] + seq_str[:start_pos]

    return SeqRecord(
        Seq(rotated_seq),
        id=record.id,
        description=record.description,
    )


def main():
    parser = argparse.ArgumentParser(
        description="Reindex sequence file to a conserved sequence"
    )
    parser.add_argument(
        "-b",
        "--reindex_bases",
        help="Bases to reindex on",
        action="store",
        required=True,
    )
    parser.add_argument(
        "-i", "--input", help="Input fasta file", action="store", required=True
    )
    parser.add_argument(
        "-o", "--output", help="Output fasta file", action="store", required=True
    )
    parser.add_argument(
        "-t",
        "--filetype",
        help="Filetype to reindex [fasta]",
        action="store",
        default="fasta",
    )

    args = parser.parse_args()

    modified_records = []

    for record in SeqIO.parse(args.input, args.filetype):
        rotated_record = reindex_circular_sequence(
            record,
            args.reindex_bases, args.input
        )
        modified_records.append(rotated_record)

    SeqIO.write(
        modified_records,
        args.output,
        args.filetype,
    )


if __name__ == "__main__":
    main()
