#!/usr/bin/env python3
# An earlier version of this script would only be able to find an anchor that began at index 0 in the subject
# This did not quite work out for genomes that start at repeats

__author__ = "Ira Zibbu"
__version__ = "0.1.0"

import argparse
import os

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord




def load_fasta(path):
    """
    Load a single-record FASTA file.
    """
    print(f"Loading FASTA: {path}")

    if not os.path.exists(path):
        raise IOError(f"FASTA file does not exist: {path}")

    try:
        record = SeqIO.read(path, "fasta")
    except Exception as e:
        print(f"Error parsing FASTA: {path}")
        raise e

    print(
        f"Loaded record {record.id} "
        f"({len(record.seq):,} bp)"
    )

    return record


def find_all_circular_matches(sequence, query):
    """
    Find all exact matches of query in a circular sequence.

    Returns:
        list[int]: start positions of matches
    """
    seq = str(sequence).upper()
    query = str(query).upper()

    query_length = len(query)
    seq_length = len(seq)

    circular_seq = seq + seq[: query_length - 1]

    matches = []

    # Only search within original coordinate space
    for i in range(seq_length):
        window = circular_seq[i : i + query_length]

        if window == query:
            matches.append(i)

    return matches


def find_anchor(subject, query, kmer_size=100):
    """
    Slide across the subject looking for a 100 bp anchor
    that occurs exactly once in the query.

    Returns:
        tuple(anchor_sequence, subject_position, query_position)
    """
    KMER_SIZE=kmer_size

    print("\nSearching for anchor...")

    subject_seq = subject.seq
    subject_length = len(subject_seq)

    if subject_length < KMER_SIZE:
        raise ValueError(
            f"Subject is shorter than k-mer size ({KMER_SIZE})"
        )

    subject_circular = subject_seq + subject_seq[: KMER_SIZE - 1]

    for subject_position in range(subject_length):

        anchor = subject_circular[
            subject_position : subject_position + KMER_SIZE
        ]

        if subject_position % 10000 == 0:
            print(
                f"Testing subject position "
                f"{subject_position:,}/{subject_length:,}"
            )

        query_matches = find_all_circular_matches(
            query.seq,
            anchor,
        )

        if len(query_matches) == 1:

            query_position = query_matches[0]

            print("\nUnique anchor found!")
            print(
                f"Subject coordinate: {subject_position}"
            )
            print(
                f"Query coordinate: {query_position}"
            )

            return (
                anchor,
                subject_position,
                query_position,
            )

    raise RuntimeError(
        f"No unique {KMER_SIZE}  bp anchor found."
    )


def rotate_query(
    query,
    subject_position,
    query_position,
):
    """
    Rotate query so that the anchor aligns to
    the same coordinate as in the subject.
    """
    print("\nRotating query...")

    seq = str(query.seq).upper()
    seq_length = len(seq)

    # If anchor is at query_position but should
    # start at subject_position, rotate accordingly.
    rotation_start = (
        query_position - subject_position
    ) % seq_length

    print(
        f"Query anchor at {query_position}"
    )
    print(
        f"Should align to subject coordinate "
        f"{subject_position}"
    )
    print(
        f"Rotation start = {rotation_start}"
    )

    rotated_seq = (
        seq[rotation_start:]
        + seq[:rotation_start]
    )

    rotated_record = SeqRecord(
        Seq(rotated_seq),
        id=query.id,
        description=query.description,
    )

    print("Rotation complete.")

    return rotated_record


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Find a unique 100 bp anchor in a subject "
            "genome and rotate a circular query genome "
            "to the same coordinate system."
        )
    )

    parser.add_argument(
        "--subject",
        help="Path to subject FASTA file",
        required=True,
    )

    parser.add_argument(
        "--query",
        help="Path to query FASTA file",
        required=True,
    )

    parser.add_argument(
        "--output",
        help="Output FASTA for rotated query",
        required=True,
    )

    args = parser.parse_args()

    subject = load_fasta(
        args.subject
    )

    query = load_fasta(
        args.query
    )

    anchor, subject_position, query_position = find_anchor(
        subject,
        query,
    )

    rotated_query = rotate_query(
        query,
        subject_position,
        query_position,
    )

    print(
        f"\nWriting rotated query to: "
        f"{args.output}"
    )

    SeqIO.write(
        rotated_query,
        args.output,
        "fasta",
    )

    print("Done.")


if __name__ == "__main__":
    main()
