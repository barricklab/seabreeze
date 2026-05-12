#!/usr/bin/env python3
"""
Unit tests for auto_rotate_assembly.py
"""

__author__ = "Ira Zibbu"
__version__ = "0.0.1"


"""
imports
"""

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from src.seabreeze.auto_rotate_assembly import (
    find_anchor,
    rotate_query,
)


def test_auto_rotate_query():
    """
    Verify that a circular query genome is rotated
    to match the subject coordinate system.
    """

    subject = SeqRecord(
        Seq("AAAACCCCGGGGTTTT"),
        id="subject",
    )

    query = SeqRecord(
        Seq("GGGGTTTTAAAACCCC"),
        id="query",
    )

    anchor, subject_position, query_position = find_anchor(
        subject,
        query,
        kmer_size=4,
    )

    rotated_query = rotate_query(
        query,
        subject_position,
        query_position,
    )

    assert str(rotated_query.seq) == "AAAACCCCGGGGTTTT"
