#!/usr/bin/env python3
"""
Unit tests for reindex_assembly.py
"""

__author__ = "Ira Zibbu"
__version__ = "0.0.1"

"""
imports
"""

#!/usr/bin/env python3
"""
Unit tests for reindex_assembly.py
"""

__author__ = "Ira Zibbu"
__version__ = "0.1.0"


from Bio import SeqIO
from src.seabreeze.reindex_assembly import reindex_circular_sequence

ltee_test_data = "test/data/Ara-5_75000gen_A_ori.1.fasta"
reindex_bases = "AGCTTTTCATTC"


def test_reindexing_ltee_genome():
    """
    Verify that a circular genome is correctly rotated
    so that it starts at the requested conserved sequence.
    """
    record = SeqIO.read(ltee_test_data, "fasta")

    reindexed_record = reindex_circular_sequence(
        record,
        reindex_bases,
    )

    # Sequence now starts at the requested motif
    assert str(reindexed_record.seq[: len(reindex_bases)]) == reindex_bases

    # Sequence length should be preserved after rotation
    assert len(reindexed_record.seq) == len(record.seq)

    # Rotation should preserve total sequence content
    assert sorted(str(reindexed_record.seq)) == sorted(str(record.seq))
