from pathlib import Path

from itdiscover.alleles import (
    CanonicalInsertionAllele,
    canonicalize_insertion,
    canonicalize_insertion_allele,
)
from itdiscover.insertions import Insertion


def test_equivalent_insertion_placements_have_one_canonical_allele() -> None:
    reference = "AAACCCGGGTTT"
    alternate = "AAACCCGGGCCCGGGTTT"
    alleles = {
        canonicalize_insertion(
            start=insertion_index - 1,
            sequence=alternate[insertion_index : insertion_index + 6],
            reference=reference,
        )
        for insertion_index in range(len(reference) + 1)
        if (
            alternate[:insertion_index] == reference[:insertion_index]
            and alternate[insertion_index + 6 :] == reference[insertion_index:]
        )
    }

    assert alleles == {
        CanonicalInsertionAllele(start=2, sequence="CCCGGG")
    }


def test_true_flt3_repeat_placements_have_one_canonical_allele() -> None:
    reference_path = (
        Path(__file__).parent / "data" / "synthetic_flt3" / "reference.fa"
    )
    reference = "".join(reference_path.read_text().splitlines()[1:])
    copied = reference[79:94]
    alternate = reference[:79] + copied + reference[79:]
    placements = [
        canonicalize_insertion(
            start=insertion_index - 1,
            sequence=alternate[insertion_index : insertion_index + len(copied)],
            reference=reference,
        )
        for insertion_index in range(len(reference) + 1)
        if (
            alternate[:insertion_index] == reference[:insertion_index]
            and alternate[insertion_index + len(copied) :]
            == reference[insertion_index:]
        )
    ]

    assert len(placements) == 16
    assert len(set(placements)) == 1
    assert all(
        allele.alternate_sequence(reference) == alternate for allele in placements
    )


def test_canonicalization_is_idempotent() -> None:
    reference = "AAACCCGGGTTT"
    raw = Insertion(
        read_id="read",
        fragment_id="fragment",
        start=8,
        sequence="CCCGGG",
        direction="forward",
    )

    once = canonicalize_insertion_allele(raw, reference)
    twice = canonicalize_insertion_allele(once.as_insertion(raw), reference)

    assert once == twice


def test_distinct_alt_sequences_remain_distinct() -> None:
    reference = "AAACCCGGGTTTGGGCCCAAA"
    first = canonicalize_insertion(
        start=2,
        sequence="CCCGGG",
        reference=reference,
    )
    second = canonicalize_insertion(
        start=11,
        sequence="GGGCCC",
        reference=reference,
    )

    assert first != second
    assert first.alternate_sequence(reference) != second.alternate_sequence(reference)


def test_partial_observation_is_not_shifted() -> None:
    allele = canonicalize_insertion(
        start=8,
        sequence="CCCGGG",
        reference="AAACCCGGGTTT",
        trailing=True,
    )

    assert allele == CanonicalInsertionAllele(
        start=8,
        sequence="CCCGGG",
        trailing=True,
    )
