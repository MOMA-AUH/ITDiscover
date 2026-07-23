"""Canonical sequence-allele representations used for call identity."""

from dataclasses import dataclass

from .insertions import Insertion
from .sequences import validate_sequence


@dataclass(frozen=True, order=True)
class CanonicalInsertionAllele:
    """A deterministic representation of an observed insertion allele.

    ``start`` follows :class:`~itdiscover.insertions.Insertion`: it is the
    zero-based reference base immediately before the insertion, with ``-1``
    denoting an insertion before the first reference base.

    Complete insertions are shifted to the left-most equivalent position.
    During a shift the inserted sequence is rotated so applying the normalized
    allele to the reference reconstructs exactly the same alternate sequence.
    Partial read-edge observations are not shifted because they do not
    establish a complete alternate allele.
    """

    start: int
    sequence: str
    trailing: bool = False

    def __post_init__(self) -> None:
        validate_sequence(self.sequence, field_name="inserted_sequence")

    @property
    def insertion_index(self) -> int:
        """Return the zero-based inter-base index at which sequence is inserted."""
        return self.start + 1

    def alternate_sequence(self, reference: str) -> str:
        """Apply this insertion to ``reference`` and return the ALT sequence."""
        _validate_location(self.start, reference)
        index = self.insertion_index
        return reference[:index] + self.sequence + reference[index:]

    def as_insertion(self, template: Insertion) -> Insertion:
        """Return an ``Insertion`` carrying template provenance and this allele."""
        return Insertion(
            read_id=template.read_id,
            fragment_id=template.fragment_id,
            start=self.start,
            sequence=self.sequence,
            direction=template.direction,
            trailing=self.trailing,
        )


def canonicalize_insertion_allele(
    insertion: Insertion,
    reference: str,
) -> CanonicalInsertionAllele:
    """Return the canonical allele for an observed insertion."""
    return canonicalize_insertion(
        start=insertion.start,
        sequence=insertion.sequence,
        reference=reference,
        trailing=insertion.trailing,
    )


def canonicalize_insertion(
    *,
    start: int,
    sequence: str,
    reference: str,
    trailing: bool = False,
) -> CanonicalInsertionAllele:
    """Left-normalize an insertion while preserving its exact ALT sequence."""
    validate_sequence(reference, field_name="reference")
    validate_sequence(sequence, field_name="inserted_sequence")
    _validate_location(start, reference)

    if trailing:
        return CanonicalInsertionAllele(
            start=start,
            sequence=sequence,
            trailing=True,
        )

    insertion_index = start + 1
    normalized_sequence = sequence
    while (
        insertion_index > 0
        and normalized_sequence
        and normalized_sequence[-1] == reference[insertion_index - 1]
    ):
        normalized_sequence = (
            normalized_sequence[-1] + normalized_sequence[:-1]
        )
        insertion_index -= 1

    return CanonicalInsertionAllele(
        start=insertion_index - 1,
        sequence=normalized_sequence,
        trailing=False,
    )


def _validate_location(start: int, reference: str) -> None:
    if start < -1 or start >= len(reference):
        raise ValueError("insertion start is outside the reference")
