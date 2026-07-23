"""Insertion models and extraction from read-to-reference alignments."""

from dataclasses import dataclass, field
from typing import Literal

from .sequences import VALID_ALIGNMENT_CHARS, validate_sequence

Direction = Literal["forward", "reverse"]

@dataclass(frozen=True)
class Alignment:
    """A read aligned to a wild-type amplicon reference."""

    read_id: str
    fragment_id: str
    read_sequence: str
    aligned_read: str
    aligned_reference: str
    direction: Direction
    aligned_qualities: tuple[int | None, ...] = field(default=(), compare=False)
    score: float | None = field(default=None, compare=False)
    is_ambiguous: bool = field(default=False, compare=False)

    def __post_init__(self) -> None:
        if len(self.aligned_read) != len(self.aligned_reference):
            raise ValueError("aligned_read and aligned_reference must have equal length")
        if self.aligned_qualities and len(self.aligned_qualities) != len(
            self.aligned_read
        ):
            raise ValueError(
                "aligned_qualities must be empty or match the aligned sequence length"
            )
        if not self.fragment_id:
            raise ValueError("fragment_id is required")
        validate_sequence(self.read_sequence, field_name="read_sequence")
        validate_sequence(
            self.aligned_read,
            valid_chars=VALID_ALIGNMENT_CHARS,
            field_name="aligned_read",
        )
        validate_sequence(
            self.aligned_reference,
            valid_chars=VALID_ALIGNMENT_CHARS,
            field_name="aligned_reference",
        )


@dataclass(frozen=True)
class InsertionEvidenceFilter:
    """Base-quality safeguards for insertion evidence."""

    min_junction_quality: int = 30
    junction_flank_size: int = 3

    def __post_init__(self) -> None:
        if self.min_junction_quality < 0:
            raise ValueError("min_junction_quality must not be negative")
        if self.junction_flank_size < 1:
            raise ValueError("junction_flank_size must be at least 1")


@dataclass(frozen=True)
class Insertion:
    """An insertion relative to the wild-type amplicon reference."""

    read_id: str
    fragment_id: str
    start: int
    sequence: str
    direction: Direction
    trailing: bool = False

    def __post_init__(self) -> None:
        if not self.fragment_id:
            raise ValueError("fragment_id is required")
        validate_sequence(self.sequence)

    @property
    def length(self) -> int:
        """Return the inserted sequence length."""
        return len(self.sequence)


def extract_insertions(
    alignment: Alignment,
    *,
    min_length: int = 6,
    require_in_frame: bool = True,
    evidence_filter: InsertionEvidenceFilter | None = None,
) -> list[Insertion]:
    """Extract insertions from an aligned read/reference pair.

    Insertions are runs where the aligned reference contains gaps and the
    aligned read contains bases. Coordinates are zero-based and `start` is the
    reference base immediately before the insertion. A leading insertion uses
    `start == -1`.
    """
    if min_length < 1:
        raise ValueError("min_length must be at least 1")

    insertions: list[Insertion] = []
    ref_pos = -1
    i = 0

    while i < len(alignment.aligned_reference):
        read_base = alignment.aligned_read[i]
        ref_base = alignment.aligned_reference[i]

        if ref_base != "-":
            ref_pos += 1
            i += 1
            continue

        if read_base == "-":
            i += 1
            continue

        start = ref_pos
        insert_bases: list[str] = []
        insert_start_index = i
        while (
            i < len(alignment.aligned_reference)
            and alignment.aligned_reference[i] == "-"
            and alignment.aligned_read[i] != "-"
        ):
            insert_bases.append(alignment.aligned_read[i])
            i += 1

        sequence = "".join(insert_bases)
        trailing = insert_start_index == 0 or i == len(alignment.aligned_reference)
        if _passes_insertion_filters(
            sequence,
            min_length=min_length,
            trailing=trailing,
            require_in_frame=require_in_frame,
            alignment=alignment,
            insert_start_index=insert_start_index,
            insert_end_index=i,
            evidence_filter=evidence_filter,
        ):
            insertions.append(
                Insertion(
                    read_id=alignment.read_id,
                    fragment_id=alignment.fragment_id,
                    start=start,
                    sequence=sequence,
                    direction=alignment.direction,
                    trailing=trailing,
                )
            )

    return insertions


def _passes_insertion_filters(
    sequence: str,
    *,
    min_length: int,
    trailing: bool,
    require_in_frame: bool,
    alignment: Alignment,
    insert_start_index: int,
    insert_end_index: int,
    evidence_filter: InsertionEvidenceFilter | None,
) -> bool:
    if len(sequence) < min_length:
        return False
    if "N" in sequence:
        return False
    # Trailing insertions may be partial observations clipped by the read edge,
    # so only fully internal insertions are required to be in-frame here.
    if require_in_frame and not trailing and len(sequence) % 3 != 0:
        return False
    if evidence_filter is not None:
        if not _passes_junction_quality(
            alignment,
            insert_start_index,
            insert_end_index,
            evidence_filter,
        ):
            return False
    return True


def _passes_junction_quality(
    alignment: Alignment,
    insert_start_index: int,
    insert_end_index: int,
    evidence_filter: InsertionEvidenceFilter,
) -> bool:
    # Manually constructed alignments from the public low-level API may not
    # carry qualities. FASTQ-derived alignments always do.
    if not alignment.aligned_qualities:
        return True

    flank = evidence_filter.junction_flank_size
    quality_indices = list(range(insert_start_index, insert_end_index))
    quality_indices.extend(
        _nearest_read_base_indices(alignment, insert_start_index - 1, -1, flank)
    )
    quality_indices.extend(
        _nearest_read_base_indices(alignment, insert_end_index, 1, flank)
    )
    if len(quality_indices) < (insert_end_index - insert_start_index) + 2 * flank:
        return False
    qualities = [alignment.aligned_qualities[index] for index in quality_indices]
    return all(
        quality is not None and quality >= evidence_filter.min_junction_quality
        for quality in qualities
    )


def _nearest_read_base_indices(
    alignment: Alignment,
    start: int,
    step: int,
    count: int,
) -> list[int]:
    indices: list[int] = []
    index = start
    while 0 <= index < len(alignment.aligned_read) and len(indices) < count:
        if alignment.aligned_read[index] != "-":
            indices.append(index)
        index += step
    return indices
