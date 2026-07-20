"""Read-to-reference alignment."""

from dataclasses import dataclass

from Bio.Align import PairwiseAligner

from .insertions import Alignment
from .reads import SequencingRead
from .sequences import validate_sequence


@dataclass(frozen=True)
class AlignmentScoring:
    """Pairwise alignment scoring parameters."""

    match: float = 5
    mismatch: float = -10
    gap_open: float = -25
    gap_extend: float = -1
    end_gap: float = 0


@dataclass(frozen=True)
class AlignmentEvidenceFilter:
    """Minimum evidence required for an alignment to enter calling."""

    min_identity: float = 0.90
    min_on_target_fraction: float = 0.80
    min_score: float | None = None
    reject_ambiguous: bool = False

    def __post_init__(self) -> None:
        if not 0 <= self.min_identity <= 1:
            raise ValueError("min_identity must be between 0 and 1")
        if not 0 <= self.min_on_target_fraction <= 1:
            raise ValueError("min_on_target_fraction must be between 0 and 1")


def align_read_to_reference(
    read: SequencingRead,
    reference: str,
    *,
    scoring: AlignmentScoring = AlignmentScoring(),
) -> Alignment:
    """Globally align an oriented read to the WT reference."""
    validate_sequence(reference, field_name="reference")
    aligner = _build_aligner(scoring)
    pairwise_alignments = aligner.align(reference, read.sequence)
    pairwise_alignment = pairwise_alignments[0]
    aligned_reference, aligned_read = _aligned_strings_from_coordinates(
        reference,
        read.sequence,
        pairwise_alignment.coordinates,
    )
    aligned_qualities = _aligned_qualities_from_coordinates(
        read.qualities,
        pairwise_alignment.coordinates,
    )
    return Alignment(
        read_id=read.read_id,
        fragment_id=read.fragment_id,
        read_sequence=read.sequence,
        aligned_read=aligned_read,
        aligned_reference=aligned_reference,
        direction=read.direction,
        aligned_qualities=aligned_qualities,
        score=float(pairwise_alignment.score),
        is_ambiguous=len(pairwise_alignments) > 1,
    )


def passes_alignment_evidence_filters(
    alignment: Alignment,
    filters: AlignmentEvidenceFilter,
) -> bool:
    """Return whether an alignment is sufficiently on-target and unambiguous."""
    paired = [
        (read_base, ref_base)
        for read_base, ref_base in zip(
            alignment.aligned_read,
            alignment.aligned_reference,
            strict=True,
        )
        if read_base != "-" and ref_base != "-"
    ]
    if not paired:
        return False
    identity = sum(read == ref for read, ref in paired) / len(paired)
    denominator = min(
        len(alignment.read_sequence),
        sum(base != "-" for base in alignment.aligned_reference),
    )
    on_target_fraction = len(paired) / denominator if denominator else 0.0
    if identity < filters.min_identity:
        return False
    if on_target_fraction < filters.min_on_target_fraction:
        return False
    if filters.min_score is not None and (
        alignment.score is None or alignment.score < filters.min_score
    ):
        return False
    if filters.reject_ambiguous and alignment.is_ambiguous:
        return False
    return True


def _build_aligner(scoring: AlignmentScoring) -> PairwiseAligner:
    aligner = PairwiseAligner(mode="global")
    aligner.match_score = scoring.match
    aligner.mismatch_score = scoring.mismatch
    aligner.open_gap_score = scoring.gap_open
    aligner.extend_gap_score = scoring.gap_extend
    aligner.end_gap_score = scoring.end_gap
    return aligner


def _aligned_strings_from_coordinates(
    reference: str,
    read_sequence: str,
    coordinates,
) -> tuple[str, str]:
    # PairwiseAligner stores alignments as coordinate breakpoints rather than
    # gapped strings. For example, aligning reference "AAAGGGTTT" to read
    # "AAACCCGGGTTT" may produce:
    #
    #   target            0 AAA---GGGTTT  9
    #                     0 |||---|||||| 12
    #   query             0 AAACCCGGGTTT 12
    #
    # with coordinates:
    #
    #   [[0, 3, 3, 9],
    #    [0, 3, 6, 12]]
    #
    # Each adjacent coordinate pair is one step through the alignment. A step
    # that advances only the read is an insertion relative to the reference, so
    # this function emits gaps in the aligned reference for that step.
    aligned_reference_parts: list[str] = []
    aligned_read_parts: list[str] = []

    for index in range(coordinates.shape[1] - 1):
        ref_start, ref_end = coordinates[0, index], coordinates[0, index + 1]
        read_start, read_end = coordinates[1, index], coordinates[1, index + 1]
        ref_span = ref_end - ref_start
        read_span = read_end - read_start

        if ref_span == read_span:
            aligned_reference_parts.append(reference[ref_start:ref_end])
            aligned_read_parts.append(read_sequence[read_start:read_end])
        elif ref_span == 0:
            aligned_reference_parts.append("-" * read_span)
            aligned_read_parts.append(read_sequence[read_start:read_end])
        elif read_span == 0:
            aligned_reference_parts.append(reference[ref_start:ref_end])
            aligned_read_parts.append("-" * ref_span)
        else:
            raise ValueError("alignment contains an unsupported coordinate step")

    return "".join(aligned_reference_parts), "".join(aligned_read_parts)


def _aligned_qualities_from_coordinates(
    qualities: tuple[int, ...],
    coordinates,
) -> tuple[int | None, ...]:
    aligned: list[int | None] = []
    for index in range(coordinates.shape[1] - 1):
        ref_start, ref_end = coordinates[0, index], coordinates[0, index + 1]
        read_start, read_end = coordinates[1, index], coordinates[1, index + 1]
        ref_span = ref_end - ref_start
        read_span = read_end - read_start
        if ref_span == read_span or ref_span == 0:
            aligned.extend(qualities[read_start:read_end])
        elif read_span == 0:
            aligned.extend([None] * ref_span)
        else:
            raise ValueError("alignment contains an unsupported coordinate step")
    return tuple(aligned)
