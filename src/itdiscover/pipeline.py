"""High-level in-memory ITD calling pipeline."""

from collections.abc import Iterable

from .alignment import (
    AlignmentEvidenceFilter,
    AlignmentScoring,
    align_read_to_reference,
    passes_alignment_evidence_filters,
)
from .calls import ITDCall, ITDFilter, call_exact_itds, call_fuzzy_itds
from .reads import Fragment, ReadTrimSettings, preprocess_fragments
from .insertions import InsertionEvidenceFilter
from .sequences import validate_sequence


def call_exact_itds_from_fragments(
    fragments: Iterable[Fragment],
    reference: str,
    *,
    min_read_length: int = 100,
    min_mean_quality: float = 30,
    min_insert_length: int = 6,
    min_tandem_length: int | None = None,
    require_in_frame: bool = True,
    trimming: ReadTrimSettings | None = None,
    filters: ITDFilter = ITDFilter(),
    scoring: AlignmentScoring = AlignmentScoring(),
    alignment_filters: AlignmentEvidenceFilter = AlignmentEvidenceFilter(),
    insertion_filters: InsertionEvidenceFilter = InsertionEvidenceFilter(),
) -> list[ITDCall]:
    """Call exact-match ITDs from paired-end sequencing fragments."""
    validate_sequence(reference, field_name="reference")
    processed_reads = preprocess_fragments(
        fragments,
        min_length=min_read_length,
        min_mean_quality=min_mean_quality,
        trimming=trimming,
    )
    alignments = [
        align_read_to_reference(
            read,
            reference,
            scoring=scoring,
            detect_ambiguous_events=alignment_filters.reject_ambiguous,
        )
        for read in processed_reads
    ]
    alignments = [
        alignment
        for alignment in alignments
        if passes_alignment_evidence_filters(alignment, alignment_filters)
    ]
    return call_exact_itds(
        alignments,
        reference,
        min_insert_length=min_insert_length,
        min_tandem_length=min_tandem_length,
        require_in_frame=require_in_frame,
        filters=filters,
        evidence_filter=insertion_filters,
    )


def call_fuzzy_itds_from_fragments(
    fragments: Iterable[Fragment],
    reference: str,
    *,
    max_mismatches: int,
    min_read_length: int = 100,
    min_mean_quality: float = 30,
    min_insert_length: int = 6,
    min_tandem_length: int | None = None,
    require_in_frame: bool = True,
    trimming: ReadTrimSettings | None = None,
    filters: ITDFilter = ITDFilter(),
    scoring: AlignmentScoring = AlignmentScoring(),
    alignment_filters: AlignmentEvidenceFilter = AlignmentEvidenceFilter(),
    insertion_filters: InsertionEvidenceFilter = InsertionEvidenceFilter(),
) -> list[ITDCall]:
    """Call fuzzy-match ITDs from paired-end sequencing fragments."""
    validate_sequence(reference, field_name="reference")
    processed_reads = preprocess_fragments(
        fragments,
        min_length=min_read_length,
        min_mean_quality=min_mean_quality,
        trimming=trimming,
    )
    alignments = [
        align_read_to_reference(
            read,
            reference,
            scoring=scoring,
            detect_ambiguous_events=alignment_filters.reject_ambiguous,
        )
        for read in processed_reads
    ]
    alignments = [
        alignment
        for alignment in alignments
        if passes_alignment_evidence_filters(alignment, alignment_filters)
    ]
    return call_fuzzy_itds(
        alignments,
        reference,
        max_mismatches=max_mismatches,
        min_insert_length=min_insert_length,
        min_tandem_length=min_tandem_length,
        require_in_frame=require_in_frame,
        filters=filters,
        evidence_filter=insertion_filters,
    )
