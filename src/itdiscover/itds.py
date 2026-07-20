"""FLT3 internal tandem duplication classification."""

from dataclasses import dataclass
from typing import Literal

from .insertions import Insertion
from .sequences import validate_sequence

TandemOrientation = Literal["upstream", "downstream"]


@dataclass(frozen=True)
class ITD:
    """An insertion classified as an internal tandem duplication."""

    insertion: Insertion
    tandem_start: int
    tandem_sequence: str
    orientation: TandemOrientation
    spacer_prefix: str = ""
    spacer_suffix: str = ""

    def __post_init__(self) -> None:
        expected_orientation = _tandem_orientation(
            insertion_start=self.insertion.start,
            tandem_start=self.tandem_start,
            tandem_length=len(self.tandem_sequence),
        )
        if self.orientation != expected_orientation:
            raise ValueError(
                "orientation is inconsistent with the tandem interval and "
                "insertion site"
            )

    @property
    def tandem_end(self) -> int:
        """Return the inclusive end coordinate of the duplicated WT segment."""
        return self.tandem_start + len(self.tandem_sequence) - 1

    @property
    def length(self) -> int:
        """Return the duplicated sequence length."""
        return len(self.tandem_sequence)

    @property
    def spacer_sequence(self) -> str:
        """Return the combined spacer sequence flanking the copied tract."""
        return f"{self.spacer_prefix}{self.spacer_suffix}"

    @property
    def spacer_length(self) -> int:
        """Return the total spacer length."""
        return len(self.spacer_prefix) + len(self.spacer_suffix)

    @property
    def is_partial_observation(self) -> bool:
        """Return whether the ITD reaches a read edge and may be incomplete.

        A read-edge insertion lacks sequence on one side of the event.  Its
        observed copied tract is useful evidence, but it cannot establish the
        full ITD length or sequence without other reconstruction evidence.
        """
        return self.insertion.trailing


@dataclass(frozen=True)
class TandemSimilarity:
    """Similarity between an inserted sequence and a WT tandem window."""

    insertion: Insertion
    tandem_start: int
    tandem_sequence: str
    mismatches: int

    @property
    def matches(self) -> int:
        """Return the number of matching bases in the best window."""
        return len(self.tandem_sequence) - self.mismatches

    @property
    def identity(self) -> float:
        """Return the fraction of matching bases in the best window."""
        if not self.tandem_sequence:
            return 0.0
        return self.matches / len(self.tandem_sequence)


@dataclass(frozen=True)
class _TandemMatch:
    insertion_start: int
    insertion_end: int
    tandem_start: int
    tandem_sequence: str
    mismatches: int


def classify_exact_itd(
    insertion: Insertion,
    reference: str,
    *,
    min_tandem_length: int = 6,
) -> ITD | None:
    """Classify an insertion as an adjacent exact tandem duplication.

    The copied reference tract must be immediately upstream or downstream of
    the insertion site. Extra inserted bases may flank that copied tract and
    are represented as spacer sequence.
    """
    return classify_fuzzy_itd(
        insertion,
        reference,
        max_mismatches=0,
        min_tandem_length=min_tandem_length,
    )


def score_tandem_similarity(
    insertion: Insertion,
    reference: str,
) -> TandemSimilarity | None:
    """Score how well an inserted sequence matches any WT tandem window."""
    _validate_reference(reference)
    match = _best_fuzzy_tandem_match(insertion.sequence, reference)
    if match is None:
        return None
    return TandemSimilarity(
        insertion=insertion,
        tandem_start=match.tandem_start,
        tandem_sequence=match.tandem_sequence,
        mismatches=match.mismatches,
    )


def classify_fuzzy_itd(
    insertion: Insertion,
    reference: str,
    *,
    max_mismatches: int,
    min_tandem_length: int = 6,
) -> ITD | None:
    """Classify an insertion as a fuzzy-match tandem duplication with spacers."""
    if max_mismatches < 0:
        raise ValueError("max_mismatches must not be negative")
    if min_tandem_length < 1:
        raise ValueError("min_tandem_length must be at least 1")

    _validate_reference(reference)
    match = _best_adjacent_copied_match(
        insertion,
        reference,
        max_mismatches=max_mismatches,
        min_copied_length=min_tandem_length,
    )
    if match is None:
        return None
    return _itd_from_match(insertion, match)


def _validate_reference(reference: str) -> None:
    validate_sequence(reference, field_name="reference")


def _best_adjacent_copied_match(
    insertion: Insertion,
    reference: str,
    *,
    max_mismatches: int,
    min_copied_length: int,
) -> _TandemMatch | None:
    sequence = insertion.sequence
    if not sequence:
        return None

    best_match: _TandemMatch | None = None
    best_key: tuple[int, int, int, int] | None = None

    for insertion_start in range(len(sequence)):
        for insertion_end in range(insertion_start + 1, len(sequence) + 1):
            copied_length = insertion_end - insertion_start
            if copied_length < min_copied_length or copied_length > len(reference):
                continue
            observed = sequence[insertion_start:insertion_end]
            for tandem_start in _adjacent_tandem_starts(insertion.start, copied_length):
                if tandem_start < 0 or tandem_start + copied_length > len(reference):
                    continue
                tandem_reference = reference[
                    tandem_start : tandem_start + copied_length
                ]
                mismatches = _mismatch_count(observed, tandem_reference)
                if mismatches > max_mismatches:
                    continue
                matches = copied_length - mismatches
                key = (-matches, mismatches, tandem_start, insertion_start)
                if best_key is None or key < best_key:
                    best_key = key
                    best_match = _TandemMatch(
                        insertion_start=insertion_start,
                        insertion_end=insertion_end,
                        tandem_start=tandem_start,
                        tandem_sequence=tandem_reference,
                        mismatches=mismatches,
                    )

    return best_match


def _best_fuzzy_tandem_match(
    sequence: str,
    reference: str,
) -> _TandemMatch | None:
    if not sequence:
        return None

    best_match: _TandemMatch | None = None
    best_key: tuple[int, int] | None = None

    if len(sequence) > len(reference):
        return None

    for tandem_start in range(len(reference) - len(sequence) + 1):
        tandem_sequence = reference[tandem_start : tandem_start + len(sequence)]
        mismatches = _mismatch_count(sequence, tandem_sequence)
        key = (mismatches, tandem_start)
        if best_match is None or best_key is None or key < best_key:
            best_key = key
            best_match = _TandemMatch(
                insertion_start=0,
                insertion_end=len(sequence),
                tandem_start=tandem_start,
                tandem_sequence=tandem_sequence,
                mismatches=mismatches,
            )

    return best_match


def _is_adjacent_tandem_start(insertion: Insertion, tandem_start: int) -> bool:
    sequence_length = len(insertion.sequence)
    return tandem_start in (
        insertion.start - sequence_length + 1,
        insertion.start + 1,
    )


def _adjacent_tandem_starts(insertion_start: int, copied_length: int) -> tuple[int, int]:
    return (
        insertion_start - copied_length + 1,
        insertion_start + 1,
    )


def _adjacent_tandem_candidates(
    insertion: Insertion,
    reference: str,
) -> list[tuple[TandemOrientation, str, int]]:
    sequence = insertion.sequence
    candidates: list[tuple[TandemOrientation, str, int]] = []

    upstream_start = insertion.start - len(sequence) + 1
    if upstream_start >= 0:
        upstream_sequence = reference[upstream_start : insertion.start + 1]
        candidates.append(
            (
                "upstream",
                upstream_sequence,
                _mismatch_count(sequence, upstream_sequence),
            )
        )

    downstream_start = insertion.start + 1
    downstream_end = downstream_start + len(sequence)
    if downstream_end <= len(reference):
        downstream_sequence = reference[downstream_start:downstream_end]
        candidates.append(
            (
                "downstream",
                downstream_sequence,
                _mismatch_count(sequence, downstream_sequence),
            )
        )

    return candidates


def _itd_from_match(insertion: Insertion, match: _TandemMatch) -> ITD:
    return ITD(
        insertion=insertion,
        tandem_start=match.tandem_start,
        tandem_sequence=match.tandem_sequence,
        orientation=_tandem_orientation(
            insertion_start=insertion.start,
            tandem_start=match.tandem_start,
            tandem_length=len(match.tandem_sequence),
        ),
        spacer_prefix=insertion.sequence[: match.insertion_start],
        spacer_suffix=insertion.sequence[match.insertion_end :],
    )


def _tandem_orientation(
    *,
    insertion_start: int,
    tandem_start: int,
    tandem_length: int,
) -> TandemOrientation:
    tandem_end = tandem_start + tandem_length - 1
    if tandem_end == insertion_start:
        return "upstream"
    if tandem_start == insertion_start + 1:
        return "downstream"
    raise ValueError("tandem interval is not adjacent to the insertion site")


def _mismatch_count(observed: str, expected: str) -> int:
    return sum(
        1 for observed_base, expected_base in zip(observed, expected, strict=True)
        if observed_base != expected_base
    )
