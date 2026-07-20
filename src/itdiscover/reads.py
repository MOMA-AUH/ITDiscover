"""Sequencing read preprocessing utilities."""

from collections.abc import Iterable
from dataclasses import dataclass

from .insertions import Direction
from .sequences import reverse_complement, validate_sequence


@dataclass(frozen=True)
class SequencingRead:
    """A sequencing read oriented to the forward WT reference."""

    read_id: str
    fragment_id: str
    sequence: str
    qualities: tuple[int, ...]
    direction: Direction

    def __post_init__(self) -> None:
        if len(self.sequence) != len(self.qualities):
            raise ValueError("sequence and qualities must have equal length")
        if not self.fragment_id:
            raise ValueError("fragment_id is required")
        validate_sequence(self.sequence)

    @property
    def length(self) -> int:
        """Return the read sequence length."""
        return len(self.sequence)

    @property
    def mean_quality(self) -> float:
        """Return the mean Phred quality score."""
        if not self.qualities:
            return 0.0
        return sum(self.qualities) / len(self.qualities)


@dataclass(frozen=True)
class Fragment:
    """A paired-end sequencing fragment."""

    fragment_id: str
    forward_read: SequencingRead
    reverse_read: SequencingRead

    def __post_init__(self) -> None:
        if not self.fragment_id:
            raise ValueError("fragment_id is required")
        if self.forward_read.fragment_id != self.fragment_id:
            raise ValueError("forward read fragment_id does not match fragment")
        if self.reverse_read.fragment_id != self.fragment_id:
            raise ValueError("reverse read fragment_id does not match fragment")
        if self.forward_read.direction != "forward":
            raise ValueError("forward_read must have forward direction")
        if self.reverse_read.direction != "reverse":
            raise ValueError("reverse_read must have reverse direction")

    @property
    def reads(self) -> tuple[SequencingRead, SequencingRead]:
        """Return the fragment's oriented reads."""
        return (self.forward_read, self.reverse_read)


@dataclass(frozen=True)
class ReadTrimSettings:
    """Optional primer sequences in their raw sequencing-read orientation.

    ``forward_primer`` is removed from the 5' end of R1.  ``reverse_primer``
    is supplied as it appears at the 5' end of raw R2; after R2 has been
    reverse-complemented into reference orientation, its reverse complement is
    removed from the 3' end.
    """

    forward_primer: str | None = None
    reverse_primer: str | None = None

    def __post_init__(self) -> None:
        for field_name, sequence in (
            ("forward_primer", self.forward_primer),
            ("reverse_primer", self.reverse_primer),
        ):
            if sequence is None:
                continue
            if not sequence:
                raise ValueError(f"{field_name} must not be empty")
            validate_sequence(sequence, field_name=field_name)


@dataclass(frozen=True)
class PreprocessingMetrics:
    """Auditable read and fragment counts retained through preprocessing."""

    input_fragment_count: int
    input_read_count: int
    primer_retained_forward_reads: int | None
    primer_retained_reverse_reads: int | None
    primer_failed_read_count: int
    length_failed_read_count: int
    quality_failed_read_count: int
    passing_read_count: int
    passing_forward_read_count: int
    passing_reverse_read_count: int
    usable_fragment_count: int


@dataclass(frozen=True)
class PreprocessingResult:
    """Passing oriented reads together with preprocessing metrics."""

    reads: tuple[SequencingRead, ...]
    metrics: PreprocessingMetrics


def orient_read(
    *,
    read_id: str,
    fragment_id: str,
    sequence: str,
    qualities: Iterable[int],
    direction: Direction,
) -> SequencingRead:
    """Return a read oriented to the forward WT reference.

    Reverse-direction reads are reverse-complemented, and their qualities are
    reversed to keep base qualities aligned with the oriented sequence.
    """
    qualities = tuple(qualities)
    if direction == "forward":
        return SequencingRead(
            read_id=read_id,
            fragment_id=fragment_id,
            sequence=sequence,
            qualities=qualities,
            direction=direction,
        )
    return SequencingRead(
        read_id=read_id,
        fragment_id=fragment_id,
        sequence=reverse_complement(sequence),
        qualities=tuple(reversed(qualities)),
        direction=direction,
    )


def trim_terminal_ns(read: SequencingRead) -> SequencingRead:
    """Trim ambiguous N bases from both ends of a read."""
    trimmed, _ = _trim_terminal_ns_with_stats(read)
    return trimmed


def trim_primers(
    read: SequencingRead,
    trimming: ReadTrimSettings | None,
) -> SequencingRead | None:
    """Trim configured raw-read primer sequences from an oriented read."""
    trimmed, _ = _trim_terminal_ns_with_stats(read)
    if trimming is None:
        return trimmed

    if trimmed.direction == "forward":
        trimmed = _trim_prefix_step(trimmed, trimming.forward_primer)
    else:
        reverse_primer = (
            reverse_complement(trimming.reverse_primer)
            if trimming.reverse_primer is not None
            else None
        )
        trimmed = _trim_suffix_step(trimmed, reverse_primer)
    if trimmed is None:
        return None
    return _trim_terminal_ns_with_stats(trimmed)[0]


def passes_read_filters(
    read: SequencingRead,
    *,
    min_length: int = 100,
    min_mean_quality: float = 30,
) -> bool:
    """Return whether a read passes length and mean-quality filters."""
    if min_length < 0:
        raise ValueError("min_length must not be negative")
    if min_mean_quality < 0:
        raise ValueError("min_mean_quality must not be negative")
    return read.length >= min_length and read.mean_quality >= min_mean_quality


def preprocess_reads(
    reads: Iterable[SequencingRead],
    *,
    min_length: int = 100,
    min_mean_quality: float = 30,
    trimming: ReadTrimSettings | None = None,
) -> list[SequencingRead]:
    """Trim and filter oriented reads."""
    passing_reads: list[SequencingRead] = []
    for read in reads:
        trimmed = trim_primers(read, trimming)
        if trimmed is None:
            continue
        if trimmed.length < min_length or trimmed.mean_quality < min_mean_quality:
            continue
        passing_reads.append(trimmed)
    return passing_reads


def preprocess_fragments(
    fragments: Iterable[Fragment],
    *,
    min_length: int = 100,
    min_mean_quality: float = 30,
    trimming: ReadTrimSettings | None = None,
) -> list[SequencingRead]:
    """Return passing oriented reads from paired-end fragments."""
    reads = [read for fragment in fragments for read in fragment.reads]
    return preprocess_reads(
        reads,
        min_length=min_length,
        min_mean_quality=min_mean_quality,
        trimming=trimming,
    )


def preprocess_fragments_with_metrics(
    fragments: Iterable[Fragment],
    *,
    min_length: int = 100,
    min_mean_quality: float = 30,
    trimming: ReadTrimSettings | None = None,
) -> PreprocessingResult:
    """Preprocess paired fragments while retaining rejection and direction counts."""
    if min_length < 0:
        raise ValueError("min_length must not be negative")
    if min_mean_quality < 0:
        raise ValueError("min_mean_quality must not be negative")

    fragments = list(fragments)
    passing_reads: list[SequencingRead] = []
    primer_failed = 0
    length_failed = 0
    quality_failed = 0
    primer_retained = {"forward": 0, "reverse": 0}

    for fragment in fragments:
        for read in fragment.reads:
            trimmed = trim_primers(read, trimming)
            primer_is_configured = trimming is not None and (
                trimming.forward_primer is not None
                if read.direction == "forward"
                else trimming.reverse_primer is not None
            )
            if trimmed is None:
                primer_failed += 1
                continue
            if primer_is_configured:
                primer_retained[read.direction] += 1
            if trimmed.length < min_length:
                length_failed += 1
                continue
            if trimmed.mean_quality < min_mean_quality:
                quality_failed += 1
                continue
            passing_reads.append(trimmed)

    passing_fragment_ids = {read.fragment_id for read in passing_reads}
    forward_primer_configured = (
        trimming is not None and trimming.forward_primer is not None
    )
    reverse_primer_configured = (
        trimming is not None and trimming.reverse_primer is not None
    )
    metrics = PreprocessingMetrics(
        input_fragment_count=len(fragments),
        input_read_count=2 * len(fragments),
        primer_retained_forward_reads=(
            primer_retained["forward"] if forward_primer_configured else None
        ),
        primer_retained_reverse_reads=(
            primer_retained["reverse"] if reverse_primer_configured else None
        ),
        primer_failed_read_count=primer_failed,
        length_failed_read_count=length_failed,
        quality_failed_read_count=quality_failed,
        passing_read_count=len(passing_reads),
        passing_forward_read_count=sum(
            read.direction == "forward" for read in passing_reads
        ),
        passing_reverse_read_count=sum(
            read.direction == "reverse" for read in passing_reads
        ),
        usable_fragment_count=len(passing_fragment_ids),
    )
    return PreprocessingResult(reads=tuple(passing_reads), metrics=metrics)


def _trim_terminal_ns_with_stats(
    read: SequencingRead,
) -> tuple[SequencingRead, int]:
    start = 0
    end = len(read.sequence)

    while start < end and read.sequence[start] == "N":
        start += 1
    while end > start and read.sequence[end - 1] == "N":
        end -= 1

    trimmed = SequencingRead(
        read_id=read.read_id,
        fragment_id=read.fragment_id,
        sequence=read.sequence[start:end],
        qualities=read.qualities[start:end],
        direction=read.direction,
    )
    return trimmed, len(read.sequence) - len(trimmed.sequence)


def _trim_prefix_step(
    read: SequencingRead,
    motif: str | None,
) -> SequencingRead | None:
    if motif is None:
        return read
    index = read.sequence.find(motif)
    if index < 0:
        return None
    return SequencingRead(
        read_id=read.read_id,
        fragment_id=read.fragment_id,
        sequence=read.sequence[index + len(motif) :],
        qualities=read.qualities[index + len(motif) :],
        direction=read.direction,
    )


def _trim_suffix_step(
    read: SequencingRead,
    motif: str | None,
) -> SequencingRead | None:
    if motif is None:
        return read
    index = read.sequence.rfind(motif)
    if index < 0:
        return None
    return SequencingRead(
        read_id=read.read_id,
        fragment_id=read.fragment_id,
        sequence=read.sequence[:index],
        qualities=read.qualities[:index],
        direction=read.direction,
    )
