//! Locations: `SequenceLocation`, `SequenceOffsetLocation`, `RelativeSequenceLocation` and the
//! polymorphic `Location`.

use crate::error::CoordinateError;
use crate::model::identifier::RefgetAccession;
use crate::model::meta::{Meta, impl_entity};
use crate::model::primitives::{AnchorOrientation, IntOrRange, IriOr, SequenceString};
use crate::model::sequence::SequenceReference;

/// A VRS `SequenceLocation`: an inter-residue interval `[start, end)` on a
/// [`SequenceReference`].
///
/// Invariants enforced at construction (and on deserialization):
/// * at least one of `start` / `end` is present (a location with only one coordinate is a
///   half-open "extends left/right" location used by adjacencies and termini);
/// * no coordinate or range bound is negative;
/// * `start <= end` unless the inline sequence reference is explicitly `circular: true`.
#[derive(Clone, Debug, PartialEq)]
pub struct SequenceLocation {
    sequence_reference: Option<IriOr<SequenceReference>>,
    start: Option<IntOrRange>,
    end: Option<IntOrRange>,
    sequence: Option<SequenceString>,
    meta: Option<Box<Meta>>,
}

impl SequenceLocation {
    /// A closed interval `[start, end)` on a linear sequence.
    ///
    /// # Errors
    /// [`CoordinateError::StartAfterEnd`] if `start > end`; [`CoordinateError::Negative`] for
    /// negative coordinates.
    pub fn new(
        sequence_reference: impl Into<IriOr<SequenceReference>>,
        start: impl Into<IntOrRange>,
        end: impl Into<IntOrRange>,
    ) -> Result<Self, CoordinateError> {
        Self::from_parts(
            Some(sequence_reference.into()),
            Some(start.into()),
            Some(end.into()),
        )
    }

    /// A location defined only by `start`, extending rightwards (increasing coordinates).
    ///
    /// # Errors
    /// [`CoordinateError::Negative`] for negative coordinates.
    pub fn starting_at(
        sequence_reference: impl Into<IriOr<SequenceReference>>,
        start: impl Into<IntOrRange>,
    ) -> Result<Self, CoordinateError> {
        Self::from_parts(Some(sequence_reference.into()), Some(start.into()), None)
    }

    /// A location defined only by `end`, extending leftwards (decreasing coordinates).
    ///
    /// # Errors
    /// [`CoordinateError::Negative`] for negative coordinates.
    pub fn ending_at(
        sequence_reference: impl Into<IriOr<SequenceReference>>,
        end: impl Into<IntOrRange>,
    ) -> Result<Self, CoordinateError> {
        Self::from_parts(Some(sequence_reference.into()), None, Some(end.into()))
    }

    /// General constructor accepting every schema-valid combination.
    ///
    /// # Errors
    /// [`CoordinateError::MissingCoordinates`] if both coordinates are absent;
    /// [`CoordinateError::Negative`] for negative values; [`CoordinateError::StartAfterEnd`]
    /// if `start > end` on a sequence not flagged circular.
    pub fn from_parts(
        sequence_reference: Option<IriOr<SequenceReference>>,
        start: Option<IntOrRange>,
        end: Option<IntOrRange>,
    ) -> Result<Self, CoordinateError> {
        if start.is_none() && end.is_none() {
            return Err(CoordinateError::MissingCoordinates);
        }
        let start = start.map(IntOrRange::require_non_negative).transpose()?;
        let end = end.map(IntOrRange::require_non_negative).transpose()?;
        let circular = sequence_reference
            .as_ref()
            .and_then(IriOr::as_object)
            .is_some_and(SequenceReference::is_circular);
        if !circular
            && let (Some(s), Some(e)) = (start, end)
            && let (Some(lo), Some(hi)) = (s.lower_bound(), e.upper_bound())
            && lo > hi
        {
            return Err(CoordinateError::StartAfterEnd { start: lo, end: hi });
        }
        Ok(Self {
            sequence_reference,
            start,
            end,
            sequence: None,
            meta: None,
        })
    }

    /// Attach the literal reference sequence at this location (the VCF `REF` allele).
    #[must_use]
    pub fn with_sequence(mut self, sequence: SequenceString) -> Self {
        self.sequence = Some(sequence);
        self
    }

    /// `sequenceReference`.
    #[inline]
    pub fn sequence_reference(&self) -> Option<&IriOr<SequenceReference>> {
        self.sequence_reference.as_ref()
    }

    /// The inline sequence reference, if present and not an IRI.
    #[inline]
    pub fn inline_sequence_reference(&self) -> Option<&SequenceReference> {
        self.sequence_reference.as_ref().and_then(IriOr::as_object)
    }

    /// The RefGet accession of the inline sequence reference, if any.
    #[inline]
    pub fn refget_accession(&self) -> Option<&RefgetAccession> {
        self.inline_sequence_reference()
            .map(SequenceReference::refget_accession)
    }

    /// `start`.
    #[inline]
    pub fn start(&self) -> Option<IntOrRange> {
        self.start
    }

    /// `end`.
    #[inline]
    pub fn end(&self) -> Option<IntOrRange> {
        self.end
    }

    /// The literal reference sequence, if carried.
    #[inline]
    pub fn sequence(&self) -> Option<&SequenceString> {
        self.sequence.as_ref()
    }

    /// Exact `(start, end)` when both coordinates are plain integers.
    #[inline]
    pub fn exact_interval(&self) -> Option<(i64, i64)> {
        Some((self.start?.as_int()?, self.end?.as_int()?))
    }

    /// Length of the interval when both coordinates are exact and `start <= end`.
    pub fn exact_length(&self) -> Option<u64> {
        self.exact_interval()
            .and_then(|(s, e)| (e >= s).then(|| (e - s) as u64))
    }

    /// Replace the coordinates, re-validating (used by normalization).
    pub(crate) fn with_coordinates(
        &self,
        start: Option<IntOrRange>,
        end: Option<IntOrRange>,
    ) -> Result<Self, CoordinateError> {
        let mut out = Self::from_parts(self.sequence_reference.clone(), start, end)?;
        out.meta.clone_from(&self.meta);
        // `sequence` describes the old interval; it is intentionally dropped.
        Ok(out)
    }

    pub(crate) fn from_all_parts(
        sequence_reference: Option<IriOr<SequenceReference>>,
        start: Option<IntOrRange>,
        end: Option<IntOrRange>,
        sequence: Option<SequenceString>,
        meta: Option<Box<Meta>>,
    ) -> Result<Self, CoordinateError> {
        let mut out = Self::from_parts(sequence_reference, start, end)?;
        out.sequence = sequence;
        out.meta = meta;
        Ok(out)
    }
}
impl_entity!(SequenceLocation);

/// A VRS `SequenceOffsetLocation` (draft): a location on a mapped sequence expressed as an
/// offset from an anchor position — typically an intronic position relative to an exon
/// boundary on a transcript.
#[derive(Clone, Debug, PartialEq)]
pub struct SequenceOffsetLocation {
    sequence_reference: IriOr<SequenceReference>,
    anchor: i64,
    anchor_orientation: AnchorOrientation,
    offset_start: Option<IntOrRange>,
    offset_end: Option<IntOrRange>,
    meta: Option<Box<Meta>>,
}

impl SequenceOffsetLocation {
    /// An offset interval `[offset_start, offset_end)` relative to `anchor`.
    ///
    /// # Errors
    /// [`CoordinateError::Negative`] if the anchor is negative.
    pub fn new(
        sequence_reference: impl Into<IriOr<SequenceReference>>,
        anchor: i64,
        anchor_orientation: AnchorOrientation,
        offset_start: impl Into<IntOrRange>,
        offset_end: impl Into<IntOrRange>,
    ) -> Result<Self, CoordinateError> {
        Self::from_parts(
            sequence_reference.into(),
            anchor,
            anchor_orientation,
            Some(offset_start.into()),
            Some(offset_end.into()),
        )
    }

    /// General constructor (offsets optional, as in the schema).
    ///
    /// # Errors
    /// [`CoordinateError::Negative`] if the anchor is negative.
    pub fn from_parts(
        sequence_reference: IriOr<SequenceReference>,
        anchor: i64,
        anchor_orientation: AnchorOrientation,
        offset_start: Option<IntOrRange>,
        offset_end: Option<IntOrRange>,
    ) -> Result<Self, CoordinateError> {
        if anchor < 0 {
            return Err(CoordinateError::Negative(anchor));
        }
        Ok(Self {
            sequence_reference,
            anchor,
            anchor_orientation,
            offset_start,
            offset_end,
            meta: None,
        })
    }

    /// `sequenceReference`.
    #[inline]
    pub fn sequence_reference(&self) -> &IriOr<SequenceReference> {
        &self.sequence_reference
    }

    /// `anchor`: the inter-residue position the offsets are measured from.
    #[inline]
    pub fn anchor(&self) -> i64 {
        self.anchor
    }

    /// `anchorOrientation`.
    #[inline]
    pub fn anchor_orientation(&self) -> AnchorOrientation {
        self.anchor_orientation
    }

    /// `offsetStart`.
    #[inline]
    pub fn offset_start(&self) -> Option<IntOrRange> {
        self.offset_start
    }

    /// `offsetEnd`.
    #[inline]
    pub fn offset_end(&self) -> Option<IntOrRange> {
        self.offset_end
    }

    /// The largest absolute offset magnitude among the defined bounds of `offsetStart` and
    /// `offsetEnd` (the quantity compared by relative-allele anchor selection).
    pub fn max_offset_magnitude(&self) -> Option<u64> {
        [self.offset_start, self.offset_end]
            .into_iter()
            .flatten()
            .flat_map(|o| [o.lower_bound(), o.upper_bound()])
            .flatten()
            .map(i64::unsigned_abs)
            .max()
    }

    pub(crate) fn set_meta(&mut self, meta: Option<Box<Meta>>) {
        self.meta = meta;
    }
}
impl_entity!(SequenceOffsetLocation);

/// A VRS `RelativeSequenceLocation` (draft): an absolute base location paired with its
/// position relative to an anchor on a mapped sequence.
#[derive(Clone, Debug, PartialEq)]
pub struct RelativeSequenceLocation {
    base_sequence_location: IriOr<SequenceLocation>,
    mapped_sequence_location: IriOr<SequenceOffsetLocation>,
    meta: Option<Box<Meta>>,
}

impl RelativeSequenceLocation {
    /// Pair a base location with its mapped offset location.
    pub fn new(
        base_sequence_location: impl Into<IriOr<SequenceLocation>>,
        mapped_sequence_location: impl Into<IriOr<SequenceOffsetLocation>>,
    ) -> Self {
        Self {
            base_sequence_location: base_sequence_location.into(),
            mapped_sequence_location: mapped_sequence_location.into(),
            meta: None,
        }
    }

    /// `baseSequenceLocation`.
    #[inline]
    pub fn base_sequence_location(&self) -> &IriOr<SequenceLocation> {
        &self.base_sequence_location
    }

    /// `mappedSequenceLocation`.
    #[inline]
    pub fn mapped_sequence_location(&self) -> &IriOr<SequenceOffsetLocation> {
        &self.mapped_sequence_location
    }

    pub(crate) fn set_meta(&mut self, meta: Option<Box<Meta>>) {
        self.meta = meta;
    }
}
impl_entity!(RelativeSequenceLocation);

/// A VRS `Location`: any contiguous segment of a biological sequence.
///
/// The (draft, rare) relative variant is boxed so that `Location` stays the size of a
/// `SequenceLocation`.
#[derive(Clone, Debug, PartialEq)]
pub enum Location {
    /// An interval on a sequence reference.
    Sequence(SequenceLocation),
    /// A base location with a mapped relative position (draft).
    RelativeSequence(Box<RelativeSequenceLocation>),
}

impl Location {
    /// The sequence location, if this is one.
    #[inline]
    pub fn as_sequence_location(&self) -> Option<&SequenceLocation> {
        match self {
            Self::Sequence(l) => Some(l),
            Self::RelativeSequence(_) => None,
        }
    }

    /// The relative sequence location, if this is one.
    #[inline]
    pub fn as_relative_sequence_location(&self) -> Option<&RelativeSequenceLocation> {
        match self {
            Self::Sequence(_) => None,
            Self::RelativeSequence(l) => Some(l),
        }
    }

    /// The VRS class name of the concrete location.
    pub const fn type_name(&self) -> &'static str {
        match self {
            Self::Sequence(_) => "SequenceLocation",
            Self::RelativeSequence(_) => "RelativeSequenceLocation",
        }
    }
}

impl From<SequenceLocation> for Location {
    fn from(l: SequenceLocation) -> Self {
        Self::Sequence(l)
    }
}
impl From<RelativeSequenceLocation> for Location {
    fn from(l: RelativeSequenceLocation) -> Self {
        Self::RelativeSequence(Box::new(l))
    }
}
impl From<Box<RelativeSequenceLocation>> for Location {
    fn from(l: Box<RelativeSequenceLocation>) -> Self {
        Self::RelativeSequence(l)
    }
}
impl From<SequenceLocation> for IriOr<Location> {
    fn from(l: SequenceLocation) -> Self {
        IriOr::Object(Location::Sequence(l))
    }
}
impl From<RelativeSequenceLocation> for IriOr<Location> {
    fn from(l: RelativeSequenceLocation) -> Self {
        IriOr::Object(Location::from(l))
    }
}
