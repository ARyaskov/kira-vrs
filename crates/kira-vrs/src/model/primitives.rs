//! Primitive VRS value types: sequence strings, ranges, IRI references and closed value sets.

use std::fmt;
use std::str::FromStr;

use crate::error::{CoordinateError, ModelError, SequenceStringError};

// ---------------------------------------------------------------------------------------------
// SequenceString
// ---------------------------------------------------------------------------------------------

/// Inline capacity of [`SequenceString`]; sequences up to this length never touch the heap.
const INLINE_CAP: usize = 22;

/// A VRS `sequenceString`: residues in `[A-Z*-]`, in conventional order (5'→3' or N→C).
///
/// Sequences up to 22 residues (every SNV, most indels) are stored inline without any heap
/// allocation; longer sequences are boxed.  Comparison, hashing and ordering are by content.
#[derive(Clone)]
pub struct SequenceString(Repr);

#[derive(Clone)]
enum Repr {
    Inline { len: u8, buf: [u8; INLINE_CAP] },
    Heap(Box<[u8]>),
}

impl SequenceString {
    /// The empty sequence (used for deletions).
    pub const EMPTY: Self = Self(Repr::Inline {
        len: 0,
        buf: [0; INLINE_CAP],
    });

    /// Validate and construct a sequence string.
    ///
    /// # Errors
    /// Returns [`SequenceStringError::InvalidResidue`] if any byte is outside `[A-Z*-]`.
    pub fn new(s: &str) -> Result<Self, SequenceStringError> {
        Self::from_bytes(s.as_bytes())
    }

    /// Validate and construct a sequence string from raw bytes.
    ///
    /// # Errors
    /// Returns [`SequenceStringError::InvalidResidue`] if any byte is outside `[A-Z*-]`.
    pub fn from_bytes(bytes: &[u8]) -> Result<Self, SequenceStringError> {
        validate_residues(bytes)?;
        Ok(Self::from_valid_bytes(bytes))
    }

    /// Construct from bytes already known to be valid residues (crate-internal fast path).
    pub(crate) fn from_valid_bytes(bytes: &[u8]) -> Self {
        debug_assert!(validate_residues(bytes).is_ok());
        if bytes.len() <= INLINE_CAP {
            let mut buf = [0u8; INLINE_CAP];
            buf[..bytes.len()].copy_from_slice(bytes);
            Self(Repr::Inline {
                len: bytes.len() as u8,
                buf,
            })
        } else {
            Self(Repr::Heap(bytes.into()))
        }
    }

    /// The residues as bytes.
    #[inline]
    pub fn as_bytes(&self) -> &[u8] {
        match &self.0 {
            Repr::Inline { len, buf } => &buf[..usize::from(*len)],
            Repr::Heap(b) => b,
        }
    }

    /// The residues as a string slice (always ASCII).
    #[inline]
    pub fn as_str(&self) -> &str {
        // Residues are validated ASCII, so this cannot fail.
        std::str::from_utf8(self.as_bytes()).unwrap_or_default()
    }

    /// Number of residues.
    #[inline]
    pub fn len(&self) -> usize {
        self.as_bytes().len()
    }

    /// `true` for the empty sequence.
    #[inline]
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }
}

#[inline]
pub(crate) const fn is_residue(b: u8) -> bool {
    b.is_ascii_uppercase() || b == b'*' || b == b'-'
}

/// Validate that every byte is a residue in `[A-Z*-]`.
pub(crate) fn validate_residues(bytes: &[u8]) -> Result<(), SequenceStringError> {
    match bytes.iter().position(|&b| !is_residue(b)) {
        None => Ok(()),
        Some(offset) => Err(SequenceStringError::InvalidResidue {
            byte: bytes[offset],
            offset,
        }),
    }
}

impl PartialEq for SequenceString {
    fn eq(&self, other: &Self) -> bool {
        self.as_bytes() == other.as_bytes()
    }
}
impl Eq for SequenceString {}
impl PartialOrd for SequenceString {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}
impl Ord for SequenceString {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        self.as_bytes().cmp(other.as_bytes())
    }
}
impl std::hash::Hash for SequenceString {
    fn hash<H: std::hash::Hasher>(&self, state: &mut H) {
        self.as_bytes().hash(state);
    }
}
impl fmt::Debug for SequenceString {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "SequenceString({:?})", self.as_str())
    }
}
impl fmt::Display for SequenceString {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(self.as_str())
    }
}
impl Default for SequenceString {
    fn default() -> Self {
        Self::EMPTY
    }
}
impl AsRef<[u8]> for SequenceString {
    fn as_ref(&self) -> &[u8] {
        self.as_bytes()
    }
}
impl AsRef<str> for SequenceString {
    fn as_ref(&self) -> &str {
        self.as_str()
    }
}
impl FromStr for SequenceString {
    type Err = SequenceStringError;
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        Self::new(s)
    }
}
impl TryFrom<&str> for SequenceString {
    type Error = SequenceStringError;
    fn try_from(s: &str) -> Result<Self, Self::Error> {
        Self::new(s)
    }
}
impl TryFrom<String> for SequenceString {
    type Error = SequenceStringError;
    fn try_from(s: String) -> Result<Self, Self::Error> {
        Self::new(&s)
    }
}
impl TryFrom<&[u8]> for SequenceString {
    type Error = SequenceStringError;
    fn try_from(b: &[u8]) -> Result<Self, Self::Error> {
        Self::from_bytes(b)
    }
}

// ---------------------------------------------------------------------------------------------
// Range and IntOrRange
// ---------------------------------------------------------------------------------------------

/// Sentinel for an absent lower bound. Never a valid bound value.
const NO_MIN: i64 = i64::MIN;
/// Sentinel for an absent upper bound. Never a valid bound value.
const NO_MAX: i64 = i64::MAX;

/// A VRS `Range`: an inclusive range bounded by one or two integers (`[min, max]`, either of
/// which may be unbounded, but not both).
///
/// Stored in 16 bytes; the extreme values `i64::MIN` / `i64::MAX` are reserved and rejected.
#[derive(Clone, Copy, PartialEq, Eq, Hash)]
pub struct Range {
    lo: i64,
    hi: i64,
}

impl Range {
    /// Construct a range from optional bounds.
    ///
    /// # Errors
    /// [`CoordinateError::UnboundedRange`] if both are `None`;
    /// [`CoordinateError::InvertedRange`] if `min > max`;
    /// [`CoordinateError::OutOfRange`] for the reserved extreme values.
    pub fn new(min: Option<i64>, max: Option<i64>) -> Result<Self, CoordinateError> {
        let lo = match min {
            Some(v) => checked_bound(v)?,
            None => NO_MIN,
        };
        let hi = match max {
            Some(v) => checked_bound(v)?,
            None => NO_MAX,
        };
        match (min, max) {
            (None, None) => Err(CoordinateError::UnboundedRange),
            (Some(a), Some(b)) if a > b => Err(CoordinateError::InvertedRange { min: a, max: b }),
            _ => Ok(Self { lo, hi }),
        }
    }

    /// A definite range `[min, max]`.
    ///
    /// # Errors
    /// See [`Range::new`].
    pub fn bounded(min: i64, max: i64) -> Result<Self, CoordinateError> {
        Self::new(Some(min), Some(max))
    }

    /// An indefinite range `[min, null]` ("at least `min`").
    ///
    /// # Errors
    /// See [`Range::new`].
    pub fn at_least(min: i64) -> Result<Self, CoordinateError> {
        Self::new(Some(min), None)
    }

    /// An indefinite range `[null, max]` ("at most `max`").
    ///
    /// # Errors
    /// See [`Range::new`].
    pub fn at_most(max: i64) -> Result<Self, CoordinateError> {
        Self::new(None, Some(max))
    }

    /// Lower bound, if bounded.
    #[inline]
    pub fn min(&self) -> Option<i64> {
        (self.lo != NO_MIN).then_some(self.lo)
    }

    /// Upper bound, if bounded.
    #[inline]
    pub fn max(&self) -> Option<i64> {
        (self.hi != NO_MAX).then_some(self.hi)
    }

    /// `true` when both bounds are present.
    #[inline]
    pub fn is_definite(&self) -> bool {
        self.lo != NO_MIN && self.hi != NO_MAX
    }

    /// `true` if any bound is negative.
    pub(crate) fn has_negative_bound(&self) -> bool {
        self.min().is_some_and(|v| v < 0) || self.max().is_some_and(|v| v < 0)
    }
}

fn checked_bound(v: i64) -> Result<i64, CoordinateError> {
    if v == NO_MIN || v == NO_MAX {
        Err(CoordinateError::OutOfRange(i128::from(v)))
    } else {
        Ok(v)
    }
}

impl fmt::Debug for Range {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "Range[{:?}, {:?}]", self.min(), self.max())
    }
}

impl fmt::Display for Range {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str("[")?;
        match self.min() {
            Some(v) => write!(f, "{v}")?,
            None => f.write_str("null")?,
        }
        f.write_str(", ")?;
        match self.max() {
            Some(v) => write!(f, "{v}")?,
            None => f.write_str("null")?,
        }
        f.write_str("]")
    }
}

/// An integer or a [`Range`] — the VRS type of `start`, `end`, `copies`, `length` and offsets.
#[derive(Clone, Copy, PartialEq, Eq, Hash, Debug)]
pub enum IntOrRange {
    /// An exact integer.
    Int(i64),
    /// An inclusive range.
    Range(Range),
}

impl IntOrRange {
    /// The exact integer, if not a range.
    #[inline]
    pub fn as_int(&self) -> Option<i64> {
        match self {
            Self::Int(v) => Some(*v),
            Self::Range(_) => None,
        }
    }

    /// The range, if not an exact integer.
    #[inline]
    pub fn as_range(&self) -> Option<&Range> {
        match self {
            Self::Int(_) => None,
            Self::Range(r) => Some(r),
        }
    }

    /// Lower bound (the value itself for an integer; `None` if unbounded below).
    #[inline]
    pub fn lower_bound(&self) -> Option<i64> {
        match self {
            Self::Int(v) => Some(*v),
            Self::Range(r) => r.min(),
        }
    }

    /// Upper bound (the value itself for an integer; `None` if unbounded above).
    #[inline]
    pub fn upper_bound(&self) -> Option<i64> {
        match self {
            Self::Int(v) => Some(*v),
            Self::Range(r) => r.max(),
        }
    }

    /// `true` if the value or any bound is negative.
    pub(crate) fn has_negative(&self) -> bool {
        match self {
            Self::Int(v) => *v < 0,
            Self::Range(r) => r.has_negative_bound(),
        }
    }

    /// Reject negative values, reporting the first negative bound.
    pub(crate) fn require_non_negative(self) -> Result<Self, CoordinateError> {
        let negative = [self.lower_bound(), self.upper_bound()]
            .into_iter()
            .flatten()
            .find(|v| *v < 0);
        match negative {
            Some(v) => Err(CoordinateError::Negative(v)),
            None => Ok(self),
        }
    }
}

impl From<i64> for IntOrRange {
    fn from(v: i64) -> Self {
        Self::Int(v)
    }
}
impl From<i32> for IntOrRange {
    fn from(v: i32) -> Self {
        Self::Int(i64::from(v))
    }
}
impl From<u32> for IntOrRange {
    fn from(v: u32) -> Self {
        Self::Int(i64::from(v))
    }
}
impl TryFrom<u64> for IntOrRange {
    type Error = CoordinateError;
    fn try_from(v: u64) -> Result<Self, Self::Error> {
        i64::try_from(v)
            .ok()
            .filter(|v| *v != NO_MAX)
            .map(Self::Int)
            .ok_or(CoordinateError::OutOfRange(i128::from(v)))
    }
}
impl TryFrom<usize> for IntOrRange {
    type Error = CoordinateError;
    fn try_from(v: usize) -> Result<Self, Self::Error> {
        Self::try_from(v as u64)
    }
}
impl From<Range> for IntOrRange {
    fn from(r: Range) -> Self {
        Self::Range(r)
    }
}

impl fmt::Display for IntOrRange {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Int(v) => write!(f, "{v}"),
            Self::Range(r) => write!(f, "{r}"),
        }
    }
}

// ---------------------------------------------------------------------------------------------
// IRI references
// ---------------------------------------------------------------------------------------------

/// A gkm-core `iriReference`: an IRI or relative reference (RFC 3987 §2.1) used to point at
/// an object defined elsewhere, e.g. `ga4gh:SL.4t6JnYWqHwYw9WzBT_lmWBb3tLQNalkT` or
/// `locations.json#/0`.
#[derive(Clone, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct Iri(Box<str>);

impl Iri {
    /// Wrap a string as an IRI reference. No syntactic validation is performed beyond what the
    /// schema requires (any string).
    pub fn new(iri: impl Into<Box<str>>) -> Self {
        Self(iri.into())
    }

    /// The IRI text.
    #[inline]
    pub fn as_str(&self) -> &str {
        &self.0
    }

    /// If this IRI is a GA4GH computed identifier (`ga4gh:<prefix>.<32-char digest>`), return
    /// the digest component. This is how IRIs are represented in digest serialization.
    pub fn ga4gh_digest(&self) -> Option<&str> {
        let rest = self.0.strip_prefix("ga4gh:")?;
        let dot = rest.find('.')?;
        let (prefix, digest) = (&rest[..dot], &rest[dot + 1..]);
        if prefix.is_empty() || prefix.contains('.') {
            return None;
        }
        crate::model::identifier::is_digest_str(digest).then_some(digest)
    }
}

impl fmt::Debug for Iri {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "Iri({:?})", &*self.0)
    }
}
impl fmt::Display for Iri {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(&self.0)
    }
}
impl From<&str> for Iri {
    fn from(s: &str) -> Self {
        Self::new(s)
    }
}
impl From<String> for Iri {
    fn from(s: String) -> Self {
        Self::new(s)
    }
}
impl AsRef<str> for Iri {
    fn as_ref(&self) -> &str {
        &self.0
    }
}

/// A property that holds either an inline object or an [`Iri`] reference to one.
///
/// VRS allows most nested identifiable objects (locations, alleles, sequence references) to
/// be replaced by an IRI, typically a GA4GH computed identifier.
#[derive(Clone, PartialEq, Eq, Hash, Debug)]
pub enum IriOr<T> {
    /// A reference to an object defined elsewhere.
    Iri(Iri),
    /// The object itself.
    Object(T),
}

impl<T> IriOr<T> {
    /// The inline object, if present.
    #[inline]
    pub fn as_object(&self) -> Option<&T> {
        match self {
            Self::Object(t) => Some(t),
            Self::Iri(_) => None,
        }
    }

    /// Mutable access to the inline object, if present.
    #[inline]
    pub fn as_object_mut(&mut self) -> Option<&mut T> {
        match self {
            Self::Object(t) => Some(t),
            Self::Iri(_) => None,
        }
    }

    /// The IRI, if this is a reference.
    #[inline]
    pub fn as_iri(&self) -> Option<&Iri> {
        match self {
            Self::Iri(i) => Some(i),
            Self::Object(_) => None,
        }
    }

    /// Consume, returning the inline object if present.
    #[inline]
    pub fn into_object(self) -> Option<T> {
        match self {
            Self::Object(t) => Some(t),
            Self::Iri(_) => None,
        }
    }

    /// Map the inline object.
    pub fn map<U>(self, f: impl FnOnce(T) -> U) -> IriOr<U> {
        match self {
            Self::Object(t) => IriOr::Object(f(t)),
            Self::Iri(i) => IriOr::Iri(i),
        }
    }
}

impl<T> From<T> for IriOr<T> {
    fn from(t: T) -> Self {
        Self::Object(t)
    }
}
/// Implements `From<Iri>` for `IriOr<T>` for concrete `T` (a blanket impl would overlap
/// `From<T>` when `T = Iri`).
macro_rules! iri_or_from_iri {
    ($($ty:ty),+ $(,)?) => {
        $(
            impl From<Iri> for IriOr<$ty> {
                fn from(i: Iri) -> Self {
                    Self::Iri(i)
                }
            }
        )+
    };
}
pub(crate) use iri_or_from_iri;

// ---------------------------------------------------------------------------------------------
// Closed value sets
// ---------------------------------------------------------------------------------------------

macro_rules! value_set {
    (
        $(#[$meta:meta])*
        $name:ident, $set_name:literal, {
            $( $(#[$vmeta:meta])* $variant:ident => $text:literal ),+ $(,)?
        }
    ) => {
        $(#[$meta])*
        #[derive(Clone, Copy, PartialEq, Eq, Hash, Debug, PartialOrd, Ord)]
        pub enum $name {
            $( $(#[$vmeta])* $variant, )+
        }

        impl $name {
            /// All members of the value set, in schema order.
            pub const ALL: &'static [Self] = &[ $( Self::$variant, )+ ];

            /// The wire representation of this value.
            #[inline]
            pub const fn as_str(self) -> &'static str {
                match self {
                    $( Self::$variant => $text, )+
                }
            }

            /// Parse the wire representation.
            ///
            /// # Errors
            /// [`ModelError::UnknownEnumValue`] for strings outside the value set.
            pub fn parse(s: &str) -> Result<Self, ModelError> {
                match s {
                    $( $text => Ok(Self::$variant), )+
                    _ => Err(ModelError::UnknownEnumValue {
                        value_set: $set_name,
                        value: s.to_owned(),
                    }),
                }
            }
        }

        impl FromStr for $name {
            type Err = ModelError;
            fn from_str(s: &str) -> Result<Self, Self::Err> {
                Self::parse(s)
            }
        }

        impl fmt::Display for $name {
            fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
                f.write_str(self.as_str())
            }
        }
    };
}

value_set! {
    /// Interpretation of the residue characters of a sequence (`SequenceReference.residueAlphabet`).
    ResidueAlphabet, "residueAlphabet", {
        /// Amino acid character set.
        AminoAcid => "aa",
        /// Nucleic acid character set.
        NucleicAcid => "na",
    }
}

value_set! {
    /// RefSeq molecule types (`SequenceReference.moleculeType`).
    MoleculeType, "moleculeType", {
        /// Genomic DNA.
        Genomic => "genomic",
        /// RNA.
        Rna => "RNA",
        /// Messenger RNA.
        MRna => "mRNA",
        /// Protein.
        Protein => "protein",
    }
}

value_set! {
    /// EFO-derived copy number change categories (`CopyNumberChange.copyChange`).
    CopyChange, "copyChange", {
        /// EFO:0030069 complete genomic deletion.
        CompleteGenomicLoss => "complete genomic loss",
        /// EFO:0020073 high-level copy number loss.
        HighLevelLoss => "high-level loss",
        /// EFO:0030068 low-level copy number loss.
        LowLevelLoss => "low-level loss",
        /// EFO:0030067 copy number loss.
        Loss => "loss",
        /// EFO:0030064 regional base ploidy.
        RegionalBasePloidy => "regional base ploidy",
        /// EFO:0030070 copy number gain.
        Gain => "gain",
        /// EFO:0030071 low-level copy number gain.
        LowLevelGain => "low-level gain",
        /// EFO:0030072 high-level copy number gain.
        HighLevelGain => "high-level gain",
    }
}

value_set! {
    /// Nomenclatures accepted by `Expression.syntax`.
    Syntax, "syntax", {
        /// HGVS coding DNA.
        HgvsC => "hgvs.c",
        /// HGVS protein.
        HgvsP => "hgvs.p",
        /// HGVS genomic.
        HgvsG => "hgvs.g",
        /// HGVS mitochondrial.
        HgvsM => "hgvs.m",
        /// HGVS non-coding.
        HgvsN => "hgvs.n",
        /// HGVS RNA.
        HgvsR => "hgvs.r",
        /// ISCN cytogenetic nomenclature.
        Iscn => "iscn",
        /// gnomAD variant identifier.
        Gnomad => "gnomad",
        /// NCBI SPDI.
        Spdi => "spdi",
    }
}

value_set! {
    /// Orientation of a `TraversalBlock` component.
    Orientation, "orientation", {
        /// Forward strand.
        Forward => "forward",
        /// Reverse complement.
        ReverseComplement => "reverse_complement",
    }
}

value_set! {
    /// Which side of a discontinuous anchor a `SequenceOffsetLocation` is measured from.
    AnchorOrientation, "anchorOrientation", {
        /// The side immediately preceding the anchor in sequence-reference coordinate order.
        Left => "left",
        /// The side immediately following the anchor in sequence-reference coordinate order.
        Right => "right",
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn sequence_string_inline_and_heap() {
        let s = SequenceString::new("ACGT").unwrap();
        assert_eq!(s.as_str(), "ACGT");
        assert_eq!(s.len(), 4);
        assert!(matches!(s.0, Repr::Inline { .. }));
        let long = "A".repeat(100);
        let l = SequenceString::new(&long).unwrap();
        assert!(matches!(l.0, Repr::Heap(_)));
        assert_eq!(l.as_str(), long);
        assert_eq!(std::mem::size_of::<SequenceString>(), 24);
    }

    #[test]
    fn sequence_string_rejects_lowercase() {
        let err = SequenceString::new("acgt").unwrap_err();
        assert_eq!(
            err,
            SequenceStringError::InvalidResidue {
                byte: b'a',
                offset: 0
            }
        );
        assert!(SequenceString::new("AC*-N").is_ok());
        assert!(SequenceString::new("").is_ok());
    }

    #[test]
    fn range_rules() {
        assert_eq!(Range::new(None, None), Err(CoordinateError::UnboundedRange));
        assert_eq!(
            Range::bounded(5, 3),
            Err(CoordinateError::InvertedRange { min: 5, max: 3 })
        );
        let r = Range::at_least(3).unwrap();
        assert_eq!((r.min(), r.max()), (Some(3), None));
        assert!(!r.is_definite());
        assert_eq!(std::mem::size_of::<Range>(), 16);
        assert_eq!(std::mem::size_of::<Option<IntOrRange>>(), 24);
    }

    #[test]
    fn iri_ga4gh_digest() {
        let iri = Iri::new("ga4gh:SL.4t6JnYWqHwYw9WzBT_lmWBb3tLQNalkT");
        assert_eq!(iri.ga4gh_digest(), Some("4t6JnYWqHwYw9WzBT_lmWBb3tLQNalkT"));
        assert_eq!(Iri::new("refseq:NM_000551.3").ga4gh_digest(), None);
        assert_eq!(Iri::new("ga4gh:SL.short").ga4gh_digest(), None);
    }

    #[test]
    fn value_sets_round_trip() {
        for v in CopyChange::ALL {
            assert_eq!(CopyChange::parse(v.as_str()).unwrap(), *v);
        }
        assert!(Syntax::parse("hgvs.x").is_err());
        assert_eq!(MoleculeType::MRna.as_str(), "mRNA");
    }
}
