//! Error types.
//!
//! Every fallible operation in this crate returns a specific error enum so that callers can
//! distinguish *why* something failed: malformed coordinates, an invalid sequence string, a bad
//! identifier, a structural model violation, a serialization problem, a missing reference
//! sequence, or an unsupported specification feature.  The umbrella [`Error`] type is provided
//! for callers that want a single error type; all specific errors convert into it.

/// Errors raised while constructing or validating genomic coordinates and ranges.
#[derive(Debug, Clone, PartialEq, Eq, thiserror::Error)]
#[non_exhaustive]
pub enum CoordinateError {
    /// A coordinate (or range bound) was negative; VRS coordinates are inter-residue
    /// positions starting at 0.
    #[error("coordinate {0} is negative; VRS inter-residue coordinates start at 0")]
    Negative(i64),
    /// A `Range` had `null` for both bounds.
    #[error("a Range must have at least one integer bound")]
    UnboundedRange,
    /// A `Range` whose minimum exceeds its maximum.
    #[error("Range minimum {min} is greater than maximum {max}")]
    InvertedRange {
        /// The offending lower bound.
        min: i64,
        /// The offending upper bound.
        max: i64,
    },
    /// `start` is after `end` on a linear sequence.
    #[error("start ({start}) is greater than end ({end}) on a linear sequence")]
    StartAfterEnd {
        /// The start coordinate (lower bound if a range).
        start: i64,
        /// The end coordinate (upper bound if a range).
        end: i64,
    },
    /// A `SequenceLocation` with neither `start` nor `end`.
    #[error("a SequenceLocation must define at least one of start and end")]
    MissingCoordinates,
    /// A value is too large to be represented (bounds are `i64`; the extreme values
    /// `i64::MIN` and `i64::MAX` are reserved).
    #[error("integer value {0} is out of the supported range")]
    OutOfRange(i128),
}

/// Errors raised while validating a VRS `sequenceString` (`^[A-Z*\-]*$`).
#[derive(Debug, Clone, PartialEq, Eq, thiserror::Error)]
#[non_exhaustive]
pub enum SequenceStringError {
    /// A byte that is not an upper-case ASCII letter, `*` or `-`.
    #[error("invalid residue byte {byte:#04x} at offset {offset}; expected [A-Z*-]")]
    InvalidResidue {
        /// Offending byte.
        byte: u8,
        /// Zero-based offset in the input.
        offset: usize,
    },
}

/// Errors raised while parsing identifiers: RefGet accessions, digests and GA4GH CURIEs.
#[derive(Debug, Clone, PartialEq, Eq, thiserror::Error)]
#[non_exhaustive]
pub enum IdentifierError {
    /// A RefGet accession that is not `SQ.` followed by 32 base64url characters.
    #[error("invalid RefGet accession {0:?}; expected `SQ.` followed by 32 base64url characters")]
    InvalidRefgetAccession(String),
    /// A digest that is not exactly 32 base64url characters.
    #[error("invalid sha512t24u digest {0:?}; expected 32 base64url characters")]
    InvalidDigest(String),
    /// A GA4GH identifier that is not `ga4gh:<prefix>.<digest>`.
    #[error("invalid GA4GH identifier {0:?}; expected `ga4gh:<type prefix>.<32-char digest>`")]
    InvalidGa4ghIdentifier(String),
    /// An unknown VRS type prefix.
    #[error("unknown VRS type prefix {0:?}")]
    UnknownTypePrefix(String),
}

/// Errors raised by structural rules of the VRS information model that are not expressible
/// through the type system alone.
#[derive(Debug, Clone, PartialEq, Eq, thiserror::Error)]
#[non_exhaustive]
pub enum ModelError {
    /// A collection had fewer members than the schema requires.
    #[error("{class}.{property} requires at least {min} items, got {actual}")]
    TooFewItems {
        /// Class name (e.g. `CisPhasedBlock`).
        class: &'static str,
        /// Property name (e.g. `members`).
        property: &'static str,
        /// Minimum required by the schema.
        min: usize,
        /// Number supplied.
        actual: usize,
    },
    /// An adjoined sequence of an `Adjacency` defined both `start` and `end`.
    #[error("an Adjacency adjoined sequence must not define both start and end")]
    AdjoinedSequenceHasStartAndEnd,
    /// A string was not a member of a closed VRS value set.
    #[error("{value:?} is not a valid {value_set}")]
    UnknownEnumValue {
        /// Name of the value set (e.g. `copyChange`).
        value_set: &'static str,
        /// The offending value.
        value: String,
    },
    /// A JSON `type` discriminator did not match the expected class.
    #[error("expected type {expected:?}, found {found:?}")]
    TypeMismatch {
        /// The class the caller expected.
        expected: &'static str,
        /// The class named in the data.
        found: String,
    },
    /// A count (copies, repeat subunit length, length) was negative.
    #[error("{property} must be non-negative, got {value}")]
    NegativeCount {
        /// Property name.
        property: &'static str,
        /// Offending value.
        value: i64,
    },
}

/// Errors raised while fetching reference sequence data from a [`SequenceProvider`].
///
/// [`SequenceProvider`]: crate::normalize::SequenceProvider
#[derive(Debug, thiserror::Error)]
#[non_exhaustive]
pub enum SequenceError {
    /// The provider has no sequence with this RefGet accession.
    #[error("unknown sequence {0}")]
    UnknownSequence(String),
    /// The requested interval is not covered by the provider.
    #[error("interval [{start}, {end}) is outside the available range of {accession}")]
    OutOfBounds {
        /// RefGet accession.
        accession: String,
        /// Requested start (inter-residue).
        start: u64,
        /// Requested end (inter-residue).
        end: u64,
    },
    /// The provider returned residues that are not valid `sequenceString` bytes.
    #[error("provider returned an invalid sequence: {0}")]
    InvalidSequence(#[from] SequenceStringError),
    /// Any other backend failure (I/O, network, database).
    #[error("sequence backend error: {0}")]
    Backend(#[source] Box<dyn std::error::Error + Send + Sync + 'static>),
}

impl SequenceError {
    /// Wrap an arbitrary backend error.
    pub fn backend<E>(err: E) -> Self
    where
        E: std::error::Error + Send + Sync + 'static,
    {
        Self::Backend(Box::new(err))
    }
}

/// A VRS feature that this implementation does not (yet) support in the requested context.
#[derive(Debug, Clone, PartialEq, Eq, thiserror::Error)]
#[error("unsupported VRS feature: {0}")]
pub struct UnsupportedError(pub String);

impl UnsupportedError {
    /// Create an unsupported-feature error.
    pub fn new(what: impl Into<String>) -> Self {
        Self(what.into())
    }
}

/// Errors raised by normalization.
#[derive(Debug, thiserror::Error)]
#[non_exhaustive]
pub enum NormalizeError {
    /// Reference sequence data could not be obtained.
    #[error(transparent)]
    Sequence(#[from] SequenceError),
    /// The resulting object violated a coordinate invariant (reported rather than panicking).
    #[error(transparent)]
    Coordinate(#[from] CoordinateError),
    /// The resulting object violated a model invariant.
    #[error(transparent)]
    Model(#[from] ModelError),
    /// The object uses a feature the normalizer cannot handle (e.g. circular sequences).
    #[error(transparent)]
    Unsupported(#[from] UnsupportedError),
}

/// Errors raised when producing legacy (VRS 1.3) digests.
#[derive(Debug, Clone, PartialEq, Eq, thiserror::Error)]
#[non_exhaustive]
pub enum LegacyDigestError {
    /// VRS 1.3 serialization needs an inline `SequenceReference`.
    #[error("VRS 1.3 serialization requires an inline SequenceReference with a RefGet accession")]
    MissingSequenceReference,
    /// VRS 1.3 had no null coordinates.
    #[error("VRS 1.3 serialization requires both start and end")]
    MissingCoordinates,
    /// VRS 1.3 alleles carried a literal sequence.
    #[error("VRS 1.3 serialization requires a literal state sequence")]
    MissingStateSequence,
    /// A state type that did not exist in VRS 1.3.
    #[error("VRS 1.3 does not support {0} states")]
    UnsupportedState(&'static str),
    /// A location type that did not exist in VRS 1.3.
    #[error("VRS 1.3 does not support {0} locations")]
    UnsupportedLocation(&'static str),
}

/// Coarse classification of a JSON error (mirrors `serde_json::error::Category`).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum JsonErrorKind {
    /// Malformed JSON text.
    Syntax,
    /// Well-formed JSON that does not fit the VRS model (missing field, bad enum value, ...).
    Data,
    /// Unexpected end of input.
    Eof,
    /// I/O failure.
    Io,
}

/// JSON (de)serialization error.
#[derive(Debug, thiserror::Error)]
#[error(transparent)]
pub struct JsonError(#[from] pub serde_json::Error);

impl JsonError {
    /// Classify the error.
    pub fn kind(&self) -> JsonErrorKind {
        match self.0.classify() {
            serde_json::error::Category::Syntax => JsonErrorKind::Syntax,
            serde_json::error::Category::Data => JsonErrorKind::Data,
            serde_json::error::Category::Eof => JsonErrorKind::Eof,
            serde_json::error::Category::Io => JsonErrorKind::Io,
        }
    }
}

/// Umbrella error type covering every failure mode of this crate.
#[derive(Debug, thiserror::Error)]
#[non_exhaustive]
pub enum Error {
    /// Invalid genomic coordinates.
    #[error(transparent)]
    Coordinate(#[from] CoordinateError),
    /// Invalid sequence string.
    #[error(transparent)]
    SequenceString(#[from] SequenceStringError),
    /// Invalid identifier.
    #[error(transparent)]
    Identifier(#[from] IdentifierError),
    /// Model-level rule violation.
    #[error(transparent)]
    Model(#[from] ModelError),
    /// Reference sequence access failure.
    #[error(transparent)]
    Sequence(#[from] SequenceError),
    /// Normalization failure.
    #[error(transparent)]
    Normalize(#[from] NormalizeError),
    /// Legacy digest failure.
    #[error(transparent)]
    LegacyDigest(#[from] LegacyDigestError),
    /// JSON error.
    #[error(transparent)]
    Json(#[from] JsonError),
    /// Unsupported feature.
    #[error(transparent)]
    Unsupported(#[from] UnsupportedError),
}

impl From<serde_json::Error> for Error {
    fn from(e: serde_json::Error) -> Self {
        Self::Json(JsonError(e))
    }
}

/// Convenience alias using the umbrella [`Error`].
pub type Result<T, E = Error> = std::result::Result<T, E>;
