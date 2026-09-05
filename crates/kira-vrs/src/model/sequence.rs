//! `SequenceReference` and the `SequenceExpression` family.

use crate::error::ModelError;
use crate::model::identifier::RefgetAccession;
use crate::model::meta::{Meta, impl_entity};
use crate::model::primitives::{IntOrRange, MoleculeType, ResidueAlphabet, SequenceString};

/// A VRS `SequenceReference`: a sequence identified by its RefGet accession, optionally
/// annotated with alphabet, molecule type, circularity and the literal sequence.
///
/// Only `refgetAccession` and `type` contribute to digests; `id` (e.g. `NC_000001.11`) and the
/// other annotations are decorative.
#[derive(Clone, Debug, PartialEq)]
pub struct SequenceReference {
    refget_accession: RefgetAccession,
    residue_alphabet: Option<ResidueAlphabet>,
    molecule_type: Option<MoleculeType>,
    circular: Option<bool>,
    sequence: Option<SequenceString>,
    meta: Option<Box<Meta>>,
}

impl SequenceReference {
    /// A reference to the sequence with the given RefGet accession.
    pub fn new(refget_accession: RefgetAccession) -> Self {
        Self {
            refget_accession,
            residue_alphabet: None,
            molecule_type: None,
            circular: None,
            sequence: None,
            meta: None,
        }
    }

    /// Parse a `SQ.…` accession and build a reference.
    ///
    /// # Errors
    /// [`crate::error::IdentifierError::InvalidRefgetAccession`] for malformed input.
    pub fn parse(refget_accession: &str) -> Result<Self, crate::error::IdentifierError> {
        RefgetAccession::parse(refget_accession).map(Self::new)
    }

    /// Set `residueAlphabet`.
    #[must_use]
    pub fn with_residue_alphabet(mut self, alphabet: ResidueAlphabet) -> Self {
        self.residue_alphabet = Some(alphabet);
        self
    }

    /// Set `moleculeType`.
    #[must_use]
    pub fn with_molecule_type(mut self, molecule_type: MoleculeType) -> Self {
        self.molecule_type = Some(molecule_type);
        self
    }

    /// Set `circular`.
    #[must_use]
    pub fn with_circular(mut self, circular: bool) -> Self {
        self.circular = Some(circular);
        self
    }

    /// Set the literal `sequence`.
    #[must_use]
    pub fn with_sequence(mut self, sequence: SequenceString) -> Self {
        self.sequence = Some(sequence);
        self
    }

    /// The RefGet accession.
    #[inline]
    pub fn refget_accession(&self) -> &RefgetAccession {
        &self.refget_accession
    }

    /// `residueAlphabet`.
    #[inline]
    pub fn residue_alphabet(&self) -> Option<ResidueAlphabet> {
        self.residue_alphabet
    }

    /// `moleculeType`.
    #[inline]
    pub fn molecule_type(&self) -> Option<MoleculeType> {
        self.molecule_type
    }

    /// `circular`.
    #[inline]
    pub fn circular(&self) -> Option<bool> {
        self.circular
    }

    /// `true` only when the reference is explicitly flagged circular.
    #[inline]
    pub fn is_circular(&self) -> bool {
        self.circular == Some(true)
    }

    /// The literal sequence, if carried.
    #[inline]
    pub fn sequence(&self) -> Option<&SequenceString> {
        self.sequence.as_ref()
    }

    pub(crate) fn from_parts(
        refget_accession: RefgetAccession,
        residue_alphabet: Option<ResidueAlphabet>,
        molecule_type: Option<MoleculeType>,
        circular: Option<bool>,
        sequence: Option<SequenceString>,
        meta: Option<Box<Meta>>,
    ) -> Self {
        Self {
            refget_accession,
            residue_alphabet,
            molecule_type,
            circular,
            sequence,
            meta,
        }
    }
}
impl_entity!(SequenceReference);

/// A VRS `LiteralSequenceExpression`: an explicit sequence.
#[derive(Clone, Debug, PartialEq)]
pub struct LiteralSequenceExpression {
    sequence: SequenceString,
    meta: Option<Box<Meta>>,
}

impl LiteralSequenceExpression {
    /// Wrap a sequence.
    pub fn new(sequence: SequenceString) -> Self {
        Self {
            sequence,
            meta: None,
        }
    }

    /// Validate a residue string and wrap it.
    ///
    /// # Errors
    /// [`crate::error::SequenceStringError`] for invalid residues.
    pub fn parse(sequence: &str) -> Result<Self, crate::error::SequenceStringError> {
        SequenceString::new(sequence).map(Self::new)
    }

    /// The sequence.
    #[inline]
    pub fn sequence(&self) -> &SequenceString {
        &self.sequence
    }

    /// Consume, returning the sequence.
    pub fn into_sequence(self) -> SequenceString {
        self.sequence
    }

    pub(crate) fn from_parts(sequence: SequenceString, meta: Option<Box<Meta>>) -> Self {
        Self { sequence, meta }
    }
}
impl_entity!(LiteralSequenceExpression);

/// A VRS `ReferenceLengthExpression`: a sequence derived by circularly repeating a
/// `repeatSubunitLength`-residue subunit of the reference at the allele's location until
/// `length` residues are produced.
#[derive(Clone, Debug, PartialEq)]
pub struct ReferenceLengthExpression {
    length: IntOrRange,
    repeat_subunit_length: i64,
    sequence: Option<SequenceString>,
    meta: Option<Box<Meta>>,
}

impl ReferenceLengthExpression {
    /// Create an expression of `length` residues built from a `repeat_subunit_length` subunit.
    ///
    /// # Errors
    /// [`ModelError::NegativeCount`] if either count is negative.
    pub fn new(
        length: impl Into<IntOrRange>,
        repeat_subunit_length: i64,
    ) -> Result<Self, ModelError> {
        let length = length.into();
        if length.has_negative() {
            return Err(ModelError::NegativeCount {
                property: "length",
                value: length.lower_bound().or(length.upper_bound()).unwrap_or(-1),
            });
        }
        if repeat_subunit_length < 0 {
            return Err(ModelError::NegativeCount {
                property: "repeatSubunitLength",
                value: repeat_subunit_length,
            });
        }
        Ok(Self {
            length,
            repeat_subunit_length,
            sequence: None,
            meta: None,
        })
    }

    /// Attach the literal sequence encoded by the expression (decorative).
    #[must_use]
    pub fn with_sequence(mut self, sequence: SequenceString) -> Self {
        self.sequence = Some(sequence);
        self
    }

    /// `length`.
    #[inline]
    pub fn length(&self) -> IntOrRange {
        self.length
    }

    /// `repeatSubunitLength`.
    #[inline]
    pub fn repeat_subunit_length(&self) -> i64 {
        self.repeat_subunit_length
    }

    /// The literal sequence, if carried.
    #[inline]
    pub fn sequence(&self) -> Option<&SequenceString> {
        self.sequence.as_ref()
    }

    pub(crate) fn from_parts(
        length: IntOrRange,
        repeat_subunit_length: i64,
        sequence: Option<SequenceString>,
        meta: Option<Box<Meta>>,
    ) -> Result<Self, ModelError> {
        Self::new(length, repeat_subunit_length).map(|mut e| {
            e.sequence = sequence;
            e.meta = meta;
            e
        })
    }
}
impl_entity!(ReferenceLengthExpression);

/// A VRS `LengthExpression` (draft): a sequence known only by its length.
#[derive(Clone, Debug, PartialEq)]
pub struct LengthExpression {
    length: Option<IntOrRange>,
    meta: Option<Box<Meta>>,
}

impl LengthExpression {
    /// A sequence of the given length.
    ///
    /// # Errors
    /// [`ModelError::NegativeCount`] if the length is negative.
    pub fn new(length: impl Into<IntOrRange>) -> Result<Self, ModelError> {
        Self::from_parts(Some(length.into()), None)
    }

    /// A sequence of unknown length.
    pub fn unknown() -> Self {
        Self {
            length: None,
            meta: None,
        }
    }

    /// `length`.
    #[inline]
    pub fn length(&self) -> Option<IntOrRange> {
        self.length
    }

    pub(crate) fn from_parts(
        length: Option<IntOrRange>,
        meta: Option<Box<Meta>>,
    ) -> Result<Self, ModelError> {
        if let Some(l) = length
            && l.has_negative()
        {
            return Err(ModelError::NegativeCount {
                property: "length",
                value: l.lower_bound().or(l.upper_bound()).unwrap_or(-1),
            });
        }
        Ok(Self { length, meta })
    }
}
impl_entity!(LengthExpression);

/// A VRS `SequenceExpression`: one of the three ways to express a sequence state.
#[derive(Clone, Debug, PartialEq)]
pub enum SequenceExpression {
    /// An explicit sequence.
    Literal(LiteralSequenceExpression),
    /// A sequence derived from the reference by repetition.
    ReferenceLength(ReferenceLengthExpression),
    /// A sequence known only by length (draft).
    Length(LengthExpression),
}

impl SequenceExpression {
    /// Shorthand for a literal sequence expression.
    ///
    /// # Errors
    /// [`crate::error::SequenceStringError`] for invalid residues.
    pub fn literal(sequence: &str) -> Result<Self, crate::error::SequenceStringError> {
        LiteralSequenceExpression::parse(sequence).map(Self::Literal)
    }

    /// The literal expression, if this is one.
    #[inline]
    pub fn as_literal(&self) -> Option<&LiteralSequenceExpression> {
        match self {
            Self::Literal(l) => Some(l),
            _ => None,
        }
    }

    /// The reference-length expression, if this is one.
    #[inline]
    pub fn as_reference_length(&self) -> Option<&ReferenceLengthExpression> {
        match self {
            Self::ReferenceLength(r) => Some(r),
            _ => None,
        }
    }

    /// The length expression, if this is one.
    #[inline]
    pub fn as_length(&self) -> Option<&LengthExpression> {
        match self {
            Self::Length(l) => Some(l),
            _ => None,
        }
    }

    /// The VRS class name of the concrete expression.
    pub const fn type_name(&self) -> &'static str {
        match self {
            Self::Literal(_) => "LiteralSequenceExpression",
            Self::ReferenceLength(_) => "ReferenceLengthExpression",
            Self::Length(_) => "LengthExpression",
        }
    }

    /// The literal sequence, when one is available (always for literal expressions, when
    /// carried for reference-length expressions, never for length expressions).
    pub fn sequence(&self) -> Option<&SequenceString> {
        match self {
            Self::Literal(l) => Some(l.sequence()),
            Self::ReferenceLength(r) => r.sequence(),
            Self::Length(_) => None,
        }
    }
}

impl From<LiteralSequenceExpression> for SequenceExpression {
    fn from(v: LiteralSequenceExpression) -> Self {
        Self::Literal(v)
    }
}
impl From<ReferenceLengthExpression> for SequenceExpression {
    fn from(v: ReferenceLengthExpression) -> Self {
        Self::ReferenceLength(v)
    }
}
impl From<LengthExpression> for SequenceExpression {
    fn from(v: LengthExpression) -> Self {
        Self::Length(v)
    }
}
impl From<SequenceString> for SequenceExpression {
    fn from(s: SequenceString) -> Self {
        Self::Literal(LiteralSequenceExpression::new(s))
    }
}
