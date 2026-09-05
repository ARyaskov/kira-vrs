//! The VRS 2.1 information model as idiomatic Rust types.
//!
//! The model is a *domain* model, not a transcription of the JSON Schema:
//!
//! * invariants (non-negative coordinates, `start <= end`, minimum cardinalities, closed value
//!   sets) are enforced by constructors, so a value of these types is always schema-valid;
//! * polymorphic properties (`oneOf` in the schema) are Rust enums ([`SequenceExpression`],
//!   [`Location`], [`Variation`], [`IriOr`]);
//! * hot-path types are compact: an SNV [`Allele`] with an inline [`SequenceReference`] needs
//!   no heap allocation at all, because sequences up to 22 residues are stored inline,
//!   RefGet accessions are fixed-size arrays and decorative metadata lives behind a null
//!   pointer until it is needed.
//!
//! JSON is handled by the [`crate::json`] module; digests by [`crate::digest`]; normalization by
//! [`crate::normalize`].

pub(crate) mod identifier;
pub(crate) mod location;
pub(crate) mod meta;
pub(crate) mod primitives;
pub(crate) mod sequence;
pub(crate) mod variation;

pub use identifier::{DIGEST_LEN, Digest, RefgetAccession, TypePrefix, VrsIdentifier};
pub use location::{Location, RelativeSequenceLocation, SequenceLocation, SequenceOffsetLocation};
pub use meta::{Entity, Expression, Extension, Meta};
pub use primitives::{
    AnchorOrientation, CopyChange, IntOrRange, Iri, IriOr, MoleculeType, Orientation, Range,
    ResidueAlphabet, SequenceString, Syntax,
};
pub use sequence::{
    LengthExpression, LiteralSequenceExpression, ReferenceLengthExpression, SequenceExpression,
    SequenceReference,
};
pub use variation::{
    Adjacency, Allele, CisPhasedBlock, CopyNumberChange, CopyNumberCount, DerivativeComponent,
    DerivativeMolecule, MolecularVariation, RelativeAllele, SystemicVariation, Terminus,
    TraversalBlock, Variation,
};

primitives::iri_or_from_iri!(
    SequenceReference,
    SequenceLocation,
    SequenceOffsetLocation,
    RelativeSequenceLocation,
    Location,
    Allele,
);
