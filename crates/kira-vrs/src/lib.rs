//! # kira-vrs
//!
//! A native Rust implementation of the GA4GH **Variation Representation Specification**
//! (VRS) 2.1: a strongly typed domain model, standards-compatible JSON serialization,
//! RFC 8785 digest serialization with `sha512t24u` computed identifiers, and the VRS
//! normalization algorithms.
//!
//! ```text
//! genomic variation ─▶ typed VRS model ─▶ canonical (normalized) form ─▶ stable identifier
//! ```
//!
//! The implemented specification revision is pinned in [`spec`].
//!
//! ## Quick start
//!
//! ```
//! use kira_vrs::prelude::*;
//!
//! # fn main() -> Result<(), kira_vrs::Error> {
//! // NC_000019.10:g.44908822C>T (rs7412), on the RefGet accession of GRCh38 chr19.
//! let chr19 = SequenceReference::parse("SQ.IIB53T8CNeJJdUqzn9V_JnRtQadwWCbl")?;
//! let location = SequenceLocation::new(chr19, 44_908_821, 44_908_822)?;
//! let allele = Allele::new(location, SequenceExpression::literal("T")?);
//!
//! assert_eq!(allele.identifier().to_string(), "ga4gh:VA.0AePZIWZUNsUlQTamyLrjm2HWUw2opLt");
//!
//! let json = kira_vrs::json::to_string(&allele)?;
//! let back: Allele = kira_vrs::json::from_str(&json)?;
//! assert_eq!(back, allele);
//! # Ok(()) }
//! ```
//!
//! ## Modules
//!
//! * [`model`] — the VRS classes as Rust types with enforced invariants.
//! * [`digest`] — digest serialization, `sha512t24u`, computed identifiers.
//! * [`json`] — JSON interchange conforming to the official VRS JSON Schema.
//! * [`normalize`] — allele normalization (fully-justified, reference-length encoding),
//!   plus the cis-phased-block, adjacency and relative-allele conventions.
//! * [`spec`] — the pinned upstream specification revision and maturity annotations.

#![forbid(unsafe_code)]
#![cfg_attr(docsrs, feature(doc_cfg))]

pub mod digest;
pub mod error;
pub mod json;
pub mod model;
pub mod normalize;
pub mod spec;

pub use error::{Error, Result};

/// Commonly used items.
pub mod prelude {
    pub use crate::digest::{DigestSerialize, Identifiable};
    pub use crate::model::{
        Adjacency, Allele, AnchorOrientation, CisPhasedBlock, CopyChange, CopyNumberChange,
        CopyNumberCount, DerivativeComponent, DerivativeMolecule, Digest, Entity, Expression,
        Extension, IntOrRange, Iri, IriOr, LengthExpression, LiteralSequenceExpression, Location,
        Meta, MolecularVariation, MoleculeType, Orientation, Range, ReferenceLengthExpression,
        RefgetAccession, RelativeAllele, RelativeSequenceLocation, ResidueAlphabet,
        SequenceExpression, SequenceLocation, SequenceOffsetLocation, SequenceReference,
        SequenceString, Syntax, SystemicVariation, Terminus, TraversalBlock, TypePrefix, Variation,
        VrsIdentifier,
    };
    pub use crate::normalize::{NormalizeOptions, SequenceProvider};
}
