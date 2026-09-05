//! VRS normalization.
//!
//! Normalization rewrites a variation into the canonical form the specification defines, so
//! that the same biological variant always yields the same computed identifier. It is a
//! *semantic* operation over reference sequence data, not a JSON transformation.
//!
//! | Class | Rule (VRS 2.1.0 *Normalization* convention) | Function |
//! |---|---|---|
//! | `Allele` | fully-justified (VOCA-derived) expansion with reference-length encoding | [`normalize_allele`] |
//! | `CisPhasedBlock` | members normalized, then ordered by digest | [`normalize_cis_phased_block`] |
//! | `Adjacency` | conventional orientation of the adjoined sequences | [`normalize_adjacency`] |
//! | `RelativeAllele` | base allele normalized; anchor selection helper | [`normalize_relative_allele`], [`preferred_anchor`] |
//! | everything else | returned unchanged (as the specification requires) | [`normalize`] |
//!
//! Reference sequence access is abstracted by [`SequenceProvider`]; [`InMemorySequenceProvider`]
//! is included for tests, examples and small pipelines.
//!
//! Every algorithm is documented in `docs/normalization.md` with its specification reference,
//! complexity and edge cases.

mod allele;
mod provider;
mod structural;

pub use allele::{expand_reference_length_expression, normalize_allele};
pub use provider::{InMemorySequenceProvider, SequenceProvider};
pub use structural::{
    AnchorChoice, normalize_adjacency, normalize_cis_phased_block, normalize_relative_allele,
    preferred_anchor,
};

use crate::error::NormalizeError;
use crate::model::Variation;

/// Options controlling normalization output.
#[derive(Clone, Debug, PartialEq, Eq)]
#[non_exhaustive]
pub struct NormalizeOptions {
    /// When an allele normalizes to a `ReferenceLengthExpression`, also attach the literal
    /// `sequence` if its length does not exceed this limit (`None` = always attach,
    /// `Some(0)` = never). The literal sequence is decorative and does not affect digests.
    /// Default: `Some(50)`, matching the reference implementation.
    pub rle_sequence_limit: Option<u64>,
}

impl Default for NormalizeOptions {
    fn default() -> Self {
        Self {
            rle_sequence_limit: Some(50),
        }
    }
}

impl NormalizeOptions {
    /// Defaults (attach literal RLE sequences up to 50 residues).
    pub fn new() -> Self {
        Self::default()
    }

    /// Never attach a literal sequence to reference-length expressions (smallest output).
    #[must_use]
    pub fn without_rle_sequence(mut self) -> Self {
        self.rle_sequence_limit = Some(0);
        self
    }

    /// Set the literal-sequence limit.
    #[must_use]
    pub fn with_rle_sequence_limit(mut self, limit: Option<u64>) -> Self {
        self.rle_sequence_limit = limit;
        self
    }

    pub(crate) fn attach_rle_sequence(&self, length: u64) -> bool {
        match self.rle_sequence_limit {
            None => true,
            Some(0) => false,
            Some(limit) => length <= limit,
        }
    }
}

/// Normalize any [`Variation`].
///
/// Classes without normalization rules are returned unchanged, as the specification requires.
///
/// # Errors
/// [`NormalizeError::Sequence`] if reference data cannot be fetched;
/// [`NormalizeError::Unsupported`] for circular sequence references.
pub fn normalize<P: SequenceProvider + ?Sized>(
    variation: &Variation,
    provider: &P,
    options: &NormalizeOptions,
) -> Result<Variation, NormalizeError> {
    Ok(match variation {
        Variation::Allele(a) => Variation::Allele(normalize_allele(a, provider, options)?),
        Variation::CisPhasedBlock(b) => {
            Variation::CisPhasedBlock(normalize_cis_phased_block(b, provider, options)?)
        }
        Variation::Adjacency(a) => Variation::from(normalize_adjacency(a)),
        Variation::RelativeAllele(r) => {
            Variation::from(normalize_relative_allele(r, provider, options)?)
        }
        other => other.clone(),
    })
}
