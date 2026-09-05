//! The upstream specification revision implemented by this crate, and the maturity level of
//! each class, as machine-readable constants.
//!
//! The same values are recorded in `Cargo.toml` under `[package.metadata.vrs]` and in
//! `docs/spec-revision.md`. They are updated only by a deliberate re-sync
//! (`scripts/sync-upstream.sh`), never at build time.

/// VRS specification version implemented.
pub const VRS_VERSION: &str = "2.1.0";
/// Git tag of the VRS repository.
pub const VRS_TAG: &str = "2.1.0";
/// Exact commit of <https://github.com/ga4gh/vrs> the model, digests and normalization follow.
pub const VRS_REVISION: &str = "cf33bfa7618011087655d5a5898e518c9d96bcdb";
/// Date of that commit.
pub const VRS_REVISION_DATE: &str = "2026-09-01";
/// The JSON Schema `$id` base of the pinned release.
pub const VRS_SCHEMA_BASE: &str = "https://w3id.org/ga4gh/schema/vrs/2.1.0/json/";

/// gkm-core (GA4GH Genomic Knowledge Model core) version the VRS release depends on.
pub const GKM_CORE_VERSION: &str = "1.2.0";
/// Exact commit of <https://github.com/ga4gh/gkm-core>.
pub const GKM_CORE_REVISION: &str = "91abbb7d0f8f05a183303853c121abd76b8b765a";
/// The JSON Schema `$id` base of gkm-core.
pub const GKM_CORE_SCHEMA_BASE: &str = "https://w3id.org/ga4gh/schema/gkm-core/1.2.0/json/";

/// GA4GH maturity levels of specification features.
#[derive(Clone, Copy, PartialEq, Eq, Hash, Debug, PartialOrd, Ord)]
pub enum Maturity {
    /// Under active development; may change in any patch release.
    Draft,
    /// Stable enough for implementation; may change in a minor release.
    TrialUse,
    /// Normative.
    Normative,
}

impl Maturity {
    /// The label used in the schema (`draft`, `trial use`, `normative`).
    pub const fn as_str(self) -> &'static str {
        match self {
            Self::Draft => "draft",
            Self::TrialUse => "trial use",
            Self::Normative => "normative",
        }
    }
}

/// Maturity of every VRS 2.1.0 class, from the `maturity` annotation in the source schema.
pub const CLASS_MATURITY: &[(&str, Maturity)] = &[
    ("Allele", Maturity::TrialUse),
    ("RelativeAllele", Maturity::Draft),
    ("CisPhasedBlock", Maturity::TrialUse),
    ("Adjacency", Maturity::TrialUse),
    ("Terminus", Maturity::Draft),
    ("DerivativeMolecule", Maturity::Draft),
    ("TraversalBlock", Maturity::Draft),
    ("CopyNumberCount", Maturity::TrialUse),
    ("CopyNumberChange", Maturity::Draft),
    ("SequenceLocation", Maturity::TrialUse),
    ("RelativeSequenceLocation", Maturity::Draft),
    ("SequenceOffsetLocation", Maturity::Draft),
    ("SequenceReference", Maturity::TrialUse),
    ("LiteralSequenceExpression", Maturity::TrialUse),
    ("ReferenceLengthExpression", Maturity::TrialUse),
    ("LengthExpression", Maturity::Draft),
    ("Expression", Maturity::TrialUse),
    ("Range", Maturity::TrialUse),
    ("sequenceString", Maturity::TrialUse),
];

/// Look up the maturity of a class by its VRS name.
pub fn class_maturity(class: &str) -> Option<Maturity> {
    CLASS_MATURITY
        .iter()
        .find(|(c, _)| *c == class)
        .map(|(_, m)| *m)
}
