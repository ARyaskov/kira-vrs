//! Legacy VRS 1.3 digests and identifiers (implementation-specific extension).
//!
//! VRS 2 identifiers are not stable across major versions: every VRS 1.x identifier differs
//! from its VRS 2.x counterpart. Databases annotated with 1.3 identifiers need a bridge, so
//! this module reproduces the 1.3 serialization for the two classes that had equivalents,
//! `SequenceLocation` (1.3 prefix `VSL`) and `Allele` (`VA`), exactly as the reference
//! implementation does for its `as_version="1.3"` option. The upstream validation vectors
//! (`ga4gh_1_3_serialize` / `ga4gh_1_3_identify`) cover both.
//!
//! This is **not** part of VRS 2.1; it is an interoperability aid.

use crate::digest::sha512t24u;
use crate::error::LegacyDigestError;
use crate::model::{Allele, Digest, IntOrRange, IriOr, SequenceExpression, SequenceLocation};

use super::jcs::{write_i64, write_quoted_ascii};

/// A VRS 1.3 computed identifier (`ga4gh:VSL.…` / `ga4gh:VA.…`).
#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub struct LegacyIdentifier {
    prefix: &'static str,
    digest: Digest,
}

impl LegacyIdentifier {
    /// The 1.3 type prefix (`VSL` or `VA`).
    pub fn prefix(&self) -> &'static str {
        self.prefix
    }

    /// The digest.
    pub fn digest(&self) -> Digest {
        self.digest
    }
}

impl std::fmt::Display for LegacyIdentifier {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "ga4gh:{}.{}", self.prefix, self.digest)
    }
}

fn write_1_3_coordinate(out: &mut Vec<u8>, v: IntOrRange) {
    match v {
        IntOrRange::Int(i) => {
            out.extend_from_slice(br#"{"type":"Number","value":"#);
            write_i64(out, i);
            out.push(b'}');
        }
        IntOrRange::Range(r) => match (r.min(), r.max()) {
            (Some(min), Some(max)) => {
                out.extend_from_slice(br#"{"max":"#);
                write_i64(out, max);
                out.extend_from_slice(br#","min":"#);
                write_i64(out, min);
                out.extend_from_slice(br#","type":"DefiniteRange"}"#);
            }
            (Some(min), None) => {
                out.extend_from_slice(br#"{"comparator":">=","type":"IndefiniteRange","value":"#);
                write_i64(out, min);
                out.push(b'}');
            }
            (None, Some(max)) => {
                out.extend_from_slice(br#"{"comparator":"<=","type":"IndefiniteRange","value":"#);
                write_i64(out, max);
                out.push(b'}');
            }
            (None, None) => unreachable!("Range invariant: at least one bound"),
        },
    }
}

/// VRS 1.3 digest serialization of a sequence location.
///
/// # Errors
/// The location must have an inline `SequenceReference` and both coordinates.
pub fn sequence_location_1_3_serialization(
    location: &SequenceLocation,
) -> Result<Vec<u8>, LegacyDigestError> {
    let accession = location
        .refget_accession()
        .ok_or(LegacyDigestError::MissingSequenceReference)?;
    let (Some(start), Some(end)) = (location.start(), location.end()) else {
        return Err(LegacyDigestError::MissingCoordinates);
    };
    let mut out = Vec::with_capacity(192);
    out.extend_from_slice(br#"{"interval":{"end":"#);
    write_1_3_coordinate(&mut out, end);
    out.extend_from_slice(br#","start":"#);
    write_1_3_coordinate(&mut out, start);
    out.extend_from_slice(br#","type":"SequenceInterval"},"sequence_id":"#);
    write_quoted_ascii(&mut out, accession.digest_str().as_bytes());
    out.extend_from_slice(br#","type":"SequenceLocation"}"#);
    Ok(out)
}

/// VRS 1.3 digest of a sequence location.
///
/// # Errors
/// See [`sequence_location_1_3_serialization`].
pub fn sequence_location_1_3_digest(
    location: &SequenceLocation,
) -> Result<Digest, LegacyDigestError> {
    sequence_location_1_3_serialization(location).map(|b| sha512t24u(&b))
}

/// VRS 1.3 identifier (`ga4gh:VSL.…`) of a sequence location.
///
/// # Errors
/// See [`sequence_location_1_3_serialization`].
pub fn sequence_location_1_3_identifier(
    location: &SequenceLocation,
) -> Result<LegacyIdentifier, LegacyDigestError> {
    Ok(LegacyIdentifier {
        prefix: "VSL",
        digest: sequence_location_1_3_digest(location)?,
    })
}

/// VRS 1.3 digest serialization of an allele.
///
/// # Errors
/// The location must be an inline `SequenceLocation` acceptable to
/// [`sequence_location_1_3_serialization`]; the state must be a literal sequence expression or
/// a reference-length expression carrying its literal `sequence`.
pub fn allele_1_3_serialization(allele: &Allele) -> Result<Vec<u8>, LegacyDigestError> {
    let location = match allele.location() {
        IriOr::Object(l) => l,
        IriOr::Iri(_) => return Err(LegacyDigestError::UnsupportedLocation("IRI")),
    };
    let sequence = match allele.state() {
        SequenceExpression::Literal(l) => l.sequence(),
        SequenceExpression::ReferenceLength(r) => r
            .sequence()
            .ok_or(LegacyDigestError::MissingStateSequence)?,
        SequenceExpression::Length(_) => {
            return Err(LegacyDigestError::UnsupportedState("LengthExpression"));
        }
    };
    let location_digest = sequence_location_1_3_digest(location)?;
    let mut out = Vec::with_capacity(160);
    out.extend_from_slice(br#"{"location":"#);
    write_quoted_ascii(&mut out, location_digest.as_bytes());
    out.extend_from_slice(br#","state":{"sequence":"#);
    write_quoted_ascii(&mut out, sequence.as_bytes());
    out.extend_from_slice(br#","type":"LiteralSequenceExpression"},"type":"Allele"}"#);
    Ok(out)
}

/// VRS 1.3 digest of an allele.
///
/// # Errors
/// See [`allele_1_3_serialization`].
pub fn allele_1_3_digest(allele: &Allele) -> Result<Digest, LegacyDigestError> {
    allele_1_3_serialization(allele).map(|b| sha512t24u(&b))
}

/// VRS 1.3 identifier (`ga4gh:VA.…`) of an allele.
///
/// # Errors
/// See [`allele_1_3_serialization`].
pub fn allele_1_3_identifier(allele: &Allele) -> Result<LegacyIdentifier, LegacyDigestError> {
    Ok(LegacyIdentifier {
        prefix: "VA",
        digest: allele_1_3_digest(allele)?,
    })
}
