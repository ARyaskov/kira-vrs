//! GA4GH computed identifiers: digest serialization, `sha512t24u` and CURIE construction.
//!
//! The pipeline defined by the VRS *Computed Identifiers* convention is
//!
//! ```text
//! typed object ──▶ digest serialization (RFC 8785 canonical JSON of the
//!                  inherent properties, nested identifiable objects replaced
//!                  by their digests, unordered arrays sorted)
//!              ──▶ SHA-512, truncated to 24 bytes, base64url  ("sha512t24u")
//!              ──▶ "ga4gh:" + type prefix + "." + digest
//! ```
//!
//! # Implementation notes
//!
//! * Serialization is written directly into a byte buffer: the set and order of inherent
//!   keys is fixed per class and known at compile time, so no dynamic JSON object is built and
//!   no runtime key sorting takes place. The only sorting is of digest arrays flagged
//!   `ordered: false` in the schema (`CisPhasedBlock.members`).
//! * Nested identifiable objects are serialized into the *tail of the same buffer*, hashed,
//!   then truncated away and replaced by the digest string, so computing an identifier for an
//!   allele allocates a single buffer.
//! * Inherent properties whose value is absent are serialized as `null`. This follows the
//!   normative validation vectors (`{"adjoinedSequences":[…],"linker":null,"type":"Adjacency"}`)
//!   and the reference implementation, even though the prose of the convention says null
//!   fields are filtered; see `docs/design.md` for the analysis.
//! * Integers are written in exact decimal form. RFC 8785 mandates ECMAScript number
//!   formatting, which coincides with exact decimal for every integer of magnitude below
//!   2^53; genomic coordinates never approach that, and the reference implementation also
//!   writes exact integers.

pub(crate) mod jcs;
pub mod legacy;
mod serialize;

use sha2::{Digest as _, Sha512};

pub use crate::model::identifier::{Digest, TypePrefix, VrsIdentifier};

/// The GA4GH truncated digest: SHA-512, left-most 24 bytes, base64url without padding.
///
/// ```
/// # use kira_vrs::digest::sha512t24u;
/// assert_eq!(sha512t24u(b"ACGT").as_str(), "aKF498dAxcJAqme6QYQ7EZ07-fiw8Kw2");
/// ```
pub fn sha512t24u(data: &[u8]) -> Digest {
    let hash = Sha512::digest(data);
    let mut raw = [0u8; 24];
    raw.copy_from_slice(&hash[..24]);
    Digest::from_raw(&raw)
}

/// Objects with a GA4GH digest serialization (every VRS class with `ga4gh.inherent`).
///
/// Implemented for identifiable classes (which are additionally [`Identifiable`]) and for
/// value classes such as `SequenceReference` and `LiteralSequenceExpression`, which are
/// serialized inline inside their parents.
pub trait DigestSerialize {
    /// Append the digest serialization of `self` to `out`.
    fn write_digest_serialization(&self, out: &mut Vec<u8>);

    /// Append the representation of `self` when nested inside another object: the digest
    /// string for identifiable objects, the inline serialization otherwise.
    fn write_nested(&self, out: &mut Vec<u8>) {
        self.write_digest_serialization(out);
    }

    /// The digest serialization as owned bytes (UTF-8 canonical JSON).
    fn digest_serialization(&self) -> Vec<u8> {
        let mut out = Vec::with_capacity(192);
        self.write_digest_serialization(&mut out);
        out
    }
}

/// GA4GH identifiable objects: classes with a type prefix, for which a computed identifier
/// can be generated.
pub trait Identifiable: DigestSerialize {
    /// The type prefix used in the identifier (`VA` for alleles, ...).
    fn type_prefix(&self) -> TypePrefix;

    /// The sha512t24u digest of the digest serialization.
    fn digest(&self) -> Digest {
        let mut out = Vec::with_capacity(192);
        self.digest_with(&mut out)
    }

    /// The digest, using `scratch` as the serialization buffer (cleared first, left holding
    /// the serialization). Reusing one buffer across a cohort avoids an allocation per
    /// object; see the `digest/identifier_reused_buffer` benchmark.
    fn digest_with(&self, scratch: &mut Vec<u8>) -> Digest {
        scratch.clear();
        self.write_digest_serialization(scratch);
        sha512t24u(scratch)
    }

    /// The GA4GH computed identifier, e.g. `ga4gh:VA.0AePZIWZUNsUlQTamyLrjm2HWUw2opLt`.
    fn identifier(&self) -> VrsIdentifier {
        VrsIdentifier::new(self.type_prefix(), self.digest())
    }

    /// The identifier, using `scratch` as the serialization buffer (see [`digest_with`]).
    ///
    /// [`digest_with`]: Identifiable::digest_with
    fn identifier_with(&self, scratch: &mut Vec<u8>) -> VrsIdentifier {
        VrsIdentifier::new(self.type_prefix(), self.digest_with(scratch))
    }
}

/// Serialize `obj` into the tail of `out`, hash it, drop the tail and write the quoted digest.
///
/// This is how nested identifiable objects appear in their parent's serialization, without
/// allocating a second buffer.
pub(crate) fn write_nested_digest<T: DigestSerialize + ?Sized>(out: &mut Vec<u8>, obj: &T) {
    let mark = out.len();
    obj.write_digest_serialization(out);
    let digest = sha512t24u(&out[mark..]);
    out.truncate(mark);
    jcs::write_quoted_ascii(out, digest.as_bytes());
}

/// Compute the digest of a nested object using `out` as scratch space (restored on return).
pub(crate) fn nested_digest<T: DigestSerialize + ?Sized>(out: &mut Vec<u8>, obj: &T) -> Digest {
    let mark = out.len();
    obj.write_digest_serialization(out);
    let digest = sha512t24u(&out[mark..]);
    out.truncate(mark);
    digest
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::model::*;

    fn rs7412() -> Allele {
        let reference = SequenceReference::parse("SQ.IIB53T8CNeJJdUqzn9V_JnRtQadwWCbl").unwrap();
        let location = SequenceLocation::new(reference, 44_908_821, 44_908_822).unwrap();
        Allele::new(location, SequenceExpression::literal("T").unwrap())
    }

    #[test]
    fn spec_example_allele() {
        let allele = rs7412();
        let ser = String::from_utf8(allele.digest_serialization()).unwrap();
        assert_eq!(
            ser,
            r#"{"location":"wIlaGykfwHIpPY2Fcxtbx4TINbbODFVz","state":{"sequence":"T","type":"LiteralSequenceExpression"},"type":"Allele"}"#
        );
        assert_eq!(
            allele.identifier().to_string(),
            "ga4gh:VA.0AePZIWZUNsUlQTamyLrjm2HWUw2opLt"
        );
        assert_eq!(
            allele.sequence_location().unwrap().identifier().to_string(),
            "ga4gh:SL.wIlaGykfwHIpPY2Fcxtbx4TINbbODFVz"
        );
    }

    #[test]
    fn decorative_metadata_does_not_change_digest() {
        let plain = rs7412();
        let decorated = rs7412()
            .with_id("my-id")
            .with_name("rs7412")
            .with_expression(Expression::new(Syntax::HgvsG, "NC_000019.10:g.44908822C>T"));
        assert_eq!(plain.digest(), decorated.digest());
    }

    #[test]
    fn iri_location_uses_digest_component() {
        let iri = Iri::new("ga4gh:SL.wIlaGykfwHIpPY2Fcxtbx4TINbbODFVz");
        let allele = Allele::new(iri, SequenceExpression::literal("T").unwrap());
        assert_eq!(
            allele.identifier().to_string(),
            "ga4gh:VA.0AePZIWZUNsUlQTamyLrjm2HWUw2opLt"
        );
    }
}
