//! Identifier value types: sha512t24u digests, RefGet accessions, VRS type prefixes and GA4GH
//! computed identifiers (`ga4gh:VA.<digest>`).

use std::fmt;
use std::str::FromStr;

use crate::error::IdentifierError;

/// Length of a base64url-encoded sha512t24u digest.
pub const DIGEST_LEN: usize = 32;

/// `true` if `b` is in the base64url alphabet (`A-Z a-z 0-9 - _`).
#[inline]
pub(crate) const fn is_base64url(b: u8) -> bool {
    b.is_ascii_alphanumeric() || b == b'-' || b == b'_'
}

/// `true` if `s` is exactly 32 base64url characters.
pub(crate) fn is_digest_str(s: &str) -> bool {
    s.len() == DIGEST_LEN && s.bytes().all(is_base64url)
}

/// A sha512t24u truncated digest: SHA-512, truncated to 24 bytes, base64url-encoded
/// (32 ASCII characters, no padding).
///
/// Stored as its 32 encoded characters so that comparison (which VRS defines in terms of
/// Unicode code points of the encoded form) and serialization are memory copies.
#[derive(Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct Digest([u8; DIGEST_LEN]);

impl Digest {
    /// Parse a 32-character base64url digest.
    ///
    /// # Errors
    /// [`IdentifierError::InvalidDigest`] if the length or alphabet is wrong.
    pub fn parse(s: &str) -> Result<Self, IdentifierError> {
        if is_digest_str(s) {
            let mut buf = [0u8; DIGEST_LEN];
            buf.copy_from_slice(s.as_bytes());
            Ok(Self(buf))
        } else {
            Err(IdentifierError::InvalidDigest(s.to_owned()))
        }
    }

    /// Encode 24 raw digest bytes.
    pub fn from_raw(raw: &[u8; 24]) -> Self {
        Self(base64url_encode_24(raw))
    }

    /// The encoded digest.
    #[inline]
    pub fn as_str(&self) -> &str {
        std::str::from_utf8(&self.0).unwrap_or_default()
    }

    /// The encoded digest as bytes.
    #[inline]
    pub fn as_bytes(&self) -> &[u8; DIGEST_LEN] {
        &self.0
    }
}

impl fmt::Debug for Digest {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "Digest({})", self.as_str())
    }
}
impl fmt::Display for Digest {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(self.as_str())
    }
}
impl FromStr for Digest {
    type Err = IdentifierError;
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        Self::parse(s)
    }
}
impl AsRef<str> for Digest {
    fn as_ref(&self) -> &str {
        self.as_str()
    }
}

const B64URL: &[u8; 64] = b"ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789-_";

/// base64url (RFC 4648 §5) encoding of exactly 24 bytes → 32 characters, no padding.
pub(crate) fn base64url_encode_24(raw: &[u8; 24]) -> [u8; DIGEST_LEN] {
    let mut out = [0u8; DIGEST_LEN];
    let (triples, _) = raw.as_chunks::<3>();
    let (quads, _) = out.as_chunks_mut::<4>();
    for (chunk, dst) in triples.iter().zip(quads.iter_mut()) {
        let n = (u32::from(chunk[0]) << 16) | (u32::from(chunk[1]) << 8) | u32::from(chunk[2]);
        dst[0] = B64URL[(n >> 18) as usize & 63];
        dst[1] = B64URL[(n >> 12) as usize & 63];
        dst[2] = B64URL[(n >> 6) as usize & 63];
        dst[3] = B64URL[n as usize & 63];
    }
    out
}

/// Length of a RefGet accession string (`SQ.` + 32).
const REFGET_LEN: usize = 3 + DIGEST_LEN;

/// A GA4GH RefGet sequence accession: `SQ.` followed by the sha512t24u digest of the
/// upper-case sequence, e.g. `SQ.IIB53T8CNeJJdUqzn9V_JnRtQadwWCbl` (GRCh38 chr19).
///
/// This is the only kind of sequence identifier VRS permits in `SequenceReference`; conventional
/// accessions such as `NC_000019.10` must be translated by an external service (see
/// `docs/vcf-to-vrs.md`).
#[derive(Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct RefgetAccession([u8; REFGET_LEN]);

impl RefgetAccession {
    /// Parse `SQ.<32 base64url chars>`.
    ///
    /// # Errors
    /// [`IdentifierError::InvalidRefgetAccession`] otherwise.
    pub fn parse(s: &str) -> Result<Self, IdentifierError> {
        match s.strip_prefix("SQ.") {
            Some(d) if is_digest_str(d) => {
                let mut buf = [0u8; REFGET_LEN];
                buf.copy_from_slice(s.as_bytes());
                Ok(Self(buf))
            }
            _ => Err(IdentifierError::InvalidRefgetAccession(s.to_owned())),
        }
    }

    /// Build an accession from a digest.
    pub fn from_digest(digest: Digest) -> Self {
        let mut buf = [0u8; REFGET_LEN];
        buf[..3].copy_from_slice(b"SQ.");
        buf[3..].copy_from_slice(digest.as_bytes());
        Self(buf)
    }

    /// Compute the RefGet accession of a sequence.
    ///
    /// RefGet digests are defined over the *upper-case* sequence; the input is hashed as given,
    /// so callers must upper-case lowercase (soft-masked) FASTA input first.
    pub fn from_sequence(sequence: &[u8]) -> Self {
        Self::from_digest(crate::digest::sha512t24u(sequence))
    }

    /// The digest component (without `SQ.`).
    pub fn digest(&self) -> Digest {
        let mut buf = [0u8; DIGEST_LEN];
        buf.copy_from_slice(&self.0[3..]);
        Digest(buf)
    }

    /// The digest component as a string slice.
    #[inline]
    pub fn digest_str(&self) -> &str {
        &self.as_str()[3..]
    }

    /// The full accession text `SQ.…`.
    #[inline]
    pub fn as_str(&self) -> &str {
        std::str::from_utf8(&self.0).unwrap_or_default()
    }
}

impl fmt::Debug for RefgetAccession {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "RefgetAccession({})", self.as_str())
    }
}
impl fmt::Display for RefgetAccession {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(self.as_str())
    }
}
impl FromStr for RefgetAccession {
    type Err = IdentifierError;
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        Self::parse(s)
    }
}
impl TryFrom<&str> for RefgetAccession {
    type Error = IdentifierError;
    fn try_from(s: &str) -> Result<Self, Self::Error> {
        Self::parse(s)
    }
}
impl AsRef<str> for RefgetAccession {
    fn as_ref(&self) -> &str {
        self.as_str()
    }
}

/// GA4GH type prefixes used in computed identifiers (`ga4gh:<prefix>.<digest>`).
#[derive(Clone, Copy, PartialEq, Eq, Hash, Debug, PartialOrd, Ord)]
#[non_exhaustive]
pub enum TypePrefix {
    /// `VA` — Allele.
    Allele,
    /// `RA` — RelativeAllele (draft).
    RelativeAllele,
    /// `CPB` — CisPhasedBlock.
    CisPhasedBlock,
    /// `AJ` — Adjacency.
    Adjacency,
    /// `TM` — Terminus (draft).
    Terminus,
    /// `DM` — DerivativeMolecule (draft).
    DerivativeMolecule,
    /// `CN` — CopyNumberCount.
    CopyNumberCount,
    /// `CX` — CopyNumberChange (draft).
    CopyNumberChange,
    /// `SL` — SequenceLocation.
    SequenceLocation,
    /// `RSL` — RelativeSequenceLocation (draft).
    RelativeSequenceLocation,
    /// `SQ` — RefGet sequence.
    Sequence,
}

impl TypePrefix {
    /// The prefix text.
    pub const fn as_str(self) -> &'static str {
        match self {
            Self::Allele => "VA",
            Self::RelativeAllele => "RA",
            Self::CisPhasedBlock => "CPB",
            Self::Adjacency => "AJ",
            Self::Terminus => "TM",
            Self::DerivativeMolecule => "DM",
            Self::CopyNumberCount => "CN",
            Self::CopyNumberChange => "CX",
            Self::SequenceLocation => "SL",
            Self::RelativeSequenceLocation => "RSL",
            Self::Sequence => "SQ",
        }
    }

    /// The VRS class name the prefix identifies.
    pub const fn class_name(self) -> &'static str {
        match self {
            Self::Allele => "Allele",
            Self::RelativeAllele => "RelativeAllele",
            Self::CisPhasedBlock => "CisPhasedBlock",
            Self::Adjacency => "Adjacency",
            Self::Terminus => "Terminus",
            Self::DerivativeMolecule => "DerivativeMolecule",
            Self::CopyNumberCount => "CopyNumberCount",
            Self::CopyNumberChange => "CopyNumberChange",
            Self::SequenceLocation => "SequenceLocation",
            Self::RelativeSequenceLocation => "RelativeSequenceLocation",
            Self::Sequence => "Sequence",
        }
    }

    /// Parse a prefix.
    ///
    /// # Errors
    /// [`IdentifierError::UnknownTypePrefix`] for anything not in the VRS prefix table.
    pub fn parse(s: &str) -> Result<Self, IdentifierError> {
        Ok(match s {
            "VA" => Self::Allele,
            "RA" => Self::RelativeAllele,
            "CPB" => Self::CisPhasedBlock,
            "AJ" => Self::Adjacency,
            "TM" => Self::Terminus,
            "DM" => Self::DerivativeMolecule,
            "CN" => Self::CopyNumberCount,
            "CX" => Self::CopyNumberChange,
            "SL" => Self::SequenceLocation,
            "RSL" => Self::RelativeSequenceLocation,
            "SQ" => Self::Sequence,
            _ => return Err(IdentifierError::UnknownTypePrefix(s.to_owned())),
        })
    }
}

impl fmt::Display for TypePrefix {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(self.as_str())
    }
}
impl FromStr for TypePrefix {
    type Err = IdentifierError;
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        Self::parse(s)
    }
}

/// A GA4GH computed identifier, a W3C CURIE of the form `ga4gh:<type prefix>.<digest>`.
#[derive(Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct VrsIdentifier {
    prefix: TypePrefix,
    digest: Digest,
}

impl VrsIdentifier {
    /// Combine a type prefix and a digest.
    pub const fn new(prefix: TypePrefix, digest: Digest) -> Self {
        Self { prefix, digest }
    }

    /// Parse `ga4gh:<prefix>.<digest>`.
    ///
    /// # Errors
    /// [`IdentifierError::InvalidGa4ghIdentifier`] for malformed input;
    /// [`IdentifierError::UnknownTypePrefix`] for unknown prefixes.
    pub fn parse(s: &str) -> Result<Self, IdentifierError> {
        let bad = || IdentifierError::InvalidGa4ghIdentifier(s.to_owned());
        let rest = s.strip_prefix("ga4gh:").ok_or_else(bad)?;
        let (prefix, digest) = rest.split_once('.').ok_or_else(bad)?;
        let digest = Digest::parse(digest).map_err(|_| bad())?;
        Ok(Self {
            prefix: TypePrefix::parse(prefix)?,
            digest,
        })
    }

    /// The type prefix.
    #[inline]
    pub const fn prefix(&self) -> TypePrefix {
        self.prefix
    }

    /// The digest.
    #[inline]
    pub const fn digest(&self) -> Digest {
        self.digest
    }

    /// Render as a CURIE string (`ga4gh:VA.…`).
    pub fn to_curie(&self) -> String {
        self.to_string()
    }

    /// Render as an [`Iri`](crate::model::Iri) for use in `IriOr` fields.
    pub fn to_iri(&self) -> crate::model::Iri {
        crate::model::Iri::new(self.to_string())
    }
}

impl fmt::Debug for VrsIdentifier {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "VrsIdentifier({self})")
    }
}
impl fmt::Display for VrsIdentifier {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "ga4gh:{}.{}", self.prefix.as_str(), self.digest.as_str())
    }
}
impl FromStr for VrsIdentifier {
    type Err = IdentifierError;
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        Self::parse(s)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn base64url_matches_known_vectors() {
        // sha512("ACGT")[..24] base64url == aKF498dAxcJAqme6QYQ7EZ07-fiw8Kw2 (VRS functions.yaml)
        let d = crate::digest::sha512t24u(b"ACGT");
        assert_eq!(d.as_str(), "aKF498dAxcJAqme6QYQ7EZ07-fiw8Kw2");
        assert_eq!(
            crate::digest::sha512t24u(b"").as_str(),
            "z4PhNX7vuL3xVChQ1m2AB9Yg5AULVxXc"
        );
    }

    #[test]
    fn refget_parse() {
        let a = RefgetAccession::parse("SQ.IIB53T8CNeJJdUqzn9V_JnRtQadwWCbl").unwrap();
        assert_eq!(a.digest_str(), "IIB53T8CNeJJdUqzn9V_JnRtQadwWCbl");
        assert!(RefgetAccession::parse("ga4gh:SQ.IIB53T8CNeJJdUqzn9V_JnRtQadwWCbl").is_err());
        assert!(RefgetAccession::parse("SQ.short").is_err());
        assert_eq!(
            RefgetAccession::from_sequence(b"ACGT").as_str(),
            "SQ.aKF498dAxcJAqme6QYQ7EZ07-fiw8Kw2"
        );
    }

    #[test]
    fn identifier_round_trip() {
        let id = VrsIdentifier::parse("ga4gh:VA.0AePZIWZUNsUlQTamyLrjm2HWUw2opLt").unwrap();
        assert_eq!(id.prefix(), TypePrefix::Allele);
        assert_eq!(id.to_string(), "ga4gh:VA.0AePZIWZUNsUlQTamyLrjm2HWUw2opLt");
        assert!(VrsIdentifier::parse("ga4gh:ZZ.0AePZIWZUNsUlQTamyLrjm2HWUw2opLt").is_err());
        assert!(VrsIdentifier::parse("refseq:NC_000001.11").is_err());
    }
}
