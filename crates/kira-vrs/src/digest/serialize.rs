//! Per-class digest serializers.
//!
//! Each implementation writes exactly the class's `ga4gh.inherent` properties, in Unicode
//! code point order of the property names, as required by RFC 8785. The key order is fixed at
//! compile time here and cross-checked against the upstream JSON Schemas by the
//! `kira-vrs-validation` crate.

use super::jcs::{
    write_i64, write_int_or_range, write_key, write_null, write_opt_int_or_range,
    write_quoted_ascii, write_str,
};
use super::{DigestSerialize, Identifiable, nested_digest, write_nested_digest};
use crate::model::identifier::{Digest, TypePrefix};
use crate::model::*;

/// Write an `IriOr<T>` that appears as a nested value: IRIs are written as their GA4GH digest
/// component when they are computed identifiers, otherwise verbatim; objects use
/// [`DigestSerialize::write_nested`].
fn write_iri_or<T: DigestSerialize>(out: &mut Vec<u8>, v: &IriOr<T>) {
    match v {
        IriOr::Iri(iri) => write_iri(out, iri),
        IriOr::Object(obj) => obj.write_nested(out),
    }
}

fn write_iri(out: &mut Vec<u8>, iri: &Iri) {
    match iri.ga4gh_digest() {
        Some(d) => write_quoted_ascii(out, d.as_bytes()),
        None => write_str(out, iri.as_str()),
    }
}

fn write_opt_iri_or<T: DigestSerialize>(out: &mut Vec<u8>, v: Option<&IriOr<T>>) {
    match v {
        None => write_null(out),
        Some(v) => write_iri_or(out, v),
    }
}

fn write_type(out: &mut Vec<u8>, name: &str) {
    write_key(out, "type");
    write_quoted_ascii(out, name.as_bytes());
    out.push(b'}');
}

// ---- value objects (inlined when nested) ---------------------------------------------------

impl DigestSerialize for SequenceReference {
    fn write_digest_serialization(&self, out: &mut Vec<u8>) {
        out.push(b'{');
        write_key(out, "refgetAccession");
        write_quoted_ascii(out, self.refget_accession().as_str().as_bytes());
        out.push(b',');
        write_type(out, "SequenceReference");
    }
}

impl DigestSerialize for LiteralSequenceExpression {
    fn write_digest_serialization(&self, out: &mut Vec<u8>) {
        out.push(b'{');
        write_key(out, "sequence");
        write_quoted_ascii(out, self.sequence().as_bytes());
        out.push(b',');
        write_type(out, "LiteralSequenceExpression");
    }
}

impl DigestSerialize for ReferenceLengthExpression {
    fn write_digest_serialization(&self, out: &mut Vec<u8>) {
        out.push(b'{');
        write_key(out, "length");
        write_int_or_range(out, self.length());
        out.push(b',');
        write_key(out, "repeatSubunitLength");
        write_i64(out, self.repeat_subunit_length());
        out.push(b',');
        write_type(out, "ReferenceLengthExpression");
    }
}

impl DigestSerialize for LengthExpression {
    fn write_digest_serialization(&self, out: &mut Vec<u8>) {
        out.push(b'{');
        write_key(out, "length");
        write_opt_int_or_range(out, self.length());
        out.push(b',');
        write_type(out, "LengthExpression");
    }
}

impl DigestSerialize for SequenceExpression {
    fn write_digest_serialization(&self, out: &mut Vec<u8>) {
        match self {
            Self::Literal(e) => e.write_digest_serialization(out),
            Self::ReferenceLength(e) => e.write_digest_serialization(out),
            Self::Length(e) => e.write_digest_serialization(out),
        }
    }
}

impl DigestSerialize for SequenceOffsetLocation {
    fn write_digest_serialization(&self, out: &mut Vec<u8>) {
        out.push(b'{');
        write_key(out, "anchor");
        write_i64(out, self.anchor());
        out.push(b',');
        write_key(out, "anchorOrientation");
        write_quoted_ascii(out, self.anchor_orientation().as_str().as_bytes());
        out.push(b',');
        write_key(out, "offsetEnd");
        write_opt_int_or_range(out, self.offset_end());
        out.push(b',');
        write_key(out, "offsetStart");
        write_opt_int_or_range(out, self.offset_start());
        out.push(b',');
        write_key(out, "sequenceReference");
        write_iri_or(out, self.sequence_reference());
        out.push(b',');
        write_type(out, "SequenceOffsetLocation");
    }
}

impl DigestSerialize for TraversalBlock {
    fn write_digest_serialization(&self, out: &mut Vec<u8>) {
        out.push(b'{');
        write_key(out, "component");
        match self.component() {
            Some(a) => write_nested_digest(out, a),
            None => write_null(out),
        }
        out.push(b',');
        write_key(out, "orientation");
        match self.orientation() {
            Some(o) => write_quoted_ascii(out, o.as_str().as_bytes()),
            None => write_null(out),
        }
        out.push(b',');
        write_type(out, "TraversalBlock");
    }
}

// ---- identifiable objects -------------------------------------------------------------------

macro_rules! identifiable {
    ($ty:ty, $prefix:expr) => {
        impl Identifiable for $ty {
            #[inline]
            fn type_prefix(&self) -> TypePrefix {
                $prefix
            }
        }
    };
}

impl DigestSerialize for SequenceLocation {
    fn write_digest_serialization(&self, out: &mut Vec<u8>) {
        out.push(b'{');
        write_key(out, "end");
        write_opt_int_or_range(out, self.end());
        out.push(b',');
        write_key(out, "sequenceReference");
        write_opt_iri_or(out, self.sequence_reference());
        out.push(b',');
        write_key(out, "start");
        write_opt_int_or_range(out, self.start());
        out.push(b',');
        write_type(out, "SequenceLocation");
    }
    fn write_nested(&self, out: &mut Vec<u8>) {
        write_nested_digest(out, self);
    }
}
identifiable!(SequenceLocation, TypePrefix::SequenceLocation);

impl DigestSerialize for RelativeSequenceLocation {
    fn write_digest_serialization(&self, out: &mut Vec<u8>) {
        out.push(b'{');
        write_key(out, "baseSequenceLocation");
        write_iri_or(out, self.base_sequence_location());
        out.push(b',');
        write_key(out, "mappedSequenceLocation");
        write_iri_or(out, self.mapped_sequence_location());
        out.push(b',');
        write_type(out, "RelativeSequenceLocation");
    }
    fn write_nested(&self, out: &mut Vec<u8>) {
        write_nested_digest(out, self);
    }
}
identifiable!(
    RelativeSequenceLocation,
    TypePrefix::RelativeSequenceLocation
);

impl DigestSerialize for Location {
    fn write_digest_serialization(&self, out: &mut Vec<u8>) {
        match self {
            Self::Sequence(l) => l.write_digest_serialization(out),
            Self::RelativeSequence(l) => l.write_digest_serialization(out),
        }
    }
    fn write_nested(&self, out: &mut Vec<u8>) {
        write_nested_digest(out, self);
    }
}
impl Identifiable for Location {
    fn type_prefix(&self) -> TypePrefix {
        match self {
            Self::Sequence(l) => l.type_prefix(),
            Self::RelativeSequence(l) => l.type_prefix(),
        }
    }
}

impl DigestSerialize for Allele {
    fn write_digest_serialization(&self, out: &mut Vec<u8>) {
        out.push(b'{');
        write_key(out, "location");
        write_iri_or(out, self.location());
        out.push(b',');
        write_key(out, "state");
        self.state().write_digest_serialization(out);
        out.push(b',');
        write_type(out, "Allele");
    }
    fn write_nested(&self, out: &mut Vec<u8>) {
        write_nested_digest(out, self);
    }
}
identifiable!(Allele, TypePrefix::Allele);

impl DigestSerialize for RelativeAllele {
    fn write_digest_serialization(&self, out: &mut Vec<u8>) {
        out.push(b'{');
        write_key(out, "baseState");
        self.base_state().write_digest_serialization(out);
        out.push(b',');
        write_key(out, "mappedState");
        self.mapped_state().write_digest_serialization(out);
        out.push(b',');
        write_key(out, "relativeLocation");
        write_iri_or(out, self.relative_location());
        out.push(b',');
        write_type(out, "RelativeAllele");
    }
    fn write_nested(&self, out: &mut Vec<u8>) {
        write_nested_digest(out, self);
    }
}
identifiable!(RelativeAllele, TypePrefix::RelativeAllele);

/// Digest strings of an `IriOr<T>` array member, for sorting unordered member lists.
///
/// Members that are non-identifier IRIs have no digest; they sort by their raw text, which is
/// what the reference implementation does with `sorted()` over the serialized strings.
enum MemberKey<'a> {
    Digest(Digest),
    Text(&'a str),
}

impl MemberKey<'_> {
    fn as_bytes(&self) -> &[u8] {
        match self {
            Self::Digest(d) => d.as_bytes(),
            Self::Text(s) => s.as_bytes(),
        }
    }
}

impl DigestSerialize for CisPhasedBlock {
    fn write_digest_serialization(&self, out: &mut Vec<u8>) {
        out.push(b'{');
        write_key(out, "members");
        // `members` is `ordered: false`: sort the member digests by code point.
        let mut keys: Vec<MemberKey<'_>> = self
            .members()
            .iter()
            .map(|m| match m {
                IriOr::Object(a) => MemberKey::Digest(nested_digest(out, a)),
                IriOr::Iri(iri) => match iri.ga4gh_digest() {
                    Some(d) => MemberKey::Text(d),
                    None => MemberKey::Text(iri.as_str()),
                },
            })
            .collect();
        keys.sort_unstable_by(|a, b| a.as_bytes().cmp(b.as_bytes()));
        out.push(b'[');
        for (i, k) in keys.iter().enumerate() {
            if i > 0 {
                out.push(b',');
            }
            match k {
                MemberKey::Digest(d) => write_quoted_ascii(out, d.as_bytes()),
                MemberKey::Text(s) => write_str(out, s),
            }
        }
        out.push(b']');
        out.push(b',');
        write_type(out, "CisPhasedBlock");
    }
    fn write_nested(&self, out: &mut Vec<u8>) {
        write_nested_digest(out, self);
    }
}
identifiable!(CisPhasedBlock, TypePrefix::CisPhasedBlock);

impl DigestSerialize for Adjacency {
    fn write_digest_serialization(&self, out: &mut Vec<u8>) {
        out.push(b'{');
        write_key(out, "adjoinedSequences");
        out.push(b'[');
        let [a, b] = self.adjoined_sequences();
        write_iri_or(out, a);
        out.push(b',');
        write_iri_or(out, b);
        out.push(b']');
        out.push(b',');
        write_key(out, "linker");
        match self.linker() {
            Some(l) => l.write_digest_serialization(out),
            None => write_null(out),
        }
        out.push(b',');
        write_type(out, "Adjacency");
    }
    fn write_nested(&self, out: &mut Vec<u8>) {
        write_nested_digest(out, self);
    }
}
identifiable!(Adjacency, TypePrefix::Adjacency);

impl DigestSerialize for Terminus {
    fn write_digest_serialization(&self, out: &mut Vec<u8>) {
        out.push(b'{');
        write_key(out, "location");
        write_iri_or(out, self.location());
        out.push(b',');
        write_type(out, "Terminus");
    }
    fn write_nested(&self, out: &mut Vec<u8>) {
        write_nested_digest(out, self);
    }
}
identifiable!(Terminus, TypePrefix::Terminus);

impl DigestSerialize for DerivativeComponent {
    fn write_digest_serialization(&self, out: &mut Vec<u8>) {
        match self {
            Self::Iri(iri) => write_iri(out, iri),
            Self::Allele(a) => a.write_nested(out),
            Self::CisPhasedBlock(c) => c.write_nested(out),
            Self::Terminus(t) => t.write_nested(out),
            Self::TraversalBlock(t) => t.write_nested(out),
        }
    }
}

impl DigestSerialize for DerivativeMolecule {
    fn write_digest_serialization(&self, out: &mut Vec<u8>) {
        out.push(b'{');
        write_key(out, "components");
        out.push(b'[');
        for (i, c) in self.components().iter().enumerate() {
            if i > 0 {
                out.push(b',');
            }
            c.write_digest_serialization(out);
        }
        out.push(b']');
        out.push(b',');
        write_type(out, "DerivativeMolecule");
    }
    fn write_nested(&self, out: &mut Vec<u8>) {
        write_nested_digest(out, self);
    }
}
identifiable!(DerivativeMolecule, TypePrefix::DerivativeMolecule);

impl DigestSerialize for CopyNumberCount {
    fn write_digest_serialization(&self, out: &mut Vec<u8>) {
        out.push(b'{');
        write_key(out, "copies");
        write_int_or_range(out, self.copies());
        out.push(b',');
        write_key(out, "location");
        write_iri_or(out, self.location());
        out.push(b',');
        write_type(out, "CopyNumberCount");
    }
    fn write_nested(&self, out: &mut Vec<u8>) {
        write_nested_digest(out, self);
    }
}
identifiable!(CopyNumberCount, TypePrefix::CopyNumberCount);

impl DigestSerialize for CopyNumberChange {
    fn write_digest_serialization(&self, out: &mut Vec<u8>) {
        out.push(b'{');
        write_key(out, "copyChange");
        write_str(out, self.copy_change().as_str());
        out.push(b',');
        write_key(out, "location");
        write_iri_or(out, self.location());
        out.push(b',');
        write_type(out, "CopyNumberChange");
    }
    fn write_nested(&self, out: &mut Vec<u8>) {
        write_nested_digest(out, self);
    }
}
identifiable!(CopyNumberChange, TypePrefix::CopyNumberChange);

macro_rules! union_digest {
    ($union:ident: $( $variant:ident ),+ $(,)?) => {
        impl DigestSerialize for $union {
            fn write_digest_serialization(&self, out: &mut Vec<u8>) {
                match self {
                    $( Self::$variant(v) => v.write_digest_serialization(out), )+
                }
            }
            fn write_nested(&self, out: &mut Vec<u8>) {
                write_nested_digest(out, self);
            }
        }
        impl Identifiable for $union {
            fn type_prefix(&self) -> TypePrefix {
                match self {
                    $( Self::$variant(v) => v.type_prefix(), )+
                }
            }
        }
    };
}

union_digest!(MolecularVariation: Allele, RelativeAllele, CisPhasedBlock, Adjacency, Terminus, DerivativeMolecule);
union_digest!(SystemicVariation: CopyNumberCount, CopyNumberChange);
union_digest!(Variation: Allele, RelativeAllele, CisPhasedBlock, Adjacency, Terminus, DerivativeMolecule, CopyNumberCount, CopyNumberChange);
