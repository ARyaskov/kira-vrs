//! `serde::Deserialize` implementations for the model.

use std::fmt;

use serde::Deserialize;
use serde::de::value::MapAccessDeserializer;
use serde::de::{self, Deserializer, MapAccess, SeqAccess, Visitor};

use super::tagged::{TagDispatch, deserialize_tagged};
use crate::model::*;

// ---- primitives ------------------------------------------------------------------------------

macro_rules! str_visitor {
    ($ty:ty, $expecting:literal, |$s:ident| $body:expr) => {
        impl<'de> Deserialize<'de> for $ty {
            fn deserialize<D: Deserializer<'de>>(d: D) -> Result<Self, D::Error> {
                struct V;
                impl Visitor<'_> for V {
                    type Value = $ty;
                    fn expecting(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
                        f.write_str($expecting)
                    }
                    fn visit_str<E: de::Error>(self, $s: &str) -> Result<Self::Value, E> {
                        $body.map_err(de::Error::custom)
                    }
                }
                d.deserialize_str(V)
            }
        }
    };
}

str_visitor!(
    SequenceString,
    "a sequenceString matching ^[A-Z*-]*$",
    |s| SequenceString::new(s)
);
str_visitor!(
    RefgetAccession,
    "a RefGet accession `SQ.<32 base64url chars>`",
    |s| { RefgetAccession::parse(s) }
);
str_visitor!(Digest, "a 32-character sha512t24u digest", |s| {
    Digest::parse(s)
});
str_visitor!(
    VrsIdentifier,
    "a GA4GH identifier `ga4gh:<prefix>.<digest>`",
    |s| { VrsIdentifier::parse(s) }
);
str_visitor!(
    Iri,
    "an IRI reference",
    |s| Ok::<_, std::convert::Infallible>(Iri::new(s))
);
str_visitor!(ResidueAlphabet, "one of \"aa\", \"na\"", |s| {
    ResidueAlphabet::parse(s)
});
str_visitor!(
    MoleculeType,
    "one of \"genomic\", \"RNA\", \"mRNA\", \"protein\"",
    |s| { MoleculeType::parse(s) }
);
str_visitor!(CopyChange, "a copyChange value", |s| CopyChange::parse(s));
str_visitor!(Syntax, "an Expression syntax", |s| Syntax::parse(s));
str_visitor!(Orientation, "\"forward\" or \"reverse_complement\"", |s| {
    Orientation::parse(s)
});
str_visitor!(AnchorOrientation, "\"left\" or \"right\"", |s| {
    AnchorOrientation::parse(s)
});

impl<'de> Deserialize<'de> for Range {
    fn deserialize<D: Deserializer<'de>>(d: D) -> Result<Self, D::Error> {
        struct V;
        impl<'de> Visitor<'de> for V {
            type Value = Range;
            fn expecting(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
                f.write_str("a Range: [min, max] with at least one integer")
            }
            fn visit_seq<A: SeqAccess<'de>>(self, mut seq: A) -> Result<Self::Value, A::Error> {
                let min: Option<i64> = seq
                    .next_element()?
                    .ok_or_else(|| de::Error::invalid_length(0, &self))?;
                let max: Option<i64> = seq
                    .next_element()?
                    .ok_or_else(|| de::Error::invalid_length(1, &self))?;
                if seq.next_element::<de::IgnoredAny>()?.is_some() {
                    return Err(de::Error::invalid_length(3, &self));
                }
                Range::new(min, max).map_err(de::Error::custom)
            }
        }
        d.deserialize_seq(V)
    }
}

impl<'de> Deserialize<'de> for IntOrRange {
    fn deserialize<D: Deserializer<'de>>(d: D) -> Result<Self, D::Error> {
        struct V;
        impl<'de> Visitor<'de> for V {
            type Value = IntOrRange;
            fn expecting(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
                f.write_str("an integer or a Range")
            }
            fn visit_i64<E: de::Error>(self, v: i64) -> Result<Self::Value, E> {
                if v == i64::MIN || v == i64::MAX {
                    return Err(de::Error::custom(
                        crate::error::CoordinateError::OutOfRange(i128::from(v)),
                    ));
                }
                Ok(IntOrRange::Int(v))
            }
            fn visit_u64<E: de::Error>(self, v: u64) -> Result<Self::Value, E> {
                IntOrRange::try_from(v).map_err(de::Error::custom)
            }
            fn visit_seq<A: SeqAccess<'de>>(self, seq: A) -> Result<Self::Value, A::Error> {
                Range::deserialize(de::value::SeqAccessDeserializer::new(seq))
                    .map(IntOrRange::Range)
            }
        }
        d.deserialize_any(V)
    }
}

impl<'de, T: Deserialize<'de>> Deserialize<'de> for IriOr<T> {
    fn deserialize<D: Deserializer<'de>>(d: D) -> Result<Self, D::Error> {
        struct V<T>(std::marker::PhantomData<T>);
        impl<'de, T: Deserialize<'de>> Visitor<'de> for V<T> {
            type Value = IriOr<T>;
            fn expecting(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
                f.write_str("an IRI reference string or an object")
            }
            fn visit_str<E: de::Error>(self, v: &str) -> Result<Self::Value, E> {
                Ok(IriOr::Iri(Iri::new(v)))
            }
            fn visit_map<A: MapAccess<'de>>(self, map: A) -> Result<Self::Value, A::Error> {
                T::deserialize(MapAccessDeserializer::new(map)).map(IriOr::Object)
            }
        }
        d.deserialize_any(V(std::marker::PhantomData))
    }
}

/// The `type` property of a concretely typed object, matched against the known class names
/// without allocating. Absent when the property was omitted (`#[serde(default)]`).
#[derive(Default)]
struct TypeTag(Option<&'static str>);

const TYPE_NAMES: &[&str] = &[
    "Allele",
    "RelativeAllele",
    "CisPhasedBlock",
    "Adjacency",
    "Terminus",
    "TraversalBlock",
    "DerivativeMolecule",
    "CopyNumberCount",
    "CopyNumberChange",
    "SequenceLocation",
    "RelativeSequenceLocation",
    "SequenceOffsetLocation",
    "SequenceReference",
    "LiteralSequenceExpression",
    "ReferenceLengthExpression",
    "LengthExpression",
];

impl<'de> Deserialize<'de> for TypeTag {
    fn deserialize<D: Deserializer<'de>>(d: D) -> Result<Self, D::Error> {
        struct V;
        impl Visitor<'_> for V {
            type Value = TypeTag;
            fn expecting(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
                f.write_str("a VRS class name")
            }
            fn visit_str<E: de::Error>(self, v: &str) -> Result<Self::Value, E> {
                TYPE_NAMES
                    .iter()
                    .find(|n| **n == v)
                    .map(|n| TypeTag(Some(n)))
                    .ok_or_else(|| de::Error::unknown_variant(v, TYPE_NAMES))
            }
        }
        d.deserialize_str(V)
    }
}

/// Check the `type` discriminator of a concretely typed object.
///
/// The JSON Schema lists `type` as required, but the official validation vectors (and the
/// reference implementation, which defaults it) omit it on nested objects whose class is
/// fixed by the containing property. Deserialization therefore accepts an absent `type` when
/// the class is already known and rejects a *wrong* one; serialization always emits it.
fn check_type<E: de::Error>(found: &TypeTag, expected: &'static str) -> Result<(), E> {
    match found.0 {
        None => Ok(()),
        Some(t) if t == expected => Ok(()),
        Some(t) => Err(de::Error::custom(crate::error::ModelError::TypeMismatch {
            expected,
            found: t.to_owned(),
        })),
    }
}

// ---- Extension / Expression -------------------------------------------------------------------

#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct ExtensionWire {
    id: Option<String>,
    #[serde(default)]
    extensions: Vec<Extension>,
    name: String,
    value: serde_json::Value,
    description: Option<String>,
}

impl<'de> Deserialize<'de> for Extension {
    fn deserialize<D: Deserializer<'de>>(d: D) -> Result<Self, D::Error> {
        let w = ExtensionWire::deserialize(d)?;
        let mut e = Extension::new(w.name, w.value);
        e.description = w.description;
        e.id = w.id;
        e.extensions = w.extensions;
        Ok(e)
    }
}

#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct ExpressionWire {
    id: Option<String>,
    #[serde(default)]
    extensions: Vec<Extension>,
    syntax: Syntax,
    value: String,
    syntax_version: Option<String>,
}

impl<'de> Deserialize<'de> for Expression {
    fn deserialize<D: Deserializer<'de>>(d: D) -> Result<Self, D::Error> {
        let w = ExpressionWire::deserialize(d)?;
        Ok(Expression::from_parts(
            w.syntax,
            w.value,
            w.syntax_version,
            w.id,
            w.extensions,
        ))
    }
}

// ---- wire structs -----------------------------------------------------------------------------

/// Declares a wire struct mirroring one JSON Schema class: the metadata properties allowed for
/// its `meta` kind, `type`, and its own properties.
macro_rules! wire {
    ($name:ident, meta = $kind:ident, { $( $(#[$fattr:meta])* $field:ident : $ty:ty ),* $(,)? }) => {
        #[derive(Deserialize)]
        #[serde(deny_unknown_fields, rename_all = "camelCase")]
        struct $name {
            #[serde(rename = "type", default)]
            type_name: TypeTag,
            id: Option<String>,
            name: Option<String>,
            description: Option<String>,
            #[serde(default)]
            aliases: Vec<String>,
            #[serde(default)]
            extensions: Vec<Extension>,
            #[serde(default, deserialize_with = "wire_digest_field", skip_serializing)]
            digest: DigestField<{ wire!(@allows_digest $kind) }>,
            #[serde(default, deserialize_with = "wire_expressions_field")]
            expressions: ExpressionsField<{ wire!(@allows_expressions $kind) }>,
            $( $(#[$fattr])* $field: $ty, )*
        }

        impl $name {
            fn take_meta(&mut self) -> Option<Box<Meta>> {
                let mut m = Meta::new();
                m.id = self.id.take();
                m.name = self.name.take();
                m.description = self.description.take();
                m.aliases = std::mem::take(&mut self.aliases);
                m.extensions = std::mem::take(&mut self.extensions);
                m.digest = self.digest.0.take();
                m.expressions = std::mem::take(&mut self.expressions.0);
                if m.is_empty() { None } else { Some(Box::new(m)) }
            }
        }
    };
    (@allows_digest entity) => { false };
    (@allows_digest identifiable) => { true };
    (@allows_digest variation) => { true };
    (@allows_expressions entity) => { false };
    (@allows_expressions identifiable) => { false };
    (@allows_expressions variation) => { true };
}

/// `digest` wire field; rejected (as an unknown property) for classes that do not have it.
#[derive(Default)]
struct DigestField<const ALLOWED: bool>(Option<Digest>);

fn wire_digest_field<'de, D: Deserializer<'de>, const ALLOWED: bool>(
    d: D,
) -> Result<DigestField<ALLOWED>, D::Error> {
    if !ALLOWED {
        return Err(de::Error::unknown_field("digest", &[]));
    }
    Option::<Digest>::deserialize(d).map(DigestField)
}

/// `expressions` wire field; rejected for non-variation classes.
#[derive(Default)]
struct ExpressionsField<const ALLOWED: bool>(Vec<Expression>);

fn wire_expressions_field<'de, D: Deserializer<'de>, const ALLOWED: bool>(
    d: D,
) -> Result<ExpressionsField<ALLOWED>, D::Error> {
    if !ALLOWED {
        return Err(de::Error::unknown_field("expressions", &[]));
    }
    Vec::<Expression>::deserialize(d).map(ExpressionsField)
}

wire!(SequenceReferenceWire, meta = entity, {
    refget_accession: RefgetAccession,
    residue_alphabet: Option<ResidueAlphabet>,
    sequence: Option<SequenceString>,
    molecule_type: Option<MoleculeType>,
    circular: Option<bool>,
});

impl<'de> Deserialize<'de> for SequenceReference {
    fn deserialize<D: Deserializer<'de>>(d: D) -> Result<Self, D::Error> {
        let mut w = SequenceReferenceWire::deserialize(d)?;
        check_type(&w.type_name, "SequenceReference")?;
        let meta = w.take_meta();
        Ok(SequenceReference::from_parts(
            w.refget_accession,
            w.residue_alphabet,
            w.molecule_type,
            w.circular,
            w.sequence,
            meta,
        ))
    }
}

wire!(LiteralSequenceExpressionWire, meta = entity, { sequence: SequenceString });

impl<'de> Deserialize<'de> for LiteralSequenceExpression {
    fn deserialize<D: Deserializer<'de>>(d: D) -> Result<Self, D::Error> {
        let mut w = LiteralSequenceExpressionWire::deserialize(d)?;
        check_type(&w.type_name, "LiteralSequenceExpression")?;
        let meta = w.take_meta();
        Ok(LiteralSequenceExpression::from_parts(w.sequence, meta))
    }
}

wire!(ReferenceLengthExpressionWire, meta = entity, {
    length: IntOrRange,
    sequence: Option<SequenceString>,
    repeat_subunit_length: i64,
});

impl<'de> Deserialize<'de> for ReferenceLengthExpression {
    fn deserialize<D: Deserializer<'de>>(d: D) -> Result<Self, D::Error> {
        let mut w = ReferenceLengthExpressionWire::deserialize(d)?;
        check_type(&w.type_name, "ReferenceLengthExpression")?;
        let meta = w.take_meta();
        ReferenceLengthExpression::from_parts(w.length, w.repeat_subunit_length, w.sequence, meta)
            .map_err(de::Error::custom)
    }
}

wire!(LengthExpressionWire, meta = entity, { length: Option<IntOrRange> });

impl<'de> Deserialize<'de> for LengthExpression {
    fn deserialize<D: Deserializer<'de>>(d: D) -> Result<Self, D::Error> {
        let mut w = LengthExpressionWire::deserialize(d)?;
        check_type(&w.type_name, "LengthExpression")?;
        let meta = w.take_meta();
        LengthExpression::from_parts(w.length, meta).map_err(de::Error::custom)
    }
}

wire!(SequenceLocationWire, meta = identifiable, {
    sequence_reference: Option<IriOr<SequenceReference>>,
    start: Option<IntOrRange>,
    end: Option<IntOrRange>,
    sequence: Option<SequenceString>,
});

impl<'de> Deserialize<'de> for SequenceLocation {
    fn deserialize<D: Deserializer<'de>>(d: D) -> Result<Self, D::Error> {
        let mut w = SequenceLocationWire::deserialize(d)?;
        check_type(&w.type_name, "SequenceLocation")?;
        let meta = w.take_meta();
        SequenceLocation::from_all_parts(w.sequence_reference, w.start, w.end, w.sequence, meta)
            .map_err(de::Error::custom)
    }
}

wire!(SequenceOffsetLocationWire, meta = entity, {
    sequence_reference: IriOr<SequenceReference>,
    anchor: i64,
    anchor_orientation: AnchorOrientation,
    offset_start: Option<IntOrRange>,
    offset_end: Option<IntOrRange>,
});

impl<'de> Deserialize<'de> for SequenceOffsetLocation {
    fn deserialize<D: Deserializer<'de>>(d: D) -> Result<Self, D::Error> {
        let mut w = SequenceOffsetLocationWire::deserialize(d)?;
        check_type(&w.type_name, "SequenceOffsetLocation")?;
        let meta = w.take_meta();
        let mut l = SequenceOffsetLocation::from_parts(
            w.sequence_reference,
            w.anchor,
            w.anchor_orientation,
            w.offset_start,
            w.offset_end,
        )
        .map_err(de::Error::custom)?;
        l.set_meta(meta);
        Ok(l)
    }
}

wire!(RelativeSequenceLocationWire, meta = identifiable, {
    base_sequence_location: IriOr<SequenceLocation>,
    mapped_sequence_location: IriOr<SequenceOffsetLocation>,
});

impl<'de> Deserialize<'de> for RelativeSequenceLocation {
    fn deserialize<D: Deserializer<'de>>(d: D) -> Result<Self, D::Error> {
        let mut w = RelativeSequenceLocationWire::deserialize(d)?;
        check_type(&w.type_name, "RelativeSequenceLocation")?;
        let meta = w.take_meta();
        let mut l =
            RelativeSequenceLocation::new(w.base_sequence_location, w.mapped_sequence_location);
        l.set_meta(meta);
        Ok(l)
    }
}

wire!(AlleleWire, meta = variation, {
    location: IriOr<SequenceLocation>,
    state: SequenceExpression,
});

impl<'de> Deserialize<'de> for Allele {
    fn deserialize<D: Deserializer<'de>>(d: D) -> Result<Self, D::Error> {
        let mut w = AlleleWire::deserialize(d)?;
        check_type(&w.type_name, "Allele")?;
        let meta = w.take_meta();
        let mut a = Allele::new(w.location, w.state);
        a.set_meta(meta);
        Ok(a)
    }
}

wire!(RelativeAlleleWire, meta = variation, {
    relative_location: IriOr<RelativeSequenceLocation>,
    base_state: SequenceExpression,
    mapped_state: SequenceExpression,
});

impl<'de> Deserialize<'de> for RelativeAllele {
    fn deserialize<D: Deserializer<'de>>(d: D) -> Result<Self, D::Error> {
        let mut w = RelativeAlleleWire::deserialize(d)?;
        check_type(&w.type_name, "RelativeAllele")?;
        let meta = w.take_meta();
        let mut a = RelativeAllele::new(w.relative_location, w.base_state, w.mapped_state);
        a.set_meta(meta);
        Ok(a)
    }
}

wire!(CisPhasedBlockWire, meta = variation, {
    members: Vec<IriOr<Allele>>,
    sequence_reference: Option<SequenceReference>,
});

impl<'de> Deserialize<'de> for CisPhasedBlock {
    fn deserialize<D: Deserializer<'de>>(d: D) -> Result<Self, D::Error> {
        let mut w = CisPhasedBlockWire::deserialize(d)?;
        check_type(&w.type_name, "CisPhasedBlock")?;
        let meta = w.take_meta();
        CisPhasedBlock::from_parts(w.members, w.sequence_reference, meta).map_err(de::Error::custom)
    }
}

wire!(AdjacencyWire, meta = variation, {
    adjoined_sequences: Vec<IriOr<Location>>,
    linker: Option<SequenceExpression>,
    homology: Option<bool>,
});

impl<'de> Deserialize<'de> for Adjacency {
    fn deserialize<D: Deserializer<'de>>(d: D) -> Result<Self, D::Error> {
        let mut w = AdjacencyWire::deserialize(d)?;
        check_type(&w.type_name, "Adjacency")?;
        let meta = w.take_meta();
        if w.adjoined_sequences.len() != 2 {
            return Err(de::Error::invalid_length(
                w.adjoined_sequences.len(),
                &"exactly 2 adjoined sequences",
            ));
        }
        Adjacency::from_parts(w.adjoined_sequences, w.linker, w.homology, meta)
            .map_err(de::Error::custom)
    }
}

wire!(TerminusWire, meta = variation, { location: IriOr<Location> });

impl<'de> Deserialize<'de> for Terminus {
    fn deserialize<D: Deserializer<'de>>(d: D) -> Result<Self, D::Error> {
        let mut w = TerminusWire::deserialize(d)?;
        check_type(&w.type_name, "Terminus")?;
        let meta = w.take_meta();
        let mut t = Terminus::new(w.location);
        t.set_meta(meta);
        Ok(t)
    }
}

wire!(TraversalBlockWire, meta = entity, {
    component: Option<Adjacency>,
    orientation: Option<Orientation>,
});

impl<'de> Deserialize<'de> for TraversalBlock {
    fn deserialize<D: Deserializer<'de>>(d: D) -> Result<Self, D::Error> {
        let mut w = TraversalBlockWire::deserialize(d)?;
        check_type(&w.type_name, "TraversalBlock")?;
        let meta = w.take_meta();
        let mut t = TraversalBlock::from_parts(w.component, w.orientation);
        t.set_meta(meta);
        Ok(t)
    }
}

wire!(DerivativeMoleculeWire, meta = variation, {
    components: Vec<DerivativeComponent>,
    circular: Option<bool>,
});

impl<'de> Deserialize<'de> for DerivativeMolecule {
    fn deserialize<D: Deserializer<'de>>(d: D) -> Result<Self, D::Error> {
        let mut w = DerivativeMoleculeWire::deserialize(d)?;
        check_type(&w.type_name, "DerivativeMolecule")?;
        let meta = w.take_meta();
        DerivativeMolecule::from_parts(w.components, w.circular, meta).map_err(de::Error::custom)
    }
}

wire!(CopyNumberCountWire, meta = variation, {
    location: IriOr<SequenceLocation>,
    copies: IntOrRange,
});

impl<'de> Deserialize<'de> for CopyNumberCount {
    fn deserialize<D: Deserializer<'de>>(d: D) -> Result<Self, D::Error> {
        let mut w = CopyNumberCountWire::deserialize(d)?;
        check_type(&w.type_name, "CopyNumberCount")?;
        let meta = w.take_meta();
        let mut c = CopyNumberCount::new(w.location, w.copies).map_err(de::Error::custom)?;
        c.set_meta(meta);
        Ok(c)
    }
}

wire!(CopyNumberChangeWire, meta = variation, {
    location: IriOr<SequenceLocation>,
    copy_change: CopyChange,
});

impl<'de> Deserialize<'de> for CopyNumberChange {
    fn deserialize<D: Deserializer<'de>>(d: D) -> Result<Self, D::Error> {
        let mut w = CopyNumberChangeWire::deserialize(d)?;
        check_type(&w.type_name, "CopyNumberChange")?;
        let meta = w.take_meta();
        let mut c = CopyNumberChange::new(w.location, w.copy_change);
        c.set_meta(meta);
        Ok(c)
    }
}

// ---- tagged unions ----------------------------------------------------------------------------

macro_rules! tagged_union {
    ($union:ident, $dispatch:ident, $expecting:literal: $( $variant:ident ),+ $(,)?) => {
        struct $dispatch;
        impl<'de> TagDispatch<'de> for $dispatch {
            type Output = $union;
            const VARIANTS: &'static [&'static str] = &[ $( stringify!($variant), )+ ];
            fn dispatch<D: Deserializer<'de>>(tag: &str, d: D) -> Result<$union, D::Error> {
                match tag {
                    $( stringify!($variant) => $variant::deserialize(d).map($union::from), )+
                    other => Err(de::Error::unknown_variant(other, Self::VARIANTS)),
                }
            }
        }
        impl<'de> Deserialize<'de> for $union {
            fn deserialize<D: Deserializer<'de>>(d: D) -> Result<Self, D::Error> {
                deserialize_tagged::<D, $dispatch>(d, $expecting)
            }
        }
    };
}

tagged_union!(Variation, VariationDispatch, "a Variation": Allele, RelativeAllele, CisPhasedBlock, Adjacency, Terminus, DerivativeMolecule, CopyNumberCount, CopyNumberChange);
tagged_union!(MolecularVariation, MolecularVariationDispatch, "a MolecularVariation": Allele, RelativeAllele, CisPhasedBlock, Adjacency, Terminus, DerivativeMolecule);
tagged_union!(SystemicVariation, SystemicVariationDispatch, "a SystemicVariation": CopyNumberCount, CopyNumberChange);

struct SequenceExpressionDispatch;
impl<'de> TagDispatch<'de> for SequenceExpressionDispatch {
    type Output = SequenceExpression;
    const VARIANTS: &'static [&'static str] = &[
        "LiteralSequenceExpression",
        "ReferenceLengthExpression",
        "LengthExpression",
    ];
    fn dispatch<D: Deserializer<'de>>(tag: &str, d: D) -> Result<SequenceExpression, D::Error> {
        match tag {
            "LiteralSequenceExpression" => {
                LiteralSequenceExpression::deserialize(d).map(SequenceExpression::Literal)
            }
            "ReferenceLengthExpression" => {
                ReferenceLengthExpression::deserialize(d).map(SequenceExpression::ReferenceLength)
            }
            "LengthExpression" => LengthExpression::deserialize(d).map(SequenceExpression::Length),
            other => Err(de::Error::unknown_variant(other, Self::VARIANTS)),
        }
    }
}
impl<'de> Deserialize<'de> for SequenceExpression {
    fn deserialize<D: Deserializer<'de>>(d: D) -> Result<Self, D::Error> {
        deserialize_tagged::<D, SequenceExpressionDispatch>(d, "a SequenceExpression")
    }
}

struct LocationDispatch;
impl<'de> TagDispatch<'de> for LocationDispatch {
    type Output = Location;
    const VARIANTS: &'static [&'static str] = &["SequenceLocation", "RelativeSequenceLocation"];
    fn dispatch<D: Deserializer<'de>>(tag: &str, d: D) -> Result<Location, D::Error> {
        match tag {
            "SequenceLocation" => SequenceLocation::deserialize(d).map(Location::Sequence),
            "RelativeSequenceLocation" => {
                RelativeSequenceLocation::deserialize(d).map(Location::from)
            }
            other => Err(de::Error::unknown_variant(other, Self::VARIANTS)),
        }
    }
}
impl<'de> Deserialize<'de> for Location {
    fn deserialize<D: Deserializer<'de>>(d: D) -> Result<Self, D::Error> {
        deserialize_tagged::<D, LocationDispatch>(d, "a Location")
    }
}

struct DerivativeComponentDispatch;
impl<'de> TagDispatch<'de> for DerivativeComponentDispatch {
    type Output = DerivativeComponent;
    const VARIANTS: &'static [&'static str] =
        &["Allele", "CisPhasedBlock", "Terminus", "TraversalBlock"];
    fn dispatch<D: Deserializer<'de>>(tag: &str, d: D) -> Result<DerivativeComponent, D::Error> {
        match tag {
            "Allele" => Allele::deserialize(d).map(DerivativeComponent::from),
            "CisPhasedBlock" => CisPhasedBlock::deserialize(d).map(DerivativeComponent::from),
            "Terminus" => Terminus::deserialize(d).map(DerivativeComponent::from),
            "TraversalBlock" => TraversalBlock::deserialize(d).map(DerivativeComponent::from),
            other => Err(de::Error::unknown_variant(other, Self::VARIANTS)),
        }
    }
}
impl<'de> Deserialize<'de> for DerivativeComponent {
    fn deserialize<D: Deserializer<'de>>(d: D) -> Result<Self, D::Error> {
        struct V;
        impl<'de> Visitor<'de> for V {
            type Value = DerivativeComponent;
            fn expecting(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
                f.write_str("an IRI reference or a DerivativeMolecule component object")
            }
            fn visit_str<E: de::Error>(self, v: &str) -> Result<Self::Value, E> {
                Ok(DerivativeComponent::Iri(Iri::new(v)))
            }
            fn visit_map<A: MapAccess<'de>>(self, map: A) -> Result<Self::Value, A::Error> {
                super::tagged::TaggedVisitor::<DerivativeComponentDispatch>::new(
                    "a DerivativeMolecule component",
                )
                .visit_map(map)
            }
        }
        d.deserialize_any(V)
    }
}
