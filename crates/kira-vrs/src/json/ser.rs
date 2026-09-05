//! `serde::Serialize` implementations for the model.

use serde::ser::{Serialize, SerializeSeq, SerializeStruct, Serializer};

use crate::model::*;

// ---- primitives ------------------------------------------------------------------------------

impl Serialize for SequenceString {
    fn serialize<S: Serializer>(&self, s: S) -> Result<S::Ok, S::Error> {
        s.serialize_str(self.as_str())
    }
}
impl Serialize for RefgetAccession {
    fn serialize<S: Serializer>(&self, s: S) -> Result<S::Ok, S::Error> {
        s.serialize_str(self.as_str())
    }
}
impl Serialize for Digest {
    fn serialize<S: Serializer>(&self, s: S) -> Result<S::Ok, S::Error> {
        s.serialize_str(self.as_str())
    }
}
impl Serialize for VrsIdentifier {
    fn serialize<S: Serializer>(&self, s: S) -> Result<S::Ok, S::Error> {
        s.collect_str(self)
    }
}
impl Serialize for Iri {
    fn serialize<S: Serializer>(&self, s: S) -> Result<S::Ok, S::Error> {
        s.serialize_str(self.as_str())
    }
}
impl Serialize for Range {
    fn serialize<S: Serializer>(&self, s: S) -> Result<S::Ok, S::Error> {
        let mut seq = s.serialize_seq(Some(2))?;
        seq.serialize_element(&self.min())?;
        seq.serialize_element(&self.max())?;
        seq.end()
    }
}
impl Serialize for IntOrRange {
    fn serialize<S: Serializer>(&self, s: S) -> Result<S::Ok, S::Error> {
        match self {
            Self::Int(v) => s.serialize_i64(*v),
            Self::Range(r) => r.serialize(s),
        }
    }
}
impl<T: Serialize> Serialize for IriOr<T> {
    fn serialize<S: Serializer>(&self, s: S) -> Result<S::Ok, S::Error> {
        match self {
            Self::Iri(i) => i.serialize(s),
            Self::Object(o) => o.serialize(s),
        }
    }
}

macro_rules! serialize_value_set {
    ($($ty:ty),+ $(,)?) => {
        $(
            impl Serialize for $ty {
                fn serialize<S: Serializer>(&self, s: S) -> Result<S::Ok, S::Error> {
                    s.serialize_str(self.as_str())
                }
            }
        )+
    };
}
serialize_value_set!(
    ResidueAlphabet,
    MoleculeType,
    CopyChange,
    Syntax,
    Orientation,
    AnchorOrientation
);

impl Serialize for Extension {
    fn serialize<S: Serializer>(&self, s: S) -> Result<S::Ok, S::Error> {
        let mut st = s.serialize_struct("Extension", 5)?;
        if let Some(id) = &self.id {
            st.serialize_field("id", id)?;
        }
        st.serialize_field("name", &self.name)?;
        st.serialize_field("value", &self.value)?;
        if let Some(d) = &self.description {
            st.serialize_field("description", d)?;
        }
        if !self.extensions.is_empty() {
            st.serialize_field("extensions", &self.extensions)?;
        }
        st.end()
    }
}

impl Serialize for Expression {
    fn serialize<S: Serializer>(&self, s: S) -> Result<S::Ok, S::Error> {
        let mut st = s.serialize_struct("Expression", 5)?;
        if let Some(id) = self.id() {
            st.serialize_field("id", id)?;
        }
        st.serialize_field("syntax", &self.syntax())?;
        st.serialize_field("value", self.value())?;
        if let Some(v) = self.syntax_version() {
            st.serialize_field("syntax_version", v)?;
        }
        if !self.extensions().is_empty() {
            st.serialize_field("extensions", self.extensions())?;
        }
        st.end()
    }
}

// ---- metadata helpers -------------------------------------------------------------------------

/// Which metadata properties a class may carry.
#[derive(Clone, Copy)]
enum MetaKind {
    /// gkm-core Entity: id, name, description, aliases, extensions.
    Entity,
    /// Entity + `digest`.
    Identifiable,
    /// Identifiable + `expressions`.
    Variation,
}

/// Emit `type`, `id` and `digest` (leading properties).
fn head<S: SerializeStruct>(
    st: &mut S,
    type_name: &'static str,
    meta: Option<&Meta>,
    kind: MetaKind,
) -> Result<(), S::Error> {
    st.serialize_field("type", type_name)?;
    if let Some(m) = meta {
        if let Some(id) = &m.id {
            st.serialize_field("id", id)?;
        }
        if !matches!(kind, MetaKind::Entity)
            && let Some(d) = &m.digest
        {
            st.serialize_field("digest", d)?;
        }
    }
    Ok(())
}

/// Emit the trailing decorative properties.
fn tail<S: SerializeStruct>(
    st: &mut S,
    meta: Option<&Meta>,
    kind: MetaKind,
) -> Result<(), S::Error> {
    let Some(m) = meta else { return Ok(()) };
    if let Some(v) = &m.name {
        st.serialize_field("name", v)?;
    }
    if let Some(v) = &m.description {
        st.serialize_field("description", v)?;
    }
    if !m.aliases.is_empty() {
        st.serialize_field("aliases", &m.aliases)?;
    }
    if !m.extensions.is_empty() {
        st.serialize_field("extensions", &m.extensions)?;
    }
    if matches!(kind, MetaKind::Variation) && !m.expressions.is_empty() {
        st.serialize_field("expressions", &m.expressions)?;
    }
    Ok(())
}

macro_rules! opt_field {
    ($st:expr, $name:literal, $value:expr) => {
        if let Some(v) = $value {
            $st.serialize_field($name, &v)?;
        }
    };
}

// ---- sequence classes ------------------------------------------------------------------------

impl Serialize for SequenceReference {
    fn serialize<S: Serializer>(&self, s: S) -> Result<S::Ok, S::Error> {
        let mut st = s.serialize_struct("SequenceReference", 10)?;
        head(&mut st, "SequenceReference", self.meta(), MetaKind::Entity)?;
        st.serialize_field("refgetAccession", self.refget_accession())?;
        opt_field!(st, "residueAlphabet", self.residue_alphabet());
        opt_field!(st, "sequence", self.sequence());
        opt_field!(st, "moleculeType", self.molecule_type());
        opt_field!(st, "circular", self.circular());
        tail(&mut st, self.meta(), MetaKind::Entity)?;
        st.end()
    }
}

impl Serialize for LiteralSequenceExpression {
    fn serialize<S: Serializer>(&self, s: S) -> Result<S::Ok, S::Error> {
        let mut st = s.serialize_struct("LiteralSequenceExpression", 7)?;
        head(
            &mut st,
            "LiteralSequenceExpression",
            self.meta(),
            MetaKind::Entity,
        )?;
        st.serialize_field("sequence", self.sequence())?;
        tail(&mut st, self.meta(), MetaKind::Entity)?;
        st.end()
    }
}

impl Serialize for ReferenceLengthExpression {
    fn serialize<S: Serializer>(&self, s: S) -> Result<S::Ok, S::Error> {
        let mut st = s.serialize_struct("ReferenceLengthExpression", 9)?;
        head(
            &mut st,
            "ReferenceLengthExpression",
            self.meta(),
            MetaKind::Entity,
        )?;
        st.serialize_field("length", &self.length())?;
        st.serialize_field("repeatSubunitLength", &self.repeat_subunit_length())?;
        opt_field!(st, "sequence", self.sequence());
        tail(&mut st, self.meta(), MetaKind::Entity)?;
        st.end()
    }
}

impl Serialize for LengthExpression {
    fn serialize<S: Serializer>(&self, s: S) -> Result<S::Ok, S::Error> {
        let mut st = s.serialize_struct("LengthExpression", 7)?;
        head(&mut st, "LengthExpression", self.meta(), MetaKind::Entity)?;
        opt_field!(st, "length", self.length());
        tail(&mut st, self.meta(), MetaKind::Entity)?;
        st.end()
    }
}

impl Serialize for SequenceExpression {
    fn serialize<S: Serializer>(&self, s: S) -> Result<S::Ok, S::Error> {
        match self {
            Self::Literal(e) => e.serialize(s),
            Self::ReferenceLength(e) => e.serialize(s),
            Self::Length(e) => e.serialize(s),
        }
    }
}

// ---- locations -------------------------------------------------------------------------------

impl Serialize for SequenceLocation {
    fn serialize<S: Serializer>(&self, s: S) -> Result<S::Ok, S::Error> {
        let mut st = s.serialize_struct("SequenceLocation", 11)?;
        head(
            &mut st,
            "SequenceLocation",
            self.meta(),
            MetaKind::Identifiable,
        )?;
        opt_field!(st, "sequenceReference", self.sequence_reference());
        opt_field!(st, "start", self.start());
        opt_field!(st, "end", self.end());
        opt_field!(st, "sequence", self.sequence());
        tail(&mut st, self.meta(), MetaKind::Identifiable)?;
        st.end()
    }
}

impl Serialize for SequenceOffsetLocation {
    fn serialize<S: Serializer>(&self, s: S) -> Result<S::Ok, S::Error> {
        let mut st = s.serialize_struct("SequenceOffsetLocation", 11)?;
        head(
            &mut st,
            "SequenceOffsetLocation",
            self.meta(),
            MetaKind::Entity,
        )?;
        st.serialize_field("sequenceReference", self.sequence_reference())?;
        st.serialize_field("anchor", &self.anchor())?;
        st.serialize_field("anchorOrientation", &self.anchor_orientation())?;
        opt_field!(st, "offsetStart", self.offset_start());
        opt_field!(st, "offsetEnd", self.offset_end());
        tail(&mut st, self.meta(), MetaKind::Entity)?;
        st.end()
    }
}

impl Serialize for RelativeSequenceLocation {
    fn serialize<S: Serializer>(&self, s: S) -> Result<S::Ok, S::Error> {
        let mut st = s.serialize_struct("RelativeSequenceLocation", 9)?;
        head(
            &mut st,
            "RelativeSequenceLocation",
            self.meta(),
            MetaKind::Identifiable,
        )?;
        st.serialize_field("baseSequenceLocation", self.base_sequence_location())?;
        st.serialize_field("mappedSequenceLocation", self.mapped_sequence_location())?;
        tail(&mut st, self.meta(), MetaKind::Identifiable)?;
        st.end()
    }
}

impl Serialize for Location {
    fn serialize<S: Serializer>(&self, s: S) -> Result<S::Ok, S::Error> {
        match self {
            Self::Sequence(l) => l.serialize(s),
            Self::RelativeSequence(l) => l.serialize(s),
        }
    }
}

// ---- variation -------------------------------------------------------------------------------

impl Serialize for Allele {
    fn serialize<S: Serializer>(&self, s: S) -> Result<S::Ok, S::Error> {
        let mut st = s.serialize_struct("Allele", 10)?;
        head(&mut st, "Allele", self.meta(), MetaKind::Variation)?;
        st.serialize_field("location", self.location())?;
        st.serialize_field("state", self.state())?;
        tail(&mut st, self.meta(), MetaKind::Variation)?;
        st.end()
    }
}

impl Serialize for RelativeAllele {
    fn serialize<S: Serializer>(&self, s: S) -> Result<S::Ok, S::Error> {
        let mut st = s.serialize_struct("RelativeAllele", 11)?;
        head(&mut st, "RelativeAllele", self.meta(), MetaKind::Variation)?;
        st.serialize_field("relativeLocation", self.relative_location())?;
        st.serialize_field("baseState", self.base_state())?;
        st.serialize_field("mappedState", self.mapped_state())?;
        tail(&mut st, self.meta(), MetaKind::Variation)?;
        st.end()
    }
}

impl Serialize for CisPhasedBlock {
    fn serialize<S: Serializer>(&self, s: S) -> Result<S::Ok, S::Error> {
        let mut st = s.serialize_struct("CisPhasedBlock", 10)?;
        head(&mut st, "CisPhasedBlock", self.meta(), MetaKind::Variation)?;
        st.serialize_field("members", self.members())?;
        opt_field!(st, "sequenceReference", self.sequence_reference());
        tail(&mut st, self.meta(), MetaKind::Variation)?;
        st.end()
    }
}

impl Serialize for Adjacency {
    fn serialize<S: Serializer>(&self, s: S) -> Result<S::Ok, S::Error> {
        let mut st = s.serialize_struct("Adjacency", 11)?;
        head(&mut st, "Adjacency", self.meta(), MetaKind::Variation)?;
        st.serialize_field("adjoinedSequences", &self.adjoined_sequences()[..])?;
        opt_field!(st, "linker", self.linker());
        opt_field!(st, "homology", self.homology());
        tail(&mut st, self.meta(), MetaKind::Variation)?;
        st.end()
    }
}

impl Serialize for Terminus {
    fn serialize<S: Serializer>(&self, s: S) -> Result<S::Ok, S::Error> {
        let mut st = s.serialize_struct("Terminus", 9)?;
        head(&mut st, "Terminus", self.meta(), MetaKind::Variation)?;
        st.serialize_field("location", self.location())?;
        tail(&mut st, self.meta(), MetaKind::Variation)?;
        st.end()
    }
}

impl Serialize for TraversalBlock {
    fn serialize<S: Serializer>(&self, s: S) -> Result<S::Ok, S::Error> {
        let mut st = s.serialize_struct("TraversalBlock", 8)?;
        head(&mut st, "TraversalBlock", self.meta(), MetaKind::Entity)?;
        opt_field!(st, "component", self.component());
        opt_field!(st, "orientation", self.orientation());
        tail(&mut st, self.meta(), MetaKind::Entity)?;
        st.end()
    }
}

impl Serialize for DerivativeComponent {
    fn serialize<S: Serializer>(&self, s: S) -> Result<S::Ok, S::Error> {
        match self {
            Self::Iri(i) => i.serialize(s),
            Self::Allele(a) => a.serialize(s),
            Self::CisPhasedBlock(c) => c.serialize(s),
            Self::Terminus(t) => t.serialize(s),
            Self::TraversalBlock(t) => t.serialize(s),
        }
    }
}

impl Serialize for DerivativeMolecule {
    fn serialize<S: Serializer>(&self, s: S) -> Result<S::Ok, S::Error> {
        let mut st = s.serialize_struct("DerivativeMolecule", 10)?;
        head(
            &mut st,
            "DerivativeMolecule",
            self.meta(),
            MetaKind::Variation,
        )?;
        st.serialize_field("components", self.components())?;
        opt_field!(st, "circular", self.circular());
        tail(&mut st, self.meta(), MetaKind::Variation)?;
        st.end()
    }
}

impl Serialize for CopyNumberCount {
    fn serialize<S: Serializer>(&self, s: S) -> Result<S::Ok, S::Error> {
        let mut st = s.serialize_struct("CopyNumberCount", 10)?;
        head(&mut st, "CopyNumberCount", self.meta(), MetaKind::Variation)?;
        st.serialize_field("location", self.location())?;
        st.serialize_field("copies", &self.copies())?;
        tail(&mut st, self.meta(), MetaKind::Variation)?;
        st.end()
    }
}

impl Serialize for CopyNumberChange {
    fn serialize<S: Serializer>(&self, s: S) -> Result<S::Ok, S::Error> {
        let mut st = s.serialize_struct("CopyNumberChange", 10)?;
        head(
            &mut st,
            "CopyNumberChange",
            self.meta(),
            MetaKind::Variation,
        )?;
        st.serialize_field("location", self.location())?;
        st.serialize_field("copyChange", &self.copy_change())?;
        tail(&mut st, self.meta(), MetaKind::Variation)?;
        st.end()
    }
}

macro_rules! serialize_union {
    ($union:ident: $( $variant:ident ),+ $(,)?) => {
        impl Serialize for $union {
            fn serialize<S: Serializer>(&self, s: S) -> Result<S::Ok, S::Error> {
                match self {
                    $( Self::$variant(v) => v.serialize(s), )+
                }
            }
        }
    };
}
serialize_union!(MolecularVariation: Allele, RelativeAllele, CisPhasedBlock, Adjacency, Terminus, DerivativeMolecule);
serialize_union!(SystemicVariation: CopyNumberCount, CopyNumberChange);
serialize_union!(Variation: Allele, RelativeAllele, CisPhasedBlock, Adjacency, Terminus, DerivativeMolecule, CopyNumberCount, CopyNumberChange);
