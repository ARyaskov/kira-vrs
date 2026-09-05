//! Variation classes: `Allele`, `RelativeAllele`, `CisPhasedBlock`, `Adjacency`, `Terminus`,
//! `DerivativeMolecule`, `CopyNumberCount`, `CopyNumberChange` and the polymorphic unions.
//!
//! In the polymorphic unions ([`Variation`], [`MolecularVariation`], [`DerivativeComponent`])
//! the large structural classes are boxed so that a `Vec<Variation>` of small variants (the
//! overwhelmingly common case) does not pay for the largest class; `From` conversions hide
//! the boxing.

use crate::error::ModelError;
use crate::model::location::{Location, RelativeSequenceLocation, SequenceLocation};
use crate::model::meta::{Meta, impl_entity, impl_variation_expressions};
use crate::model::primitives::{CopyChange, IntOrRange, Iri, IriOr, Orientation};
use crate::model::sequence::{SequenceExpression, SequenceReference};

// ---------------------------------------------------------------------------------------------
// Allele
// ---------------------------------------------------------------------------------------------

/// A VRS `Allele`: the state of a molecule at a location — the workhorse class for SNVs,
/// indels and other contiguous variation.
#[derive(Clone, Debug, PartialEq)]
pub struct Allele {
    location: IriOr<SequenceLocation>,
    state: SequenceExpression,
    meta: Option<Box<Meta>>,
}

impl Allele {
    /// An allele with the given location and state.
    pub fn new(
        location: impl Into<IriOr<SequenceLocation>>,
        state: impl Into<SequenceExpression>,
    ) -> Self {
        Self {
            location: location.into(),
            state: state.into(),
            meta: None,
        }
    }

    /// `location`.
    #[inline]
    pub fn location(&self) -> &IriOr<SequenceLocation> {
        &self.location
    }

    /// The inline sequence location, if the location is not an IRI.
    #[inline]
    pub fn sequence_location(&self) -> Option<&SequenceLocation> {
        self.location.as_object()
    }

    /// `state`.
    #[inline]
    pub fn state(&self) -> &SequenceExpression {
        &self.state
    }

    /// Decompose into location and state.
    pub fn into_parts(self) -> (IriOr<SequenceLocation>, SequenceExpression) {
        (self.location, self.state)
    }

    /// Rebuild with a different location and state, keeping metadata (used by normalization).
    pub(crate) fn rebuilt(&self, location: SequenceLocation, state: SequenceExpression) -> Self {
        Self {
            location: IriOr::Object(location),
            state,
            meta: self.meta.clone(),
        }
    }

    pub(crate) fn set_meta(&mut self, meta: Option<Box<Meta>>) {
        self.meta = meta;
    }
}
impl_entity!(Allele);
impl_variation_expressions!(Allele);

// ---------------------------------------------------------------------------------------------
// RelativeAllele
// ---------------------------------------------------------------------------------------------

/// A VRS `RelativeAllele` (draft): an allele on a mapped location relative to a base
/// location, used for intronic variants expressed against a transcript.
#[derive(Clone, Debug, PartialEq)]
pub struct RelativeAllele {
    relative_location: IriOr<RelativeSequenceLocation>,
    base_state: SequenceExpression,
    mapped_state: SequenceExpression,
    meta: Option<Box<Meta>>,
}

impl RelativeAllele {
    /// A relative allele with the given location and base / mapped states.
    pub fn new(
        relative_location: impl Into<IriOr<RelativeSequenceLocation>>,
        base_state: impl Into<SequenceExpression>,
        mapped_state: impl Into<SequenceExpression>,
    ) -> Self {
        Self {
            relative_location: relative_location.into(),
            base_state: base_state.into(),
            mapped_state: mapped_state.into(),
            meta: None,
        }
    }

    /// `relativeLocation`.
    #[inline]
    pub fn relative_location(&self) -> &IriOr<RelativeSequenceLocation> {
        &self.relative_location
    }

    /// `baseState`: the state on the base (e.g. genomic) sequence.
    #[inline]
    pub fn base_state(&self) -> &SequenceExpression {
        &self.base_state
    }

    /// `mappedState`: the state on the mapped (e.g. transcript) sequence.
    #[inline]
    pub fn mapped_state(&self) -> &SequenceExpression {
        &self.mapped_state
    }

    pub(crate) fn set_meta(&mut self, meta: Option<Box<Meta>>) {
        self.meta = meta;
    }
}
impl_entity!(RelativeAllele);
impl_variation_expressions!(RelativeAllele);

// ---------------------------------------------------------------------------------------------
// CisPhasedBlock
// ---------------------------------------------------------------------------------------------

/// A VRS `CisPhasedBlock`: alleles known to occur on the same molecule (a haplotype).
#[derive(Clone, Debug, PartialEq)]
pub struct CisPhasedBlock {
    members: Vec<IriOr<Allele>>,
    sequence_reference: Option<SequenceReference>,
    meta: Option<Box<Meta>>,
}

impl CisPhasedBlock {
    /// A block of at least two in-cis alleles.
    ///
    /// # Errors
    /// [`ModelError::TooFewItems`] with fewer than two members.
    pub fn new(members: Vec<IriOr<Allele>>) -> Result<Self, ModelError> {
        if members.len() < 2 {
            return Err(ModelError::TooFewItems {
                class: "CisPhasedBlock",
                property: "members",
                min: 2,
                actual: members.len(),
            });
        }
        Ok(Self {
            members,
            sequence_reference: None,
            meta: None,
        })
    }

    /// Set the shared `sequenceReference` on which all members lie (member locations may then
    /// omit their own reference).
    #[must_use]
    pub fn with_sequence_reference(mut self, reference: SequenceReference) -> Self {
        self.sequence_reference = Some(reference);
        self
    }

    /// `members`.
    #[inline]
    pub fn members(&self) -> &[IriOr<Allele>] {
        &self.members
    }

    /// `sequenceReference`.
    #[inline]
    pub fn sequence_reference(&self) -> Option<&SequenceReference> {
        self.sequence_reference.as_ref()
    }

    pub(crate) fn from_parts(
        members: Vec<IriOr<Allele>>,
        sequence_reference: Option<SequenceReference>,
        meta: Option<Box<Meta>>,
    ) -> Result<Self, ModelError> {
        let mut b = Self::new(members)?;
        b.sequence_reference = sequence_reference;
        b.meta = meta;
        Ok(b)
    }

    pub(crate) fn members_mut(&mut self) -> &mut Vec<IriOr<Allele>> {
        &mut self.members
    }
}
impl_entity!(CisPhasedBlock);
impl_variation_expressions!(CisPhasedBlock);

// ---------------------------------------------------------------------------------------------
// Adjacency / Terminus / DerivativeMolecule
// ---------------------------------------------------------------------------------------------

/// A VRS `Adjacency`: the junction of the end of one sequence with the start of another,
/// optionally through a linker sequence — the core structural-variation concept.
///
/// Each adjoined sequence is a half-open location defined by *either* `start` (extending
/// rightwards) *or* `end` (extending leftwards), never both.
#[derive(Clone, Debug, PartialEq)]
pub struct Adjacency {
    adjoined_sequences: [IriOr<Location>; 2],
    linker: Option<SequenceExpression>,
    homology: Option<bool>,
    meta: Option<Box<Meta>>,
}

impl Adjacency {
    /// An adjacency between two adjoined sequences.
    ///
    /// # Errors
    /// [`ModelError::AdjoinedSequenceHasStartAndEnd`] if an inline sequence location defines
    /// both `start` and `end`.
    pub fn new(
        first: impl Into<IriOr<Location>>,
        second: impl Into<IriOr<Location>>,
    ) -> Result<Self, ModelError> {
        let adjoined_sequences = [first.into(), second.into()];
        for a in &adjoined_sequences {
            if let IriOr::Object(Location::Sequence(l)) = a
                && l.start().is_some()
                && l.end().is_some()
            {
                return Err(ModelError::AdjoinedSequenceHasStartAndEnd);
            }
        }
        Ok(Self {
            adjoined_sequences,
            linker: None,
            homology: None,
            meta: None,
        })
    }

    /// Set the `linker` sequence found between the adjoined sequences.
    #[must_use]
    pub fn with_linker(mut self, linker: impl Into<SequenceExpression>) -> Self {
        self.linker = Some(linker.into());
        self
    }

    /// Set the draft `homology` flag.
    #[must_use]
    pub fn with_homology(mut self, homology: bool) -> Self {
        self.homology = Some(homology);
        self
    }

    /// `adjoinedSequences` (ordered).
    #[inline]
    pub fn adjoined_sequences(&self) -> &[IriOr<Location>; 2] {
        &self.adjoined_sequences
    }

    /// `linker`.
    #[inline]
    pub fn linker(&self) -> Option<&SequenceExpression> {
        self.linker.as_ref()
    }

    /// `homology`.
    #[inline]
    pub fn homology(&self) -> Option<bool> {
        self.homology
    }

    pub(crate) fn from_parts(
        adjoined: Vec<IriOr<Location>>,
        linker: Option<SequenceExpression>,
        homology: Option<bool>,
        meta: Option<Box<Meta>>,
    ) -> Result<Self, ModelError> {
        let [first, second]: [IriOr<Location>; 2] =
            adjoined
                .try_into()
                .map_err(|v: Vec<_>| ModelError::TooFewItems {
                    class: "Adjacency",
                    property: "adjoinedSequences",
                    min: 2,
                    actual: v.len(),
                })?;
        let mut a = Self::new(first, second)?;
        a.linker = linker;
        a.homology = homology;
        a.meta = meta;
        Ok(a)
    }

    /// Rebuild with swapped adjoined sequences, keeping other properties (normalization).
    pub(crate) fn with_adjoined(&self, adjoined: [IriOr<Location>; 2]) -> Self {
        Self {
            adjoined_sequences: adjoined,
            linker: self.linker.clone(),
            homology: self.homology,
            meta: self.meta.clone(),
        }
    }
}
impl_entity!(Adjacency);
impl_variation_expressions!(Adjacency);

/// A VRS `Terminus` (draft): the end of a molecule, described by a half-open location.
#[derive(Clone, Debug, PartialEq)]
pub struct Terminus {
    location: IriOr<Location>,
    meta: Option<Box<Meta>>,
}

impl Terminus {
    /// A terminus at the given location.
    pub fn new(location: impl Into<IriOr<Location>>) -> Self {
        Self {
            location: location.into(),
            meta: None,
        }
    }

    /// `location`.
    #[inline]
    pub fn location(&self) -> &IriOr<Location> {
        &self.location
    }

    pub(crate) fn set_meta(&mut self, meta: Option<Box<Meta>>) {
        self.meta = meta;
    }
}
impl_entity!(Terminus);
impl_variation_expressions!(Terminus);

/// A VRS `TraversalBlock` (draft): an adjacency with an orientation, used inside a
/// [`DerivativeMolecule`] to resolve strand traversal.
#[derive(Clone, Debug, PartialEq)]
pub struct TraversalBlock {
    component: Option<Adjacency>,
    orientation: Option<Orientation>,
    meta: Option<Box<Meta>>,
}

impl TraversalBlock {
    /// An oriented adjacency.
    pub fn new(component: Adjacency, orientation: Orientation) -> Self {
        Self {
            component: Some(component),
            orientation: Some(orientation),
            meta: None,
        }
    }

    /// General constructor (both properties optional, as in the schema).
    pub fn from_parts(component: Option<Adjacency>, orientation: Option<Orientation>) -> Self {
        Self {
            component,
            orientation,
            meta: None,
        }
    }

    /// `component`.
    #[inline]
    pub fn component(&self) -> Option<&Adjacency> {
        self.component.as_ref()
    }

    /// `orientation`.
    #[inline]
    pub fn orientation(&self) -> Option<Orientation> {
        self.orientation
    }

    pub(crate) fn set_meta(&mut self, meta: Option<Box<Meta>>) {
        self.meta = meta;
    }
}
impl_entity!(TraversalBlock);

/// A component of a [`DerivativeMolecule`].
#[derive(Clone, Debug, PartialEq)]
pub enum DerivativeComponent {
    /// A reference to a component defined elsewhere.
    Iri(Iri),
    /// An allele on the derived molecule.
    Allele(Allele),
    /// A haplotype on the derived molecule.
    CisPhasedBlock(CisPhasedBlock),
    /// A molecule end.
    Terminus(Box<Terminus>),
    /// An oriented adjacency.
    TraversalBlock(Box<TraversalBlock>),
}

impl DerivativeComponent {
    /// The VRS class name of the component (`"iriReference"` for IRIs).
    pub const fn type_name(&self) -> &'static str {
        match self {
            Self::Iri(_) => "iriReference",
            Self::Allele(_) => "Allele",
            Self::CisPhasedBlock(_) => "CisPhasedBlock",
            Self::Terminus(_) => "Terminus",
            Self::TraversalBlock(_) => "TraversalBlock",
        }
    }
}

impl From<Iri> for DerivativeComponent {
    fn from(v: Iri) -> Self {
        Self::Iri(v)
    }
}
impl From<Allele> for DerivativeComponent {
    fn from(v: Allele) -> Self {
        Self::Allele(v)
    }
}
impl From<CisPhasedBlock> for DerivativeComponent {
    fn from(v: CisPhasedBlock) -> Self {
        Self::CisPhasedBlock(v)
    }
}
impl From<Terminus> for DerivativeComponent {
    fn from(v: Terminus) -> Self {
        Self::Terminus(Box::new(v))
    }
}
impl From<TraversalBlock> for DerivativeComponent {
    fn from(v: TraversalBlock) -> Self {
        Self::TraversalBlock(Box::new(v))
    }
}

/// A VRS `DerivativeMolecule` (draft): a molecule assembled from ordered components of other
/// molecules, typically the product of structural variation.
#[derive(Clone, Debug, PartialEq)]
pub struct DerivativeMolecule {
    components: Vec<DerivativeComponent>,
    circular: Option<bool>,
    meta: Option<Box<Meta>>,
}

impl DerivativeMolecule {
    /// A molecule of at least two ordered components.
    ///
    /// # Errors
    /// [`ModelError::TooFewItems`] with fewer than two components.
    pub fn new(components: Vec<DerivativeComponent>) -> Result<Self, ModelError> {
        if components.len() < 2 {
            return Err(ModelError::TooFewItems {
                class: "DerivativeMolecule",
                property: "components",
                min: 2,
                actual: components.len(),
            });
        }
        Ok(Self {
            components,
            circular: None,
            meta: None,
        })
    }

    /// Set `circular`.
    #[must_use]
    pub fn with_circular(mut self, circular: bool) -> Self {
        self.circular = Some(circular);
        self
    }

    /// `components` (ordered).
    #[inline]
    pub fn components(&self) -> &[DerivativeComponent] {
        &self.components
    }

    /// `circular`.
    #[inline]
    pub fn circular(&self) -> Option<bool> {
        self.circular
    }

    pub(crate) fn from_parts(
        components: Vec<DerivativeComponent>,
        circular: Option<bool>,
        meta: Option<Box<Meta>>,
    ) -> Result<Self, ModelError> {
        let mut m = Self::new(components)?;
        m.circular = circular;
        m.meta = meta;
        Ok(m)
    }
}
impl_entity!(DerivativeMolecule);
impl_variation_expressions!(DerivativeMolecule);

// ---------------------------------------------------------------------------------------------
// Systemic variation
// ---------------------------------------------------------------------------------------------

/// A VRS `CopyNumberCount`: the absolute number of copies of a location in a system.
#[derive(Clone, Debug, PartialEq)]
pub struct CopyNumberCount {
    location: IriOr<SequenceLocation>,
    copies: IntOrRange,
    meta: Option<Box<Meta>>,
}

impl CopyNumberCount {
    /// `copies` copies of `location`.
    ///
    /// # Errors
    /// [`ModelError::NegativeCount`] for negative copy numbers.
    pub fn new(
        location: impl Into<IriOr<SequenceLocation>>,
        copies: impl Into<IntOrRange>,
    ) -> Result<Self, ModelError> {
        let copies = copies.into();
        if copies.has_negative() {
            return Err(ModelError::NegativeCount {
                property: "copies",
                value: copies.lower_bound().or(copies.upper_bound()).unwrap_or(-1),
            });
        }
        Ok(Self {
            location: location.into(),
            copies,
            meta: None,
        })
    }

    /// `location`.
    #[inline]
    pub fn location(&self) -> &IriOr<SequenceLocation> {
        &self.location
    }

    /// `copies`.
    #[inline]
    pub fn copies(&self) -> IntOrRange {
        self.copies
    }

    pub(crate) fn set_meta(&mut self, meta: Option<Box<Meta>>) {
        self.meta = meta;
    }
}
impl_entity!(CopyNumberCount);
impl_variation_expressions!(CopyNumberCount);

/// A VRS `CopyNumberChange` (draft): a categorical copy-number assessment of a location
/// relative to baseline ploidy.
#[derive(Clone, Debug, PartialEq)]
pub struct CopyNumberChange {
    location: IriOr<SequenceLocation>,
    copy_change: CopyChange,
    meta: Option<Box<Meta>>,
}

impl CopyNumberChange {
    /// A copy-number change at `location`.
    pub fn new(location: impl Into<IriOr<SequenceLocation>>, copy_change: CopyChange) -> Self {
        Self {
            location: location.into(),
            copy_change,
            meta: None,
        }
    }

    /// `location`.
    #[inline]
    pub fn location(&self) -> &IriOr<SequenceLocation> {
        &self.location
    }

    /// `copyChange`.
    #[inline]
    pub fn copy_change(&self) -> CopyChange {
        self.copy_change
    }

    pub(crate) fn set_meta(&mut self, meta: Option<Box<Meta>>) {
        self.meta = meta;
    }
}
impl_entity!(CopyNumberChange);
impl_variation_expressions!(CopyNumberChange);

// ---------------------------------------------------------------------------------------------
// Polymorphic unions
// ---------------------------------------------------------------------------------------------

/// A VRS `MolecularVariation`: variation on a contiguous molecule.
#[derive(Clone, Debug, PartialEq)]
pub enum MolecularVariation {
    /// `Allele`.
    Allele(Allele),
    /// `RelativeAllele` (draft).
    RelativeAllele(Box<RelativeAllele>),
    /// `CisPhasedBlock`.
    CisPhasedBlock(CisPhasedBlock),
    /// `Adjacency`.
    Adjacency(Box<Adjacency>),
    /// `Terminus` (draft).
    Terminus(Box<Terminus>),
    /// `DerivativeMolecule` (draft).
    DerivativeMolecule(Box<DerivativeMolecule>),
}

/// A VRS `SystemicVariation`: variation of a location across a system (genome, cell).
#[derive(Clone, Debug, PartialEq)]
pub enum SystemicVariation {
    /// `CopyNumberCount`.
    CopyNumberCount(CopyNumberCount),
    /// `CopyNumberChange` (draft).
    CopyNumberChange(CopyNumberChange),
}

/// A VRS `Variation`: the root of every variation class.
#[derive(Clone, Debug, PartialEq)]
pub enum Variation {
    /// `Allele`.
    Allele(Allele),
    /// `RelativeAllele` (draft).
    RelativeAllele(Box<RelativeAllele>),
    /// `CisPhasedBlock`.
    CisPhasedBlock(CisPhasedBlock),
    /// `Adjacency`.
    Adjacency(Box<Adjacency>),
    /// `Terminus` (draft).
    Terminus(Box<Terminus>),
    /// `DerivativeMolecule` (draft).
    DerivativeMolecule(Box<DerivativeMolecule>),
    /// `CopyNumberCount`.
    CopyNumberCount(CopyNumberCount),
    /// `CopyNumberChange` (draft).
    CopyNumberChange(CopyNumberChange),
}

macro_rules! union_impls {
    ($union:ident; inline: [$( $iv:ident ),* $(,)?]; boxed: [$( $bv:ident ),* $(,)?]) => {
        $(
            impl From<$iv> for $union {
                fn from(v: $iv) -> Self {
                    Self::$iv(v)
                }
            }
        )*
        $(
            impl From<$bv> for $union {
                fn from(v: $bv) -> Self {
                    Self::$bv(Box::new(v))
                }
            }
            impl From<Box<$bv>> for $union {
                fn from(v: Box<$bv>) -> Self {
                    Self::$bv(v)
                }
            }
        )*
        impl $union {
            /// The VRS class name of the concrete object.
            pub const fn type_name(&self) -> &'static str {
                match self {
                    $( Self::$iv(_) => stringify!($iv), )*
                    $( Self::$bv(_) => stringify!($bv), )*
                }
            }
        }
    };
}

union_impls!(MolecularVariation;
    inline: [Allele, CisPhasedBlock];
    boxed: [RelativeAllele, Adjacency, Terminus, DerivativeMolecule]);
union_impls!(SystemicVariation; inline: [CopyNumberCount, CopyNumberChange]; boxed: []);
union_impls!(Variation;
    inline: [Allele, CisPhasedBlock, CopyNumberCount, CopyNumberChange];
    boxed: [RelativeAllele, Adjacency, Terminus, DerivativeMolecule]);

impl From<MolecularVariation> for Variation {
    fn from(v: MolecularVariation) -> Self {
        match v {
            MolecularVariation::Allele(x) => Self::Allele(x),
            MolecularVariation::RelativeAllele(x) => Self::RelativeAllele(x),
            MolecularVariation::CisPhasedBlock(x) => Self::CisPhasedBlock(x),
            MolecularVariation::Adjacency(x) => Self::Adjacency(x),
            MolecularVariation::Terminus(x) => Self::Terminus(x),
            MolecularVariation::DerivativeMolecule(x) => Self::DerivativeMolecule(x),
        }
    }
}

impl From<SystemicVariation> for Variation {
    fn from(v: SystemicVariation) -> Self {
        match v {
            SystemicVariation::CopyNumberCount(x) => Self::CopyNumberCount(x),
            SystemicVariation::CopyNumberChange(x) => Self::CopyNumberChange(x),
        }
    }
}

impl Variation {
    /// The allele, if this is one.
    #[inline]
    pub fn as_allele(&self) -> Option<&Allele> {
        match self {
            Self::Allele(a) => Some(a),
            _ => None,
        }
    }

    /// `true` for the molecular variation classes.
    pub const fn is_molecular(&self) -> bool {
        !self.is_systemic()
    }

    /// `true` for the systemic variation classes.
    pub const fn is_systemic(&self) -> bool {
        matches!(self, Self::CopyNumberCount(_) | Self::CopyNumberChange(_))
    }
}
