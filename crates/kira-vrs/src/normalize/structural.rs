//! Normalization conventions for `CisPhasedBlock`, `Adjacency` and `RelativeAllele`.

use std::cmp::Ordering;

use super::{NormalizeOptions, SequenceProvider, normalize_allele};
use crate::digest::Identifiable;
use crate::error::NormalizeError;
use crate::model::{
    Adjacency, Allele, CisPhasedBlock, IriOr, Location, RelativeAllele, RelativeSequenceLocation,
    SequenceExpression, SequenceLocation, SequenceOffsetLocation,
};

/// Normalize a [`CisPhasedBlock`]: every inline member allele is normalized (using the block's
/// shared `sequenceReference` for members that omit their own), then members are ordered by
/// digest so that the JSON form is canonical (the digest serialization sorts them regardless).
///
/// # Errors
/// See [`normalize_allele`].
pub fn normalize_cis_phased_block<P: SequenceProvider + ?Sized>(
    block: &CisPhasedBlock,
    provider: &P,
    options: &NormalizeOptions,
) -> Result<CisPhasedBlock, NormalizeError> {
    let mut out = block.clone();
    for member in out.members_mut().iter_mut() {
        if let IriOr::Object(allele) = member {
            *allele = normalize_member(allele, block, provider, options)?;
        }
    }
    let mut keyed: Vec<(Vec<u8>, IriOr<Allele>)> = out
        .members_mut()
        .drain(..)
        .map(|m| {
            let key = match &m {
                IriOr::Object(a) => a.digest().as_bytes().to_vec(),
                IriOr::Iri(iri) => iri
                    .ga4gh_digest()
                    .unwrap_or(iri.as_str())
                    .as_bytes()
                    .to_vec(),
            };
            (key, m)
        })
        .collect();
    keyed.sort_by(|a, b| a.0.cmp(&b.0));
    out.members_mut().extend(keyed.into_iter().map(|(_, m)| m));
    Ok(out)
}

/// Normalize a member allele, borrowing the block's sequence reference when the member's
/// location has none; the output location keeps the member's original (absent) reference.
fn normalize_member<P: SequenceProvider + ?Sized>(
    allele: &Allele,
    block: &CisPhasedBlock,
    provider: &P,
    options: &NormalizeOptions,
) -> Result<Allele, NormalizeError> {
    let needs_shared = allele
        .sequence_location()
        .is_some_and(|l| l.sequence_reference().is_none())
        && block.sequence_reference().is_some();
    if !needs_shared {
        return normalize_allele(allele, provider, options);
    }
    let location = allele.sequence_location().expect("checked above");
    let shared = block.sequence_reference().expect("checked above").clone();
    let with_ref =
        SequenceLocation::from_parts(Some(shared.into()), location.start(), location.end())?;
    let tmp = Allele::new(with_ref, allele.state().clone());
    let normalized = normalize_allele(&tmp, provider, options)?;
    let (loc, state) = normalized.into_parts();
    let loc = loc
        .into_object()
        .expect("normalize_allele keeps inline locations");
    let stripped = SequenceLocation::from_parts(None, loc.start(), loc.end())?;
    Ok(allele.rebuilt(stripped, state))
}

/// Apply the *Adjacency Normalization* orientation convention (VRS 2.1.0):
///
/// 1. the first adjoined sequence SHOULD have forward orientation (a location defined by
///    `end`);
/// 2. the adjoined sequence accessions are equal or in ascending lexicographical order;
/// 3. the defined coordinates are in ascending numerical order.
///
/// The two orientations of an adjacency are the two orders of its adjoined sequences; the
/// order that satisfies the criteria (compared in sequence) is returned. Accessions are
/// compared by RefGet accession, the only sequence identifier that participates in digests.
///
/// Orientation is only changed when both adjoined sequences are inline `SequenceLocation`s
/// with inline sequence references and the linker (if any) is orientation-free
/// (`LengthExpression` or absent): reversing a literal linker would require reverse
/// complementation, which the specification does not describe. Ambiguity expansion for
/// homologous breakpoints is not yet specified upstream and is therefore not implemented.
pub fn normalize_adjacency(adjacency: &Adjacency) -> Adjacency {
    let [a, b] = adjacency.adjoined_sequences();
    let (Some(a_loc), Some(b_loc)) = (sequence_location(a), sequence_location(b)) else {
        return adjacency.clone();
    };
    if matches!(
        adjacency.linker(),
        Some(SequenceExpression::Literal(_) | SequenceExpression::ReferenceLength(_))
    ) {
        return adjacency.clone();
    }
    match orientation_key(a_loc, b_loc).cmp(&orientation_key(b_loc, a_loc)) {
        Ordering::Less => adjacency.with_adjoined([b.clone(), a.clone()]),
        _ => adjacency.clone(),
    }
}

fn sequence_location(v: &IriOr<Location>) -> Option<&SequenceLocation> {
    match v {
        IriOr::Object(Location::Sequence(l)) => Some(l),
        _ => None,
    }
}

/// Higher is better: (first is defined by `end`, accessions ascending, coordinates ascending).
fn orientation_key(first: &SequenceLocation, second: &SequenceLocation) -> (bool, bool, bool) {
    let forward = first.end().is_some();
    let acc = |l: &SequenceLocation| l.refget_accession().map(|a| a.as_str().to_owned());
    let accessions_ascending = match (acc(first), acc(second)) {
        (Some(a), Some(b)) => a <= b,
        _ => true,
    };
    let coord = |l: &SequenceLocation| l.end().or(l.start()).and_then(|c| c.lower_bound());
    let coordinates_ascending = match (coord(first), coord(second)) {
        (Some(a), Some(b)) => a <= b,
        _ => true,
    };
    (forward, accessions_ascending, coordinates_ascending)
}

/// Which anchor a `RelativeAllele` should persist.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum AnchorChoice {
    /// The left-anchored representation.
    Left,
    /// The right-anchored representation.
    Right,
}

/// *Relative Allele Normalization* anchor selection (VRS 2.1.0, draft): given the two
/// equivalent mapped locations of a relative allele — one expressed from the left anchor of
/// the gap and one from the right anchor — select the one whose largest offset magnitude is
/// smaller; on a tie, the left anchor.
///
/// Deriving the alternative representation requires the transcript alignment (the position
/// of the opposite exon boundary), which a `RelativeAllele` does not carry, so the caller
/// computes both candidates and this function applies the rule.
pub fn preferred_anchor(
    left: &SequenceOffsetLocation,
    right: &SequenceOffsetLocation,
) -> AnchorChoice {
    let l = left.max_offset_magnitude().unwrap_or(0);
    let r = right.max_offset_magnitude().unwrap_or(0);
    if r < l {
        AnchorChoice::Right
    } else {
        AnchorChoice::Left
    }
}

/// Normalize a [`RelativeAllele`] (draft): the base representation (`baseState` on the
/// `baseSequenceLocation`) is normalized exactly like an [`Allele`].
///
/// If normalization moves the base location (an ambiguous indel), the mapped offsets would
/// have to be re-derived from the transcript alignment, which the object does not carry; in
/// that case the input is returned unchanged. State-only changes (identity → reference-length
/// encoding, trimmed substitutions at the same span) are applied. Use [`preferred_anchor`] to
/// choose between anchor representations when the alignment is available.
///
/// # Errors
/// See [`normalize_allele`].
pub fn normalize_relative_allele<P: SequenceProvider + ?Sized>(
    relative: &RelativeAllele,
    provider: &P,
    options: &NormalizeOptions,
) -> Result<RelativeAllele, NormalizeError> {
    let IriOr::Object(rsl) = relative.relative_location() else {
        return Ok(relative.clone());
    };
    let IriOr::Object(base_loc) = rsl.base_sequence_location() else {
        return Ok(relative.clone());
    };
    let base = Allele::new(base_loc.clone(), relative.base_state().clone());
    let normalized = normalize_allele(&base, provider, options)?;
    let (loc, state) = normalized.into_parts();
    let Some(loc) = loc.into_object() else {
        return Ok(relative.clone());
    };
    if loc.start() != base_loc.start() || loc.end() != base_loc.end() {
        return Ok(relative.clone());
    }
    let mut new_rsl = RelativeSequenceLocation::new(loc, rsl.mapped_sequence_location().clone());
    if let Some(m) = crate::model::Entity::meta(rsl) {
        new_rsl = new_rsl.with_meta(m.clone());
    }
    let mut out = RelativeAllele::new(new_rsl, state, relative.mapped_state().clone());
    if let Some(m) = crate::model::Entity::meta(relative) {
        out = out.with_meta(m.clone());
    }
    Ok(out)
}
