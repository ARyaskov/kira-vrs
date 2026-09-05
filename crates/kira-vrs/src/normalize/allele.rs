//! Allele normalization (VRS 2.1.0 *Allele Normalization*, "LiteralSequenceExpression
//! Alleles").
//!
//! The algorithm is the fully-justified normalization adapted from NCBI's Variant
//! Overprecision Correction Algorithm (VOCA), extended in VRS 2 with reference-length
//! encoding and revised in VRS 2.1 to select the *smallest* repeat-subunit factor. Step
//! numbers below refer to the numbered list in the specification.
//!
//! ```text
//! 0. ref = reference[start, end), alt = state.sequence
//! 1. trim common suffix, then common prefix                     (O(n))
//! 2. both empty      → reference allele: RLE(len, len) at original location
//!    both non-empty  → substitution: Literal(alt) at trimmed location
//!    one empty       → seed = the non-empty one; continue
//! 3. roll left / roll right over the region of ambiguity          (O(width of ambiguity))
//! 4. expand: prepend / append reference to both sequences
//! 5. width 0         → unambiguous insertion: Literal(alt)
//!    deletion        → RLE(length = |alt|, repeatSubunitLength = |seed|)
//!    insertion       → smallest factor d of |seed| (d ≤ |ref|) such that alt is the
//!                      d-periodic extension of ref → RLE(|alt|, d); else Literal(alt)
//! ```
//!
//! Reference bases are fetched in windows (the interval plus context) and the window is
//! extended lazily if a roll reaches its edge, so a normalization needs one length query and
//! typically one sequence fetch. Coordinates given as indefinite ranges (`[null, x]` /
//! `[x, null]`) are normalized on their defined bound and keep their form; alleles with
//! definite ranges, non-literal states, IRI locations or unknown sequence references are
//! returned unchanged, as are alleles on circular references (which the algorithm does not
//! define).

use std::borrow::Cow;

use super::{NormalizeOptions, SequenceProvider};
use crate::error::{CoordinateError, NormalizeError, SequenceError, UnsupportedError};
use crate::model::{
    Allele, IntOrRange, IriOr, LiteralSequenceExpression, Range, ReferenceLengthExpression,
    RefgetAccession, SequenceExpression, SequenceLocation, SequenceReference, SequenceString,
};

/// How a coordinate was expressed in the input, so the output keeps the same form.
#[derive(Clone, Copy)]
enum PosKind {
    Exact,
    AtMost,
    AtLeast,
}

fn pos_kind(v: IntOrRange) -> Option<(u64, PosKind)> {
    match v {
        IntOrRange::Int(i) => Some((i as u64, PosKind::Exact)),
        IntOrRange::Range(r) => match (r.min(), r.max()) {
            (Some(_), Some(_)) | (None, None) => None,
            (None, Some(hi)) => Some((hi as u64, PosKind::AtMost)),
            (Some(lo), None) => Some((lo as u64, PosKind::AtLeast)),
        },
    }
}

fn rebuild_pos(value: u64, kind: PosKind) -> Result<IntOrRange, CoordinateError> {
    let v = i64::try_from(value).map_err(|_| CoordinateError::OutOfRange(i128::from(value)))?;
    Ok(match kind {
        PosKind::Exact => IntOrRange::Int(v),
        PosKind::AtMost => IntOrRange::Range(Range::at_most(v)?),
        PosKind::AtLeast => IntOrRange::Range(Range::at_least(v)?),
    })
}

/// The RefGet accession behind an allele's location: an inline `SequenceReference`, or an
/// IRI of the form `ga4gh:SQ.<digest>`.
fn location_accession(location: &SequenceLocation) -> Option<RefgetAccession> {
    match location.sequence_reference()? {
        IriOr::Object(r) => Some(*r.refget_accession()),
        IriOr::Iri(iri) => {
            let rest = iri.as_str().strip_prefix("ga4gh:")?;
            RefgetAccession::parse(rest).ok()
        }
    }
}

/// A window of reference sequence around the allele, extended on demand.
struct Window<'p, P: ?Sized> {
    provider: &'p P,
    accession: RefgetAccession,
    total_length: u64,
    offset: u64,
    bytes: Vec<u8>,
}

impl<P: SequenceProvider + ?Sized> Window<'_, P> {
    fn end(&self) -> u64 {
        self.offset + self.bytes.len() as u64
    }

    fn fetch(
        provider: &P,
        accession: &RefgetAccession,
        start: u64,
        end: u64,
    ) -> Result<Vec<u8>, SequenceError> {
        let bytes = provider.sequence(accession, start, end)?;
        crate::model::primitives::validate_residues(&bytes)?;
        Ok(match bytes {
            Cow::Borrowed(b) => b.to_vec(),
            Cow::Owned(v) => v,
        })
    }

    /// Fetch `[start, end)`, or exactly `[need_start, need_end)` if the provider cannot serve
    /// the wider window (providers backed by partial segments).
    fn fetch_window(
        provider: &P,
        accession: &RefgetAccession,
        start: u64,
        end: u64,
        need_start: u64,
        need_end: u64,
    ) -> Result<(u64, Vec<u8>), SequenceError> {
        match Self::fetch(provider, accession, start, end) {
            Ok(bytes) => Ok((start, bytes)),
            Err(SequenceError::OutOfBounds { .. }) if (start, end) != (need_start, need_end) => {
                Self::fetch(provider, accession, need_start, need_end).map(|b| (need_start, b))
            }
            Err(e) => Err(e),
        }
    }

    /// Ensure `[lo, hi)` is loaded, growing the window geometrically.
    fn ensure(&mut self, lo: u64, hi: u64) -> Result<(), SequenceError> {
        if lo >= self.offset && hi <= self.end() {
            return Ok(());
        }
        let grow = (self.bytes.len() as u64).max(64);
        if lo < self.offset {
            let want = self.offset.saturating_sub(grow).min(lo);
            let (got_lo, mut head) = Self::fetch_window(
                self.provider,
                &self.accession,
                want,
                self.offset,
                lo,
                self.offset,
            )?;
            head.extend_from_slice(&self.bytes);
            self.bytes = head;
            self.offset = got_lo;
        }
        if hi > self.end() {
            let want = (self.end() + grow).max(hi).min(self.total_length);
            let (_, tail) = Self::fetch_window(
                self.provider,
                &self.accession,
                self.end(),
                want,
                self.end(),
                hi,
            )?;
            self.bytes.extend_from_slice(&tail);
        }
        Ok(())
    }

    #[inline]
    fn slice(&self, lo: u64, hi: u64) -> &[u8] {
        &self.bytes[(lo - self.offset) as usize..(hi - self.offset) as usize]
    }

    #[inline]
    fn byte(&self, pos: u64) -> u8 {
        self.bytes[(pos - self.offset) as usize]
    }
}

/// Normalize an [`Allele`] using the fully-justified algorithm (see module docs).
///
/// # Errors
/// [`NormalizeError::Sequence`] if reference sequence cannot be fetched;
/// [`NormalizeError::Unsupported`] for circular references or `start > end`.
pub fn normalize_allele<P: SequenceProvider + ?Sized>(
    allele: &Allele,
    provider: &P,
    options: &NormalizeOptions,
) -> Result<Allele, NormalizeError> {
    // Only literal states are normalized; RLE and length states are already canonical or
    // carry no sequence to normalize (specification: other types are returned as-is).
    let SequenceExpression::Literal(literal) = allele.state() else {
        return Ok(allele.clone());
    };
    let Some(location) = allele.sequence_location() else {
        return Ok(allele.clone());
    };
    let Some(accession) = location_accession(location) else {
        return Ok(allele.clone());
    };
    if location
        .inline_sequence_reference()
        .is_some_and(SequenceReference::is_circular)
    {
        return Err(UnsupportedError::new("normalization of alleles on circular sequences").into());
    }
    let (Some(start), Some(end)) = (location.start(), location.end()) else {
        return Ok(allele.clone());
    };
    let (Some((start, start_kind)), Some((end, end_kind))) = (pos_kind(start), pos_kind(end))
    else {
        return Ok(allele.clone()); // definite ranges are not normalized
    };
    if start > end {
        return Err(UnsupportedError::new("normalization of alleles with start > end").into());
    }

    let total_length = provider.sequence_length(&accession)?;
    if end > total_length {
        return Err(SequenceError::OutOfBounds {
            accession: accession.to_string(),
            start,
            end,
        }
        .into());
    }
    let alt: &[u8] = literal.sequence().as_bytes();

    // Step 0: fetch the interval plus context for rolling.
    let context = (alt.len() as u64).max(end - start).max(16) * 2;
    let win_lo = start.saturating_sub(context);
    let win_hi = (end + context).min(total_length);
    let (offset, bytes) =
        Window::<P>::fetch_window(provider, &accession, win_lo, win_hi, start, end)?;
    let mut window = Window {
        provider,
        accession,
        total_length,
        offset,
        bytes,
    };

    // Step 1: trim common suffix, then common prefix.
    let reference = window.slice(start, end);
    let suffix = common_suffix(reference, alt);
    let (reference, alt) = (
        &reference[..reference.len() - suffix],
        &alt[..alt.len() - suffix],
    );
    let prefix = common_prefix(reference, alt);
    let trimmed_ref = &reference[prefix..];
    let trimmed_alt = &alt[prefix..];
    let t_start = start + prefix as u64;
    let t_end = end - suffix as u64;

    // Step 2.
    if trimmed_ref.is_empty() && trimmed_alt.is_empty() {
        // Reference allele: RLE over the original location.
        let length = end - start;
        let seq = window.slice(start, end);
        let state = rle_state(length, length, seq, options)?;
        let loc = location.with_coordinates(
            Some(rebuild_pos(start, start_kind)?),
            Some(rebuild_pos(end, end_kind)?),
        )?;
        return Ok(allele.rebuilt(loc, state));
    }
    if !trimmed_ref.is_empty() && !trimmed_alt.is_empty() {
        // Substitution.
        let state = LiteralSequenceExpression::new(SequenceString::from_valid_bytes(trimmed_alt));
        let loc = location.with_coordinates(
            Some(rebuild_pos(t_start, start_kind)?),
            Some(rebuild_pos(t_end, end_kind)?),
        )?;
        return Ok(allele.rebuilt(loc, state.into()));
    }

    // Step 3: roll. `seed` is the non-empty sequence; for a deletion it is reference-derived.
    let is_deletion = trimmed_alt.is_empty();
    let seed: Vec<u8> = if is_deletion {
        trimmed_ref.to_vec()
    } else {
        trimmed_alt.to_vec()
    };
    let seed_len = seed.len() as u64;

    let left_bound = roll_left(&mut window, &seed, t_start)?;
    let right_bound = roll_right(&mut window, &seed, t_end)?;

    // Step 4: expand to the region of ambiguity.
    window.ensure(left_bound, right_bound)?;
    let expanded_ref = window.slice(left_bound, right_bound);
    let mut expanded_alt = Vec::with_capacity(expanded_ref.len() + trimmed_alt.len());
    expanded_alt.extend_from_slice(window.slice(left_bound, t_start));
    expanded_alt.extend_from_slice(trimmed_alt);
    expanded_alt.extend_from_slice(window.slice(t_end, right_bound));

    let loc = location.with_coordinates(
        Some(rebuild_pos(left_bound, start_kind)?),
        Some(rebuild_pos(right_bound, end_kind)?),
    )?;

    // Step 5.
    let state: SequenceExpression = if expanded_ref.is_empty() {
        // Unambiguous insertion.
        LiteralSequenceExpression::new(SequenceString::from_valid_bytes(&expanded_alt)).into()
    } else if is_deletion {
        rle_state(expanded_alt.len() as u64, seed_len, &expanded_alt, options)?
    } else {
        match smallest_reference_period(seed_len, expanded_ref, &expanded_alt) {
            Some(d) => rle_state(expanded_alt.len() as u64, d, &expanded_alt, options)?,
            None => LiteralSequenceExpression::new(SequenceString::from_valid_bytes(&expanded_alt))
                .into(),
        }
    };
    Ok(allele.rebuilt(loc, state))
}

fn rle_state(
    length: u64,
    repeat_subunit_length: u64,
    sequence: &[u8],
    options: &NormalizeOptions,
) -> Result<SequenceExpression, NormalizeError> {
    let length_i =
        i64::try_from(length).map_err(|_| CoordinateError::OutOfRange(i128::from(length)))?;
    let rsl = i64::try_from(repeat_subunit_length)
        .map_err(|_| CoordinateError::OutOfRange(i128::from(repeat_subunit_length)))?;
    let mut rle = ReferenceLengthExpression::new(length_i, rsl)?;
    if options.attach_rle_sequence(length) {
        rle = rle.with_sequence(SequenceString::from_valid_bytes(sequence));
    }
    Ok(rle.into())
}

/// Length of the common prefix of two byte strings.
#[inline]
fn common_prefix(a: &[u8], b: &[u8]) -> usize {
    a.iter().zip(b).take_while(|(x, y)| x == y).count()
}

/// Length of the common suffix of two byte strings.
#[inline]
fn common_suffix(a: &[u8], b: &[u8]) -> usize {
    a.iter()
        .rev()
        .zip(b.iter().rev())
        .take_while(|(x, y)| x == y)
        .count()
}

/// Step 3a: roll the seed leftwards while the base preceding the window equals its last base.
/// Returns the left roll bound.
fn roll_left<P: SequenceProvider + ?Sized>(
    window: &mut Window<'_, P>,
    seed: &[u8],
    start: u64,
) -> Result<u64, SequenceError> {
    let mut bound = start;
    // Rotation index: the seed rotated right by k has last byte seed[(n - 1 - k) mod n].
    let n = seed.len();
    let mut k = 0usize;
    while bound > 0 {
        if bound - 1 < window.offset {
            window.ensure(bound - 1, bound)?;
        }
        let last = seed[(n - 1 + n - (k % n)) % n];
        if window.byte(bound - 1) != last {
            break;
        }
        bound -= 1;
        k += 1;
    }
    Ok(bound)
}

/// Step 3b: roll the seed rightwards while the base following the window equals its first
/// base. Returns the right roll bound. `end` is the end of the trimmed interval (equal to
/// `start` for an insertion, `start + |seed|` for a deletion).
fn roll_right<P: SequenceProvider + ?Sized>(
    window: &mut Window<'_, P>,
    seed: &[u8],
    end: u64,
) -> Result<u64, SequenceError> {
    let mut bound = end;
    let n = seed.len();
    let mut k = 0usize;
    while bound < window.total_length {
        if bound >= window.end() {
            window.ensure(bound, bound + 1)?;
        }
        let first = seed[k % n];
        if window.byte(bound) != first {
            break;
        }
        bound += 1;
        k += 1;
    }
    Ok(bound)
}

/// Step 5c: the smallest factor `d` of `seed_len` with `d <= |reference|` such that
/// `alt` equals the cyclic extension of `reference[..d]` truncated to `|alt|` — i.e. such that
/// `ReferenceLengthExpression { length: |alt|, repeatSubunitLength: d }` expands to `alt`.
fn smallest_reference_period(seed_len: u64, reference: &[u8], alt: &[u8]) -> Option<u64> {
    let ref_len = reference.len() as u64;
    factors_ascending(seed_len)
        .filter(|&d| d <= ref_len)
        .find(|&d| is_periodic_extension(&reference[..d as usize], alt))
}

/// `alt` equals `unit` repeated cyclically and truncated to `|alt|`.
fn is_periodic_extension(unit: &[u8], alt: &[u8]) -> bool {
    !unit.is_empty() && alt.iter().zip(unit.iter().cycle()).all(|(a, u)| a == u)
}

/// Factors of `n` in ascending order.
fn factors_ascending(n: u64) -> impl Iterator<Item = u64> {
    let mut small = Vec::new();
    let mut large = Vec::new();
    let mut i = 1;
    while i * i <= n {
        if n.is_multiple_of(i) {
            small.push(i);
            if i != n / i {
                large.push(n / i);
            }
        }
        i += 1;
    }
    small.into_iter().chain(large.into_iter().rev())
}

/// Expand a `ReferenceLengthExpression` into the literal sequence it denotes, given the
/// reference sequence at the allele's location.
///
/// The subunit is the first `repeatSubunitLength` residues of `reference`; it is repeated
/// (with a trailing partial copy if needed) until `length` residues are produced. This is the
/// inverse of the encoding performed by [`normalize_allele`].
///
/// # Errors
/// [`CoordinateError::OutOfRange`] if the reference is shorter than the repeat subunit or the
/// expression's length is not an exact integer.
pub fn expand_reference_length_expression(
    reference: &[u8],
    expression: &ReferenceLengthExpression,
) -> Result<SequenceString, CoordinateError> {
    let length = expression
        .length()
        .as_int()
        .ok_or(CoordinateError::OutOfRange(-1))? as usize;
    let d = expression.repeat_subunit_length() as usize;
    if d > reference.len() || (d == 0 && length > 0) {
        return Err(CoordinateError::OutOfRange(d as i128));
    }
    crate::model::primitives::validate_residues(reference)
        .map_err(|_| CoordinateError::OutOfRange(-1))?;
    let unit = &reference[..d];
    let bytes: Vec<u8> = unit.iter().copied().cycle().take(length).collect();
    Ok(SequenceString::from_valid_bytes(&bytes))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn factors() {
        assert_eq!(
            factors_ascending(12).collect::<Vec<_>>(),
            vec![1, 2, 3, 4, 6, 12]
        );
        assert_eq!(factors_ascending(1).collect::<Vec<_>>(), vec![1]);
        assert_eq!(factors_ascending(9).collect::<Vec<_>>(), vec![1, 3, 9]);
    }

    #[test]
    fn periodic_extension() {
        assert!(is_periodic_extension(b"CAG", b"CAGCAGCAGCA"));
        assert!(!is_periodic_extension(b"CA", b"CAGCAGCAGCA"));
        assert!(is_periodic_extension(b"A", b""));
    }

    #[test]
    fn trims() {
        assert_eq!(common_prefix(b"ACGT", b"ACTT"), 2);
        assert_eq!(common_suffix(b"ACGT", b"TT"), 1);
        assert_eq!(common_suffix(b"", b"TT"), 0);
    }
}
