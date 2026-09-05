//! Allele normalization against the cases the reference implementation (vrs-python) tests,
//! using the exact reference bases its recorded SeqRepo responses contain, plus VRS 2.1
//! specific cases.
//!
//! Each case builds an in-memory provider holding a partial segment of the real chromosome,
//! padded with `N` (which never matches a real base, so it cannot extend a roll).

use kira_vrs::digest::Identifiable;
use kira_vrs::model::*;
use kira_vrs::normalize::{
    InMemorySequenceProvider, NormalizeOptions, SequenceProvider,
    expand_reference_length_expression, normalize_allele, normalize_cis_phased_block,
};

const CHR1: &str = "SQ.Ya6Rs7DHhDeg7YaOSg1EoNi3U_nQ9SvO";
const CHR6: &str = "SQ.0iKlIQk2oZLoeOG9P1riRU6hvL5Ux8TV";
const CHR19: &str = "SQ.IIB53T8CNeJJdUqzn9V_JnRtQadwWCbl";
const CHRX: &str = "SQ.w0WZEvgJF0zf_P4yyTzjjv9oW1z61HHP";
const PAD: u64 = 64;

/// A provider holding `bases` at `offset` on `accession`, padded with `N` on both sides.
fn segment(
    accession: &str,
    offset: u64,
    bases: &str,
) -> (InMemorySequenceProvider, RefgetAccession) {
    let acc = RefgetAccession::parse(accession).unwrap();
    let mut seq = "N".repeat(PAD as usize);
    seq.push_str(bases);
    seq.push_str(&"N".repeat(PAD as usize));
    let mut p = InMemorySequenceProvider::new();
    p.insert_segment(acc, offset - PAD, seq.as_bytes(), 250_000_000)
        .unwrap();
    (p, acc)
}

fn allele(
    acc: RefgetAccession,
    start: impl Into<IntOrRange>,
    end: impl Into<IntOrRange>,
    alt: &str,
) -> Allele {
    let loc = SequenceLocation::new(SequenceReference::new(acc), start, end).unwrap();
    Allele::new(loc, SequenceExpression::literal(alt).unwrap())
}

fn span(a: &Allele) -> (IntOrRange, IntOrRange) {
    let l = a.sequence_location().unwrap();
    (l.start().unwrap(), l.end().unwrap())
}

fn rle(a: &Allele) -> (i64, i64, Option<String>) {
    let r = a
        .state()
        .as_reference_length()
        .expect("expected ReferenceLengthExpression");
    (
        r.length().as_int().unwrap(),
        r.repeat_subunit_length(),
        r.sequence().map(ToString::to_string),
    )
}

fn literal(a: &Allele) -> String {
    a.state()
        .as_literal()
        .expect("expected LiteralSequenceExpression")
        .sequence()
        .to_string()
}

#[test]
fn reference_allele_snv_becomes_rle() {
    let (p, acc) = segment(CHR6, 26_090_950, "C");
    let a = allele(acc, 26_090_950, 26_090_951, "C");
    let n = normalize_allele(&a, &p, &NormalizeOptions::default()).unwrap();
    assert_eq!(
        span(&n),
        (IntOrRange::Int(26_090_950), IntOrRange::Int(26_090_951))
    );
    assert_eq!(rle(&n), (1, 1, Some("C".into())));
}

#[test]
fn indefinite_range_deletion_becomes_rle_and_keeps_range_form() {
    // 155980374 T, 155980375 T, 155980376 A, 155980377 A
    let (p, acc) = segment(CHRX, 155_980_374, "TTAA");
    let loc = SequenceLocation::new(
        SequenceReference::new(acc),
        Range::at_most(155_980_375).unwrap(),
        Range::at_least(155_980_377).unwrap(),
    )
    .unwrap();
    let a = Allele::new(loc, SequenceExpression::literal("").unwrap());
    let n = normalize_allele(&a, &p, &NormalizeOptions::default().without_rle_sequence()).unwrap();
    assert_eq!(
        span(&n),
        (
            Range::at_most(155_980_375).unwrap().into(),
            Range::at_least(155_980_377).unwrap().into()
        )
    );
    assert_eq!(rle(&n), (0, 2, None));
}

#[test]
fn definite_ranges_are_left_untouched() {
    let (p, acc) = segment(CHRX, 155_980_374, "TTAA");
    let loc = SequenceLocation::new(
        SequenceReference::new(acc),
        Range::bounded(155_980_374, 155_980_375).unwrap(),
        Range::bounded(155_980_377, 155_980_378).unwrap(),
    )
    .unwrap();
    let a = Allele::new(loc, SequenceExpression::literal("").unwrap());
    let n = normalize_allele(&a, &p, &NormalizeOptions::default()).unwrap();
    assert_eq!(n, a);
}

#[test]
fn tandem_duplication_becomes_rle() {
    // 155980372 G, 373 G, 374 T, 375 T
    let (p, acc) = segment(CHRX, 155_980_372, "GGTT");
    let a = allele(acc, 155_980_373, 155_980_375, "GTGT");
    let n = normalize_allele(&a, &p, &NormalizeOptions::default()).unwrap();
    assert_eq!(span(&n), (155_980_373.into(), 155_980_375.into()));
    assert_eq!(rle(&n), (4, 2, Some("GTGT".into())));
}

#[test]
fn repeat_subunit_insertion_in_non_integer_repeat_region() {
    // 289463 T, 464 C, 465 A, 466 G, 467 C, 468 A, 469 C
    let (p, acc) = segment(CHR19, 289_463, "TCAGCAC");
    let a = allele(acc, 289_464, 289_464, "CAGCAG");
    let n = normalize_allele(&a, &p, &NormalizeOptions::default()).unwrap();
    assert_eq!(span(&n), (289_464.into(), 289_469.into()));
    assert_eq!(rle(&n), (11, 3, Some("CAGCAGCAGCA".into())));
}

#[test]
fn reference_allele_multi_base() {
    let (p, acc) = segment(CHR1, 100_210_777, "AA");
    let a = allele(acc, 100_210_777, 100_210_779, "AA");
    let n = normalize_allele(&a, &p, &NormalizeOptions::default()).unwrap();
    assert_eq!(rle(&n), (2, 2, Some("AA".into())));
    assert_eq!(span(&n), (100_210_777.into(), 100_210_779.into()));
}

#[test]
fn multi_base_substitution_stays_literal() {
    let (p, acc) = segment(CHR1, 939_145, "GA");
    let a = allele(acc, 939_145, 939_147, "TT");
    let n = normalize_allele(&a, &p, &NormalizeOptions::default()).unwrap();
    assert_eq!(n, a);
}

#[test]
fn unambiguous_insertion_in_repeat_stays_literal() {
    let (p, acc) = segment(CHR1, 236_900_412, "AA");
    let a = allele(acc, 236_900_413, 236_900_413, "CGT");
    let n = normalize_allele(&a, &p, &NormalizeOptions::default()).unwrap();
    assert_eq!(n, a);
}

#[test]
fn simple_deletion_from_non_repeat_region() {
    // 66925 A, 66926 G, 66927 A
    let (p, acc) = segment(CHR1, 66_925, "AGA");
    let a = allele(acc, 66_926, 66_927, "");
    let n = normalize_allele(&a, &p, &NormalizeOptions::default().without_rle_sequence()).unwrap();
    assert_eq!(span(&n), (66_926.into(), 66_927.into()));
    assert_eq!(rle(&n), (0, 1, None));
}

#[test]
fn microsatellite_deletion_expands_to_full_repeat() {
    // 766398 G, 399 A, 400 A, 401 T, 402 A, 403 A, 404 A, 405 T, 406 A, 407 C
    let (p, acc) = segment(CHR1, 766_398, "GAATAAATAC");
    let a = allele(acc, 766_400, 766_404, "");
    let n = normalize_allele(&a, &p, &NormalizeOptions::default().without_rle_sequence()).unwrap();
    assert_eq!(span(&n), (766_399.into(), 766_407.into()));
    assert_eq!(rle(&n), (4, 4, None));
}

#[test]
fn tandem_repeat_deletion() {
    // 930135 T, 136 C, 137 T, 138 C, 139 C, 140 T, 141 C, 142 C, 143 T, 144 G
    let (p, acc) = segment(CHR1, 930_135, "TCTCCTCCTG");
    let a = allele(acc, 930_137, 930_140, "");
    let n = normalize_allele(&a, &p, &NormalizeOptions::default().without_rle_sequence()).unwrap();
    assert_eq!(span(&n), (930_136.into(), 930_144.into()));
    assert_eq!(rle(&n), (5, 3, None));
}

/// chr1:1752907-1752937 = C + CCT×9 + C + G (the CCT microsatellite used by many upstream cases).
fn cct_region() -> (InMemorySequenceProvider, RefgetAccession) {
    segment(CHR1, 1_752_907, "CCCTCCTCCTCCTCCTCCTCCTCCTCCTCG")
}

#[test]
fn repeat_subunit_insertion_in_microsatellite() {
    let (p, acc) = cct_region();
    let a = allele(acc, 1_752_908, 1_752_908, "CCT");
    let n = normalize_allele(&a, &p, &NormalizeOptions::default().without_rle_sequence()).unwrap();
    assert_eq!(span(&n), (1_752_908.into(), 1_752_936.into()));
    assert_eq!(rle(&n), (31, 3, None));
}

#[test]
fn partial_repeat_insertions_and_deletions() {
    let (p, acc) = cct_region();
    let opts = NormalizeOptions::default().without_rle_sequence();
    let n = |start: i64, end: i64, alt: &str| {
        normalize_allele(&allele(acc, start, end, alt), &p, &opts).unwrap()
    };

    // middle insertion, 2 bp
    let x = n(1_752_915, 1_752_915, "CT");
    assert_eq!(span(&x), (1_752_915.into(), 1_752_918.into()));
    assert_eq!(rle(&x), (5, 2, None));
    // middle insertion, 4 bp
    let x = n(1_752_915, 1_752_915, "CCTC");
    assert_eq!(span(&x), (1_752_911.into(), 1_752_916.into()));
    assert_eq!(rle(&x), (9, 4, None));
    // tail insertion, 2 bp: unambiguous, stays literal
    let x = n(1_752_934, 1_752_934, "CT");
    assert_eq!(span(&x), (1_752_934.into(), 1_752_934.into()));
    assert_eq!(literal(&x), "CT");
    // tail insertion, 4 bp: rolls left by one, not reference-derived
    let x = n(1_752_934, 1_752_934, "CCTC");
    assert_eq!(span(&x), (1_752_933.into(), 1_752_934.into()));
    assert_eq!(literal(&x), "CCCTC");
    // middle deletion, 4 bp
    let x = n(1_752_912, 1_752_916, "");
    assert_eq!(span(&x), (1_752_911.into(), 1_752_916.into()));
    assert_eq!(rle(&x), (1, 4, None));
    // middle deletion, 2 bp
    let x = n(1_752_913, 1_752_915, "");
    assert_eq!(span(&x), (1_752_912.into(), 1_752_915.into()));
    assert_eq!(rle(&x), (1, 2, None));
    // tail deletion, 2 bp
    let x = n(1_752_934, 1_752_936, "");
    assert_eq!(span(&x), (1_752_933.into(), 1_752_936.into()));
    assert_eq!(rle(&x), (1, 2, None));
    // tail deletion, 4 bp
    let x = n(1_752_932, 1_752_936, "");
    assert_eq!(span(&x), (1_752_932.into(), 1_752_936.into()));
    assert_eq!(rle(&x), (0, 4, None));
}

/// VRS 2.1 selects the *smallest* repeat-subunit factor. Inserting `CACA` into a `CACACA`
/// run is reference-derived with period 2 (and also with period 4); 2.1 chooses 2. (VRS 2.0
/// and vrs-python ≤ 2.4.0-a4 choose the greatest factor, 4.)
#[test]
fn smallest_repeat_subunit_factor_per_vrs_2_1() {
    let mut p = InMemorySequenceProvider::new();
    let acc = p.insert_sequence(b"GGCACACATT").unwrap();
    let a = allele(acc, 2, 2, "CACA");
    let n = normalize_allele(&a, &p, &NormalizeOptions::default()).unwrap();
    assert_eq!(span(&n), (2.into(), 8.into()));
    assert_eq!(rle(&n), (10, 2, Some("CACACACACA".into())));
}

#[test]
fn equivalent_representations_share_an_identifier() {
    let mut p = InMemorySequenceProvider::new();
    let acc = p.insert_sequence(b"GGTATATATACC").unwrap();
    let opts = NormalizeOptions::default();
    // Insert "TA" anywhere in the TA repeat, or as "AT" one base to the right: same variant.
    let reps = [
        allele(acc, 2, 2, "TA"),
        allele(acc, 4, 4, "TA"),
        allele(acc, 3, 3, "AT"),
        allele(acc, 10, 10, "TA"),
    ];
    let ids: Vec<_> = reps
        .iter()
        .map(|a| normalize_allele(a, &p, &opts).unwrap().identifier())
        .collect();
    assert!(ids.iter().all(|i| i == &ids[0]), "{ids:?}");
    // And the VCF-style representation with an anchor base ("GTA" replacing "G" at 1).
    let vcf = allele(acc, 1, 2, "GTA");
    assert_eq!(
        normalize_allele(&vcf, &p, &opts).unwrap().identifier(),
        ids[0]
    );
}

#[test]
fn normalization_is_idempotent() {
    let (p, acc) = cct_region();
    let opts = NormalizeOptions::default();
    let a = allele(acc, 1_752_915, 1_752_915, "CCTC");
    let n1 = normalize_allele(&a, &p, &opts).unwrap();
    let n2 = normalize_allele(&n1, &p, &opts).unwrap();
    assert_eq!(n1, n2);
    // Expanding the RLE back to a literal and normalizing again gives the same result.
    let loc = n1.sequence_location().unwrap();
    let (s, e) = loc.exact_interval().unwrap();
    let reference = p.sequence(&acc, s as u64, e as u64).unwrap();
    let expanded =
        expand_reference_length_expression(&reference, n1.state().as_reference_length().unwrap())
            .unwrap();
    let literal_again = Allele::new(loc.clone(), SequenceExpression::from(expanded));
    assert_eq!(normalize_allele(&literal_again, &p, &opts).unwrap(), n1);
}

#[test]
fn iri_locations_and_non_literal_states_pass_through() {
    let p = InMemorySequenceProvider::new();
    let opts = NormalizeOptions::default();
    let a = Allele::new(
        Iri::new("ga4gh:SL.4t6JnYWqHwYw9WzBT_lmWBb3tLQNalkT"),
        SequenceExpression::literal("T").unwrap(),
    );
    assert_eq!(normalize_allele(&a, &p, &opts).unwrap(), a);
    let acc = RefgetAccession::parse(CHR1).unwrap();
    let r = Allele::new(
        SequenceLocation::new(SequenceReference::new(acc), 10, 13).unwrap(),
        ReferenceLengthExpression::new(5, 3).unwrap(),
    );
    assert_eq!(normalize_allele(&r, &p, &opts).unwrap(), r);
}

#[test]
fn ga4gh_sequence_iri_references_are_resolved() {
    let mut p = InMemorySequenceProvider::new();
    let acc = p.insert_sequence(b"ACGTTTTTGA").unwrap();
    let loc = SequenceLocation::new(Iri::new(format!("ga4gh:{acc}")), 4, 5).unwrap();
    let a = Allele::new(loc, SequenceExpression::literal("").unwrap());
    let n = normalize_allele(&a, &p, &NormalizeOptions::default()).unwrap();
    // The deleted T lies in the T5 run [3, 8): the normalized allele spans the run and
    // leaves four of its five residues.
    assert_eq!(span(&n), (3.into(), 8.into()));
    assert_eq!(rle(&n), (4, 1, Some("TTTT".into())));
}

#[test]
fn circular_references_are_rejected_and_unknown_sequences_error() {
    let mut p = InMemorySequenceProvider::new();
    let acc = p.insert_sequence(b"ACGT").unwrap();
    let circular = SequenceReference::new(acc).with_circular(true);
    let a = Allele::new(
        SequenceLocation::new(circular, 1, 2).unwrap(),
        SequenceExpression::literal("T").unwrap(),
    );
    assert!(matches!(
        normalize_allele(&a, &p, &NormalizeOptions::default()),
        Err(kira_vrs::error::NormalizeError::Unsupported(_))
    ));
    let unknown = RefgetAccession::parse(CHR6).unwrap();
    let b = allele(unknown, 1, 2, "T");
    assert!(matches!(
        normalize_allele(&b, &p, &NormalizeOptions::default()),
        Err(kira_vrs::error::NormalizeError::Sequence(_))
    ));
}

#[test]
fn cis_phased_block_members_are_normalized_and_sorted() {
    let mut p = InMemorySequenceProvider::new();
    let acc = p.insert_sequence(b"GGTATATATACCAAAAGG").unwrap();
    let reference = SequenceReference::new(acc);
    let member = |s: i64, e: i64, alt: &str| {
        Allele::new(
            SequenceLocation::from_parts(None, Some(s.into()), Some(e.into())).unwrap(),
            SequenceExpression::literal(alt).unwrap(),
        )
    };
    let block = CisPhasedBlock::new(vec![member(2, 2, "TA").into(), member(13, 14, "").into()])
        .unwrap()
        .with_sequence_reference(reference.clone());
    let n = normalize_cis_phased_block(&block, &p, &NormalizeOptions::default()).unwrap();
    // Members keep the block's implicit reference (no sequenceReference of their own).
    for m in n.members() {
        let a = m.as_object().unwrap();
        assert!(
            a.sequence_location()
                .unwrap()
                .sequence_reference()
                .is_none()
        );
        assert!(a.state().as_reference_length().is_some());
    }
    let digests: Vec<_> = n
        .members()
        .iter()
        .map(|m| m.as_object().unwrap().digest())
        .collect();
    assert!(digests.windows(2).all(|w| w[0] <= w[1]));
    // The block digest is independent of member order.
    let reversed = CisPhasedBlock::new(n.members().iter().rev().cloned().collect()).unwrap();
    assert_eq!(reversed.digest(), n.digest());
}
