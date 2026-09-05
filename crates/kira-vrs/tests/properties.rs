//! Property-based tests (proptest).

use kira_vrs::digest::{DigestSerialize, Identifiable};
use kira_vrs::json::{from_str, to_string};
use kira_vrs::model::*;
use kira_vrs::normalize::{InMemorySequenceProvider, NormalizeOptions, normalize_allele};
use proptest::prelude::*;

fn dna(len: impl Into<proptest::sample::SizeRange>) -> impl Strategy<Value = String> {
    proptest::collection::vec(prop_oneof![Just('A'), Just('C'), Just('G'), Just('T')], len)
        .prop_map(|v| v.into_iter().collect())
}

fn residues() -> impl Strategy<Value = String> {
    "[A-Z*\\-]{0,40}"
}

fn reference() -> SequenceReference {
    SequenceReference::parse("SQ.IIB53T8CNeJJdUqzn9V_JnRtQadwWCbl").unwrap()
}

proptest! {
    #[test]
    fn range_invariants(min in any::<Option<i64>>(), max in any::<Option<i64>>()) {
        if let Ok(r) = Range::new(min, max) {
            prop_assert!(r.min().is_some() || r.max().is_some());
            if let (Some(a), Some(b)) = (r.min(), r.max()) {
                prop_assert!(a <= b);
            }
            prop_assert_eq!(r.min(), min);
            prop_assert_eq!(r.max(), max);
        } else {
            let reserved = [min, max].iter().flatten().any(|v| *v == i64::MIN || *v == i64::MAX);
            let inverted = matches!((min, max), (Some(a), Some(b)) if a > b);
            prop_assert!((min.is_none() && max.is_none()) || inverted || reserved);
        }
    }

    #[test]
    fn sequence_string_round_trip(s in residues()) {
        let seq = SequenceString::new(&s).unwrap();
        prop_assert_eq!(seq.as_str(), s.as_str());
        prop_assert_eq!(seq.len(), s.len());
        let json = to_string(&seq).unwrap();
        let back: SequenceString = from_str(&json).unwrap();
        prop_assert_eq!(back, seq);
    }

    #[test]
    fn sequence_string_rejects_non_residues(s in "[^A-Z*\\-]{1,5}") {
        prop_assert!(SequenceString::new(&s).is_err());
    }

    #[test]
    fn location_start_never_after_end(start in 0i64..1_000_000, end in 0i64..1_000_000) {
        match SequenceLocation::new(reference(), start, end) {
            Ok(loc) => prop_assert!(start <= end && loc.exact_length() == Some((end - start) as u64)),
            Err(_) => prop_assert!(start > end),
        }
    }

    #[test]
    fn allele_json_round_trip_and_deterministic_identifier(
        start in 0i64..250_000_000,
        ref_len in 0i64..60,
        alt in dna(0..60usize),
    ) {
        let loc = SequenceLocation::new(reference(), start, start + ref_len).unwrap();
        let allele = Allele::new(loc, SequenceExpression::from(SequenceString::new(&alt).unwrap()));
        let json = to_string(&allele).unwrap();
        let back: Allele = from_str(&json).unwrap();
        prop_assert_eq!(&back, &allele);
        prop_assert_eq!(back.identifier(), allele.identifier());
        prop_assert_eq!(back.digest_serialization(), allele.digest_serialization());
        // Digest serialization is canonical JSON: no whitespace, fixed key order.
        let ser = String::from_utf8(allele.digest_serialization()).unwrap();
        let canonical = ser.starts_with("{\"location\":\"") && ser.ends_with("\"type\":\"Allele\"}");
        prop_assert!(canonical, "not canonical: {}", ser);
        // Polymorphic parse agrees.
        let v: Variation = from_str(&json).unwrap();
        prop_assert_eq!(v.identifier(), allele.identifier());
    }

    /// Two textual representations of the same insertion (the variant expressed at an
    /// equivalent shifted position) normalize to the same identifier, and normalization is
    /// idempotent.
    #[test]
    fn normalization_is_shift_invariant_and_idempotent(
        seq in dna(20..80usize),
        pos in 0usize..60,
        ins in dna(1..6usize),
    ) {
        let pos = pos.min(seq.len());
        let mut provider = InMemorySequenceProvider::new();
        let acc = provider.insert_sequence(seq.as_bytes()).unwrap();
        let opts = NormalizeOptions::default();
        let mk = |p: usize, alt: &str| {
            Allele::new(
                SequenceLocation::new(SequenceReference::new(acc), p as i64, p as i64).unwrap(),
                SequenceExpression::from(SequenceString::new(alt).unwrap()),
            )
        };
        let a = mk(pos, &ins);
        let n = normalize_allele(&a, &provider, &opts).unwrap();
        prop_assert_eq!(normalize_allele(&n, &provider, &opts).unwrap(), n.clone());

        // Left-shift by one if the base before `pos` equals the last inserted base.
        if pos > 0 && seq.as_bytes()[pos - 1] == *ins.as_bytes().last().unwrap() {
            let shifted = format!("{}{}", &ins[ins.len() - 1..], &ins[..ins.len() - 1]);
            let b = mk(pos - 1, &shifted);
            prop_assert_eq!(normalize_allele(&b, &provider, &opts).unwrap().identifier(), n.identifier());
        }
        // Right-shift by one if the base at `pos` equals the first inserted base.
        if pos < seq.len() && seq.as_bytes()[pos] == ins.as_bytes()[0] {
            let shifted = format!("{}{}", &ins[1..], &ins[..1]);
            let b = mk(pos + 1, &shifted);
            prop_assert_eq!(normalize_allele(&b, &provider, &opts).unwrap().identifier(), n.identifier());
        }
        // VCF-style anchored representation (ref base at pos-1 kept in both alleles).
        if pos > 0 {
            let anchored = Allele::new(
                SequenceLocation::new(SequenceReference::new(acc), pos as i64 - 1, pos as i64).unwrap(),
                SequenceExpression::from(SequenceString::new(&format!("{}{ins}", &seq[pos - 1..pos])).unwrap()),
            );
            prop_assert_eq!(normalize_allele(&anchored, &provider, &opts).unwrap().identifier(), n.identifier());
        }
    }

    /// Deletions: the same deleted sequence expressed at any equivalent position normalizes
    /// identically, and the reference-length expression expands back to the alternate.
    #[test]
    fn deletion_normalization_is_consistent(seq in dna(20..80usize), start in 0usize..60, len in 1usize..6) {
        let start = start.min(seq.len().saturating_sub(1));
        let end = (start + len).min(seq.len());
        let mut provider = InMemorySequenceProvider::new();
        let acc = provider.insert_sequence(seq.as_bytes()).unwrap();
        let opts = NormalizeOptions::default();
        let del = |s: usize, e: usize| {
            Allele::new(
                SequenceLocation::new(SequenceReference::new(acc), s as i64, e as i64).unwrap(),
                SequenceExpression::literal("").unwrap(),
            )
        };
        let n = normalize_allele(&del(start, end), &provider, &opts).unwrap();
        let rle = n.state().as_reference_length().expect("deletions normalize to RLE");
        let (ns, ne) = n.sequence_location().unwrap().exact_interval().unwrap();
        let region = &seq.as_bytes()[ns as usize..ne as usize];
        let expanded = kira_vrs::normalize::expand_reference_length_expression(region, rle).unwrap();
        // Applying the normalized deletion yields the same molecule as the input deletion.
        let mut applied_input = seq.as_bytes()[..start].to_vec();
        applied_input.extend_from_slice(&seq.as_bytes()[end..]);
        let mut applied_norm = seq.as_bytes()[..ns as usize].to_vec();
        applied_norm.extend_from_slice(expanded.as_bytes());
        applied_norm.extend_from_slice(&seq.as_bytes()[ne as usize..]);
        prop_assert_eq!(applied_input, applied_norm);
        // Right-shift by one when possible: same identifier.
        if end < seq.len() && seq.as_bytes()[end] == seq.as_bytes()[start] {
            let m = normalize_allele(&del(start + 1, end + 1), &provider, &opts).unwrap();
            prop_assert_eq!(m.identifier(), n.identifier());
        }
    }
}
