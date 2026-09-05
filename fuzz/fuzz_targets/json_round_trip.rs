//! Structure-aware fuzzing: build alleles from arbitrary parts, serialize, parse, and check
//! that identifiers and digest serializations are stable across the round trip.
#![no_main]

use arbitrary::Arbitrary;
use kira_vrs::digest::{DigestSerialize, Identifiable};
use kira_vrs::model::*;
use libfuzzer_sys::fuzz_target;

#[derive(Arbitrary, Debug)]
struct Input {
    start: u32,
    len: u8,
    alt: Vec<u8>,
    start_range: Option<(bool, bool)>,
    id: Option<String>,
    rle: Option<(u16, u8)>,
}

fuzz_target!(|input: Input| {
    let reference = SequenceReference::parse("SQ.IIB53T8CNeJJdUqzn9V_JnRtQadwWCbl").unwrap();
    let start = i64::from(input.start);
    let end = start + i64::from(input.len);
    let start_pos: IntOrRange = match input.start_range {
        Some((true, _)) => Range::at_most(start).unwrap().into(),
        Some((false, true)) => Range::at_least(start).unwrap().into(),
        _ => start.into(),
    };
    let location = SequenceLocation::new(reference, start_pos, end).unwrap();
    let state: SequenceExpression = match input.rle {
        Some((len, unit)) => ReferenceLengthExpression::new(i64::from(len), i64::from(unit)).unwrap().into(),
        None => match SequenceString::from_bytes(&input.alt) {
            Ok(s) => s.into(),
            Err(_) => return,
        },
    };
    let mut allele = Allele::new(location, state);
    if let Some(id) = input.id {
        allele = allele.with_id(id);
    }
    let json = kira_vrs::json::to_string(&allele).unwrap();
    let back: Allele = kira_vrs::json::from_str(&json).unwrap();
    assert_eq!(back, allele);
    assert_eq!(back.identifier(), allele.identifier());
    assert_eq!(back.digest_serialization(), allele.digest_serialization());
    let v: Variation = kira_vrs::json::from_str(&json).unwrap();
    assert_eq!(v.identifier(), allele.identifier());
});
