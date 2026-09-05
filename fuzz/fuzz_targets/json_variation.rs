//! Fuzz JSON deserialization of every VRS class: arbitrary bytes must never panic, and
//! whatever parses must re-serialize and re-parse to an equal object.
#![no_main]

use kira_vrs::model::{Location, SequenceExpression, Variation};
use libfuzzer_sys::fuzz_target;

fuzz_target!(|data: &[u8]| {
    if let Ok(v) = kira_vrs::json::from_slice::<Variation>(data) {
        let json = kira_vrs::json::to_string(&v).unwrap();
        let back: Variation = kira_vrs::json::from_str(&json).unwrap();
        assert_eq!(back, v);
        let _ = kira_vrs::digest::Identifiable::identifier(&v);
    }
    if let Ok(l) = kira_vrs::json::from_slice::<Location>(data) {
        let json = kira_vrs::json::to_string(&l).unwrap();
        let back: Location = kira_vrs::json::from_str(&json).unwrap();
        assert_eq!(back, l);
    }
    let _ = kira_vrs::json::from_slice::<SequenceExpression>(data);
});
