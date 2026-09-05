//! Fuzz digest serialization: for any parseable variation, the digest serialization must be
//! valid canonical JSON (parseable, no whitespace, sorted keys) and identifiers must be
//! deterministic.
#![no_main]

use kira_vrs::digest::{DigestSerialize, Identifiable};
use kira_vrs::model::Variation;
use libfuzzer_sys::fuzz_target;

fuzz_target!(|data: &[u8]| {
    let Ok(v) = kira_vrs::json::from_slice::<Variation>(data) else { return };
    let blob = v.digest_serialization();
    let text = std::str::from_utf8(&blob).expect("UTF-8");
    let value: serde_json::Value = serde_json::from_str(text).expect("canonical JSON parses");
    // Canonical JSON re-encoded compactly by serde_json (which sorts nothing) must round-trip
    // byte-for-byte when keys are already sorted and there is no whitespace.
    assert!(!text.contains('\n') && !text.contains(": "));
    if let serde_json::Value::Object(map) = &value {
        let keys: Vec<&String> = map.keys().collect();
        let mut sorted = keys.clone();
        sorted.sort();
        assert_eq!(keys, sorted);
        assert_eq!(map["type"].as_str().unwrap(), v.type_name());
    }
    assert_eq!(v.identifier(), v.clone().identifier());
    assert_eq!(v.digest().as_str().len(), 32);
});
