//! Run the official GA4GH VRS validation vectors (`validation/models.yaml`,
//! `validation/functions.yaml`) against kira-vrs.
//!
//! Every vector is executed: for each class the input JSON is deserialized into the typed
//! model, and every expected output (`ga4gh_serialize`, `ga4gh_digest`, `ga4gh_identify`, and
//! the VRS 1.3 variants) is compared byte-for-byte.

use kira_vrs::digest::{DigestSerialize, Identifiable, legacy, sha512t24u};
use kira_vrs::model::*;
use kira_vrs_validation::{model_vectors, sha512t24u_vectors};
use serde_json::Value;

#[test]
fn sha512t24u_function_vectors() {
    let vectors = sha512t24u_vectors();
    assert!(!vectors.is_empty());
    for v in vectors {
        assert_eq!(
            sha512t24u(v.input.blob.as_bytes()).as_str(),
            v.out,
            "blob {:?}",
            v.input.blob
        );
    }
}

/// Everything the harness can compute for one deserialized object.
struct Computed {
    serialize: Vec<u8>,
    digest: Option<String>,
    identify: Option<String>,
    v1_3_serialize: Option<Vec<u8>>,
    v1_3_identify: Option<String>,
}

fn identifiable<T: Identifiable>(obj: &T) -> Computed {
    Computed {
        serialize: obj.digest_serialization(),
        digest: Some(obj.digest().to_string()),
        identify: Some(obj.identifier().to_string()),
        v1_3_serialize: None,
        v1_3_identify: None,
    }
}

fn value_object<T: DigestSerialize>(obj: &T) -> Computed {
    Computed {
        serialize: obj.digest_serialization(),
        digest: None,
        identify: None,
        v1_3_serialize: None,
        v1_3_identify: None,
    }
}

fn compute(class: &str, input: &Value) -> Computed {
    let text = input.to_string();
    match class {
        "SequenceReference" => value_object(&parse::<SequenceReference>(&text)),
        "LengthExpression" => value_object(&parse::<LengthExpression>(&text)),
        "LiteralSequenceExpression" => value_object(&parse::<LiteralSequenceExpression>(&text)),
        "ReferenceLengthExpression" => value_object(&parse::<ReferenceLengthExpression>(&text)),
        "SequenceOffsetLocation" => value_object(&parse::<SequenceOffsetLocation>(&text)),
        "SequenceLocation" => {
            let loc = parse::<SequenceLocation>(&text);
            let mut c = identifiable(&loc);
            c.v1_3_serialize = legacy::sequence_location_1_3_serialization(&loc).ok();
            c.v1_3_identify = legacy::sequence_location_1_3_identifier(&loc)
                .ok()
                .map(|i| i.to_string());
            c
        }
        "RelativeSequenceLocation" => identifiable(&parse::<RelativeSequenceLocation>(&text)),
        "Allele" => {
            let allele = parse::<Allele>(&text);
            let mut c = identifiable(&allele);
            c.v1_3_serialize = legacy::allele_1_3_serialization(&allele).ok();
            c.v1_3_identify = legacy::allele_1_3_identifier(&allele)
                .ok()
                .map(|i| i.to_string());
            c
        }
        "RelativeAllele" => identifiable(&parse::<RelativeAllele>(&text)),
        "CisPhasedBlock" => identifiable(&parse::<CisPhasedBlock>(&text)),
        "Adjacency" => identifiable(&parse::<Adjacency>(&text)),
        "Terminus" => identifiable(&parse::<Terminus>(&text)),
        "DerivativeMolecule" => identifiable(&parse::<DerivativeMolecule>(&text)),
        "CopyNumberCount" => identifiable(&parse::<CopyNumberCount>(&text)),
        "CopyNumberChange" => identifiable(&parse::<CopyNumberChange>(&text)),
        other => panic!("validation vector for unknown class {other}"),
    }
}

fn parse<T: serde::de::DeserializeOwned>(text: &str) -> T {
    kira_vrs::json::from_str(text).unwrap_or_else(|e| panic!("deserialize {text}: {e}"))
}

#[test]
fn model_vectors_all_pass() {
    let mut executed = 0usize;
    let mut classes = Vec::new();
    for (class, vectors) in model_vectors() {
        classes.push(class.clone());
        for (i, vector) in vectors.iter().enumerate() {
            let label = vector
                .name
                .clone()
                .unwrap_or_else(|| format!("{class}[{i}]"));
            let computed = compute(&class, &vector.input);
            for (function, expected) in &vector.out {
                let actual: Option<String> = match function.as_str() {
                    "ga4gh_serialize" => {
                        Some(String::from_utf8(computed.serialize.clone()).unwrap())
                    }
                    "ga4gh_digest" => computed.digest.clone(),
                    "ga4gh_identify" => computed.identify.clone(),
                    "ga4gh_1_3_serialize" => computed
                        .v1_3_serialize
                        .clone()
                        .map(|b| String::from_utf8(b).unwrap()),
                    "ga4gh_1_3_identify" => computed.v1_3_identify.clone(),
                    "ga4gh_1_3_digest" => computed
                        .v1_3_identify
                        .clone()
                        .map(|id| id.rsplit('.').next().unwrap().to_owned()),
                    other => panic!("{label}: unknown validation function {other}"),
                };
                let expected: Option<String> = match expected {
                    Value::Null => None,
                    Value::String(s) => Some(s.clone()),
                    other => panic!("{label}: unexpected expectation {other}"),
                };
                assert_eq!(actual, expected, "{label}: {function}");
                executed += 1;
            }
        }
    }
    // Guard against silently running an empty or truncated suite.
    assert!(executed >= 60, "only {executed} expectations executed");
    for required in [
        "Allele",
        "SequenceLocation",
        "CisPhasedBlock",
        "Adjacency",
        "RelativeAllele",
    ] {
        assert!(
            classes.iter().any(|c| c == required),
            "missing class {required}"
        );
    }
}

/// The same vectors, parsed polymorphically as `Variation` / `Location` where applicable,
/// must produce the same identifiers.
#[test]
fn model_vectors_via_polymorphic_unions() {
    for (class, vectors) in model_vectors() {
        for vector in &vectors {
            let text = vector.input.to_string();
            let expected = vector.out.get("ga4gh_identify").and_then(Value::as_str);
            let Some(expected) = expected else { continue };
            match class.as_str() {
                "SequenceLocation" | "RelativeSequenceLocation" => {
                    let loc: Location = parse(&text);
                    assert_eq!(loc.identifier().to_string(), expected);
                    assert_eq!(loc.type_name(), class);
                }
                _ => {
                    let v: Variation = parse(&text);
                    assert_eq!(v.identifier().to_string(), expected);
                    assert_eq!(v.type_name(), class);
                }
            }
        }
    }
}

/// JSON round trip of every vector input: serialize the parsed object back to JSON, parse it
/// again, and require an equal object with an equal digest serialization.
#[test]
fn model_vectors_round_trip_json() {
    for (class, vectors) in model_vectors() {
        for vector in &vectors {
            let text = vector.input.to_string();
            macro_rules! round_trip {
                ($ty:ty) => {{
                    let a: $ty = parse(&text);
                    let json = kira_vrs::json::to_string(&a).unwrap();
                    let b: $ty = parse(&json);
                    assert_eq!(a, b, "{class} round trip");
                    assert_eq!(a.digest_serialization(), b.digest_serialization());
                }};
            }
            match class.as_str() {
                "SequenceReference" => round_trip!(SequenceReference),
                "LengthExpression" => round_trip!(LengthExpression),
                "LiteralSequenceExpression" => round_trip!(LiteralSequenceExpression),
                "ReferenceLengthExpression" => round_trip!(ReferenceLengthExpression),
                "SequenceOffsetLocation" => round_trip!(SequenceOffsetLocation),
                "SequenceLocation" | "RelativeSequenceLocation" => round_trip!(Location),
                _ => round_trip!(Variation),
            }
        }
    }
}
