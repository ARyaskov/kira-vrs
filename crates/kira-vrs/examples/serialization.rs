//! JSON interchange: serialize a typed allele to schema-conformant VRS JSON, parse VRS JSON
//! back into typed objects (including polymorphic `Variation`), and see how invalid input is
//! rejected.
//!
//! ```text
//! cargo run --example serialization
//! ```

use kira_vrs::json;
use kira_vrs::prelude::*;

fn main() -> Result<(), kira_vrs::Error> {
    let chr19 = SequenceReference::parse("SQ.IIB53T8CNeJJdUqzn9V_JnRtQadwWCbl")?;
    let allele = Allele::new(
        SequenceLocation::new(chr19, 44_908_821, 44_908_822)?,
        SequenceExpression::literal("T")?,
    )
    .with_id(Allele::identifier_placeholder());

    // Typed object -> JSON (compact or pretty).
    println!("{}", json::to_string_pretty(&allele)?);

    // JSON -> typed object, validated on the way in.
    let text = json::to_string(&allele)?;
    let parsed: Allele = json::from_str(&text)?;
    assert_eq!(parsed, allele);

    // Polymorphic parsing: `Variation` dispatches on `type`.
    let variation: Variation = json::from_str(&text)?;
    println!(
        "parsed a {} with identifier {}",
        variation.type_name(),
        variation.identifier()
    );

    // The SPDI-expansion example from the specification (a ReferenceLengthExpression state).
    let spdi = r#"{
      "type": "Allele",
      "expressions": [{"syntax": "spdi", "value": "NC_000001.11:40819438:CTCCTCCT:CTCCTCCTCCT"}],
      "location": {
        "type": "SequenceLocation",
        "sequenceReference": {"type": "SequenceReference", "refgetAccession": "SQ.Ya6Rs7DHhDeg7YaOSg1EoNi3U_nQ9SvO", "id": "NC_000001.11"},
        "start": 40819438, "end": 40819446
      },
      "state": {"type": "ReferenceLengthExpression", "length": 11, "repeatSubunitLength": 3}
    }"#;
    let expansion: Allele = json::from_str(spdi)?;
    println!(
        "SPDI expansion: {} ({})",
        expansion.identifier(),
        expansion.expressions()[0].value()
    );

    // Invalid input is rejected with a specific error, never a panic.
    for bad in [
        r#"{"type":"Allele","location":{"type":"SequenceLocation","start":10,"end":5},"state":{"type":"LiteralSequenceExpression","sequence":"T"}}"#,
        r#"{"type":"Allele","location":"ga4gh:SL.wIlaGykfwHIpPY2Fcxtbx4TINbbODFVz","state":{"type":"LiteralSequenceExpression","sequence":"acgt"}}"#,
        r#"{"type":"CopyNumberChange","location":"ga4gh:SL.wIlaGykfwHIpPY2Fcxtbx4TINbbODFVz","copyChange":"deleted"}"#,
    ] {
        match json::from_str::<Variation>(bad) {
            Ok(_) => unreachable!(),
            Err(e) => println!("rejected ({:?}): {e}", e.kind()),
        }
    }
    Ok(())
}

trait Placeholder {
    fn identifier_placeholder() -> &'static str;
}
impl Placeholder for Allele {
    fn identifier_placeholder() -> &'static str {
        "ga4gh:VA.0AePZIWZUNsUlQTamyLrjm2HWUw2opLt"
    }
}
