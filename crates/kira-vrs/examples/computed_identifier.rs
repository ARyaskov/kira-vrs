//! Computed identifiers step by step: digest serialization (RFC 8785 canonical JSON of the
//! inherent properties), `sha512t24u`, and the `ga4gh:` CURIE — reproducing the worked
//! example of the VRS *Computed Identifiers* convention.
//!
//! ```text
//! cargo run --example computed_identifier
//! ```

use kira_vrs::digest::{legacy, sha512t24u};
use kira_vrs::prelude::*;

fn main() -> Result<(), kira_vrs::Error> {
    let chr19 = SequenceReference::parse("SQ.IIB53T8CNeJJdUqzn9V_JnRtQadwWCbl")?;
    let location = SequenceLocation::new(chr19, 44_908_821, 44_908_822)?;
    let allele = Allele::new(location, SequenceExpression::literal("T")?);

    // 1. Digest serialization of the nested identifiable object (the location).
    let location = allele.sequence_location().unwrap();
    let loc_blob = location.digest_serialization();
    println!(
        "location serialization: {}",
        String::from_utf8_lossy(&loc_blob)
    );
    println!("location digest:        {}", sha512t24u(&loc_blob));

    // 2. The allele's serialization references the location *by digest*.
    let blob = allele.digest_serialization();
    println!("allele serialization:   {}", String::from_utf8_lossy(&blob));
    assert_eq!(
        String::from_utf8_lossy(&blob),
        r#"{"location":"wIlaGykfwHIpPY2Fcxtbx4TINbbODFVz","state":{"sequence":"T","type":"LiteralSequenceExpression"},"type":"Allele"}"#
    );

    // 3. sha512t24u = SHA-512, truncated to 24 bytes, base64url.
    let digest = sha512t24u(&blob);
    println!("allele digest:          {digest}");
    assert_eq!(digest, allele.digest());

    // 4. Identifier = "ga4gh:" + type prefix + "." + digest.
    let id = allele.identifier();
    println!("allele identifier:      {id}");
    assert_eq!(id.prefix(), TypePrefix::Allele);
    assert_eq!(id.to_string(), "ga4gh:VA.0AePZIWZUNsUlQTamyLrjm2HWUw2opLt");

    // Decorative metadata never changes the identifier.
    let decorated = allele.clone().with_name("rs7412").with_alias("APOE e2");
    assert_eq!(decorated.identifier(), id);

    // Legacy VRS 1.3 identifiers (implementation extension) for migrating older databases.
    println!(
        "VRS 1.3 identifier:     {}",
        legacy::allele_1_3_identifier(&allele)?
    );

    // Identifiers can be parsed back and used as IRIs in place of nested objects.
    let parsed = VrsIdentifier::parse("ga4gh:SL.wIlaGykfwHIpPY2Fcxtbx4TINbbODFVz")?;
    let by_reference = Allele::new(parsed.to_iri(), SequenceExpression::literal("T")?);
    assert_eq!(by_reference.identifier(), id);
    Ok(())
}
