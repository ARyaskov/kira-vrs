//! Structural validation: invariants are enforced by constructors, so a value of a model type
//! is always schema-valid; JSON input goes through the same constructors. This example
//! exercises the main rules and prints the resulting errors.
//!
//! The official GA4GH validation vectors are run by the `kira-vrs-validation` crate
//! (`cargo test -p kira-vrs-validation`).
//!
//! ```text
//! cargo run --example validation
//! ```

use kira_vrs::error::{CoordinateError, ModelError};
use kira_vrs::prelude::*;

fn main() {
    let chr1 = SequenceReference::parse("SQ.Ya6Rs7DHhDeg7YaOSg1EoNi3U_nQ9SvO").unwrap();

    // Coordinates.
    let err = SequenceLocation::new(chr1.clone(), 10, 5).unwrap_err();
    assert!(matches!(
        err,
        CoordinateError::StartAfterEnd { start: 10, end: 5 }
    ));
    println!("start > end:            {err}");

    let err = SequenceLocation::from_parts(Some(chr1.clone().into()), None, None).unwrap_err();
    println!("no coordinates:         {err}");

    let err = Range::new(None, None).unwrap_err();
    println!("unbounded range:        {err}");

    let err = SequenceLocation::new(chr1.clone(), -3, 5).unwrap_err();
    println!("negative coordinate:    {err}");

    // Circular sequences may wrap around the origin.
    let mito = SequenceReference::parse("SQ.k3grVkjY-hoWcCUojHw6VU6GE3MZ8Sct")
        .unwrap()
        .with_circular(true);
    let wrap = SequenceLocation::new(mito, 16_560, 10).unwrap();
    println!("circular wrap-around:   {:?}", wrap.exact_interval());

    // Sequences and identifiers.
    println!(
        "lowercase residues:     {}",
        SequenceString::new("acgt").unwrap_err()
    );
    println!(
        "bad refget accession:   {}",
        RefgetAccession::parse("NC_000001.11").unwrap_err()
    );
    println!(
        "bad identifier:         {}",
        VrsIdentifier::parse("ga4gh:VA.short").unwrap_err()
    );

    // Model rules.
    let err = CisPhasedBlock::new(vec![]).unwrap_err();
    assert!(matches!(err, ModelError::TooFewItems { min: 2, .. }));
    println!("cis-phased block:       {err}");

    let both = SequenceLocation::new(chr1.clone(), 100, 200).unwrap();
    let one = SequenceLocation::starting_at(chr1.clone(), 456).unwrap();
    let err = Adjacency::new(both, one.clone()).unwrap_err();
    println!("adjacency:              {err}");
    assert!(Adjacency::new(SequenceLocation::ending_at(chr1, 123).unwrap(), one).is_ok());

    println!(
        "copy change value set:  {}",
        CopyChange::parse("deleted").unwrap_err()
    );
    println!(
        "negative copies:        {}",
        CopyNumberCount::new(Iri::new("x"), -1).unwrap_err()
    );

    // Everything above is also enforced for JSON input.
    let err = kira_vrs::json::from_str::<Allele>(
        r#"{"type":"Allele","location":{"type":"SequenceLocation","start":10,"end":5},"state":{"type":"LiteralSequenceExpression","sequence":"T"}}"#,
    )
    .unwrap_err();
    println!("from JSON:              {err}");
}
