//! Parse the specification's example objects (`examples/*.json`) and reproduce the
//! identifiers printed in the specification text.

use kira_vrs::digest::Identifiable;
use kira_vrs::model::*;
use kira_vrs_validation::{EXAMPLES, example_text};

#[test]
fn examples_parse_into_the_declared_class() {
    for ex in EXAMPLES {
        let text = example_text(ex.file);
        let result: Result<(), String> = match ex.class {
            "Adjacency" => kira_vrs::json::from_str::<Adjacency>(&text)
                .map(|_| ())
                .map_err(|e| e.to_string()),
            "Terminus" => kira_vrs::json::from_str::<Terminus>(&text)
                .map(|_| ())
                .map_err(|e| e.to_string()),
            "DerivativeMolecule" => kira_vrs::json::from_str::<DerivativeMolecule>(&text)
                .map(|_| ())
                .map_err(|e| e.to_string()),
            "CisPhasedBlock" => kira_vrs::json::from_str::<CisPhasedBlock>(&text)
                .map(|_| ())
                .map_err(|e| e.to_string()),
            "Allele" => kira_vrs::json::from_str::<Allele>(&text)
                .map(|_| ())
                .map_err(|e| e.to_string()),
            other => panic!("unexpected class {other}"),
        };
        if ex.should_fail {
            assert!(result.is_err(), "{} should be rejected", ex.file);
        } else {
            assert!(
                result.is_ok(),
                "{} failed: {}",
                ex.file,
                result.unwrap_err()
            );
            // Also parse polymorphically and check the class.
            let v: Variation = kira_vrs::json::from_str(&text).unwrap();
            assert_eq!(v.type_name(), ex.class);
        }
    }
}

/// Identifiers quoted in the VRS 2.1.0 documentation for these examples.
#[test]
fn examples_reproduce_documented_identifiers() {
    let expansion: Allele = kira_vrs::json::from_str(&example_text("SPDI_expansion.json")).unwrap();
    assert_eq!(
        expansion.identifier().to_string(),
        "ga4gh:VA.Oop4kjdTtKcg1kiZjIJAAR3bp7qi4aNT"
    );

    let haplotype: CisPhasedBlock =
        kira_vrs::json::from_str(&example_text("simple_haplotype.json")).unwrap();
    assert_eq!(
        haplotype.identifier().to_string(),
        "ga4gh:CPB.YAWwnFF0e-T7fnuT4wRzZW4Lzg7jc-zQ"
    );

    let dm: DerivativeMolecule =
        kira_vrs::json::from_str(&example_text("sv_derivative_molecule.json")).unwrap();
    assert_eq!(dm.identifier().prefix(), TypePrefix::DerivativeMolecule);
    assert_eq!(dm.components().len(), 5);

    // Adjacency example from the Adjacency concept page (order: chr2 start, chr1 end).
    let adjacency = Adjacency::new(
        SequenceLocation::starting_at(
            SequenceReference::parse("SQ.9KdcA9ZpY1Cpvxvg8bMSLYDUpsX6GDLO").unwrap(),
            456,
        )
        .unwrap(),
        SequenceLocation::ending_at(
            SequenceReference::parse("SQ.S_KjnFVz-FE7M0W6yoaUDgYxLPc1jyWU").unwrap(),
            123,
        )
        .unwrap(),
    )
    .unwrap();
    assert_eq!(
        adjacency.identifier().to_string(),
        "ga4gh:AJ.O0IbSYyhnBAtUsR51bpdoqeSo4YaDMFo"
    );
}

/// The invalid example must be rejected for the right reason.
#[test]
fn invalid_adjacency_reports_start_and_end() {
    let err =
        kira_vrs::json::from_str::<Adjacency>(&example_text("invalid_adjacency.json")).unwrap_err();
    assert!(err.to_string().contains("both start and end"), "{err}");
}
