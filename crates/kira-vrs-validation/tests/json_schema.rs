//! Validate kira-vrs JSON output against the official VRS 2.1.0 JSON Schemas.
//!
//! The upstream schemas reference each other with root-relative `$ref`s
//! (`/ga4gh/schema/vrs/2.1.0/json/SequenceLocation`); a retriever maps those to the vendored
//! files.

use std::sync::LazyLock;

use jsonschema::{Retrieve, Uri, Validator};
use kira_vrs::model::*;
use kira_vrs_validation::{EXAMPLES, example_text, model_vectors, schema, schemas};
use serde_json::Value;

struct VendoredSchemas;

impl Retrieve for VendoredSchemas {
    fn retrieve(
        &self,
        uri: &Uri<String>,
    ) -> Result<Value, Box<dyn std::error::Error + Send + Sync>> {
        // e.g. https://w3id.org/ga4gh/schema/vrs/2.1.0/json/SequenceLocation
        let path = uri.path().as_str();
        let mut parts = path.trim_start_matches('/').split('/');
        let (
            Some("ga4gh"),
            Some("schema"),
            Some(module),
            Some(_version),
            Some("json"),
            Some(class),
        ) = (
            parts.next(),
            parts.next(),
            parts.next(),
            parts.next(),
            parts.next(),
            parts.next(),
        )
        else {
            return Err(format!("unexpected schema reference {uri}").into());
        };
        Ok(schema(module, class))
    }
}

fn validator_for(class: &str) -> Validator {
    jsonschema::options()
        .with_retriever(VendoredSchemas)
        .build(&schema("vrs", class))
        .unwrap_or_else(|e| panic!("build validator for {class}: {e}"))
}

static VALIDATORS: LazyLock<Vec<(String, Validator)>> = LazyLock::new(|| {
    schemas("vrs")
        .into_iter()
        .filter(|(name, s)| s.get("properties").is_some() && name != "Expression")
        .map(|(name, _)| {
            let v = validator_for(&name);
            (name, v)
        })
        .collect()
});

fn validator(class: &str) -> &'static Validator {
    &VALIDATORS
        .iter()
        .find(|(n, _)| n == class)
        .unwrap_or_else(|| panic!("no schema {class}"))
        .1
}

fn assert_valid(class: &str, json: &str) {
    let value: Value = serde_json::from_str(json).unwrap();
    let errors: Vec<String> = validator(class)
        .iter_errors(&value)
        .map(|e| e.to_string())
        .collect();
    assert!(
        errors.is_empty(),
        "{class} output violates the schema:\n{json}\n{errors:#?}"
    );
}

/// Every validation-vector input, re-serialized by kira-vrs, must validate against the
/// official schema of its class.
#[test]
fn reserialized_vectors_validate_against_upstream_schema() {
    for (class, vectors) in model_vectors() {
        for vector in &vectors {
            let text = vector.input.to_string();
            let json = match class.as_str() {
                "SequenceReference" => reserialize::<SequenceReference>(&text),
                "LengthExpression" => reserialize::<LengthExpression>(&text),
                "LiteralSequenceExpression" => reserialize::<LiteralSequenceExpression>(&text),
                "ReferenceLengthExpression" => reserialize::<ReferenceLengthExpression>(&text),
                "SequenceOffsetLocation" => reserialize::<SequenceOffsetLocation>(&text),
                "SequenceLocation" | "RelativeSequenceLocation" => reserialize::<Location>(&text),
                _ => reserialize::<Variation>(&text),
            };
            assert_valid(&class, &json);
        }
    }
}

fn reserialize<T: serde::de::DeserializeOwned + serde::Serialize>(text: &str) -> String {
    let obj: T = kira_vrs::json::from_str(text).unwrap();
    kira_vrs::json::to_string(&obj).unwrap()
}

/// Upstream examples validate (or fail) against the schema exactly as upstream's own test
/// suite expects, and kira-vrs's re-serialization of the valid ones validates too.
#[test]
fn examples_match_upstream_schema_expectations() {
    for ex in EXAMPLES {
        let text = example_text(ex.file);
        let value: Value = serde_json::from_str(&text).unwrap();
        let valid = validator(ex.class).is_valid(&value);
        assert_eq!(valid, !ex.should_fail, "{}", ex.file);
        if !ex.should_fail {
            let json = reserialize::<Variation>(&text);
            assert_valid(ex.class, &json);
        }
    }
}

/// Objects built through the Rust API, with decorative metadata, produce schema-valid JSON.
#[test]
fn constructed_objects_with_metadata_validate() {
    let reference = SequenceReference::parse("SQ.IIB53T8CNeJJdUqzn9V_JnRtQadwWCbl")
        .unwrap()
        .with_id("NC_000019.10")
        .with_residue_alphabet(ResidueAlphabet::NucleicAcid)
        .with_molecule_type(MoleculeType::Genomic)
        .with_circular(false);
    let location = SequenceLocation::new(reference, 44_908_821, 44_908_822)
        .unwrap()
        .with_sequence(SequenceString::new("C").unwrap())
        .with_name("rs7412 locus");
    let allele = Allele::new(location, SequenceExpression::literal("T").unwrap())
        .with_id("ga4gh:VA.0AePZIWZUNsUlQTamyLrjm2HWUw2opLt")
        .with_name("rs7412")
        .with_description("APOE e2 allele")
        .with_alias("rs7412-T")
        .with_extension(Extension::new("clinvar_id", 17864).with_description("ClinVar VCV"))
        .with_expression(
            Expression::new(Syntax::HgvsG, "NC_000019.10:g.44908822C>T")
                .with_syntax_version("21.0"),
        );
    assert_valid("Allele", &kira_vrs::json::to_string(&allele).unwrap());

    let cnc = CopyNumberCount::new(
        Iri::new("ga4gh:SL.4t6JnYWqHwYw9WzBT_lmWBb3tLQNalkT"),
        Range::at_least(3).unwrap(),
    )
    .unwrap();
    assert_valid("CopyNumberCount", &kira_vrs::json::to_string(&cnc).unwrap());

    let rle = ReferenceLengthExpression::new(11, 3)
        .unwrap()
        .with_sequence(SequenceString::new("CTCCTCCTCCT").unwrap());
    assert_valid(
        "ReferenceLengthExpression",
        &kira_vrs::json::to_string(&rle).unwrap(),
    );
}

/// The inherent-property tables hard-coded in the digest serializers must agree with the
/// `ga4gh.inherent` annotations of the upstream schemas (sorted, since RFC 8785 orders keys).
#[test]
fn inherent_property_tables_match_upstream_schema() {
    let expected: &[(&str, &[&str])] = &[
        ("Allele", &["location", "state", "type"]),
        (
            "RelativeAllele",
            &["baseState", "mappedState", "relativeLocation", "type"],
        ),
        ("CisPhasedBlock", &["members", "type"]),
        ("Adjacency", &["adjoinedSequences", "linker", "type"]),
        ("Terminus", &["location", "type"]),
        ("TraversalBlock", &["component", "orientation", "type"]),
        ("DerivativeMolecule", &["components", "type"]),
        ("CopyNumberCount", &["copies", "location", "type"]),
        ("CopyNumberChange", &["copyChange", "location", "type"]),
        (
            "SequenceLocation",
            &["end", "sequenceReference", "start", "type"],
        ),
        (
            "RelativeSequenceLocation",
            &["baseSequenceLocation", "mappedSequenceLocation", "type"],
        ),
        (
            "SequenceOffsetLocation",
            &[
                "anchor",
                "anchorOrientation",
                "offsetEnd",
                "offsetStart",
                "sequenceReference",
                "type",
            ],
        ),
        ("SequenceReference", &["refgetAccession", "type"]),
        ("LiteralSequenceExpression", &["sequence", "type"]),
        (
            "ReferenceLengthExpression",
            &["length", "repeatSubunitLength", "type"],
        ),
        ("LengthExpression", &["length", "type"]),
    ];
    for (class, keys) in expected {
        let s = schema("vrs", class);
        let mut inherent: Vec<String> = s["ga4gh"]["inherent"]
            .as_array()
            .unwrap_or_else(|| panic!("{class} has no ga4gh.inherent"))
            .iter()
            .map(|v| v.as_str().unwrap().to_owned())
            .collect();
        inherent.sort();
        assert_eq!(inherent, keys.to_vec(), "{class}");
    }
}

/// Type prefixes hard-coded in the crate must agree with the schemas.
#[test]
fn type_prefixes_match_upstream_schema() {
    for (class, prefix) in [
        ("Allele", TypePrefix::Allele),
        ("RelativeAllele", TypePrefix::RelativeAllele),
        ("CisPhasedBlock", TypePrefix::CisPhasedBlock),
        ("Adjacency", TypePrefix::Adjacency),
        ("Terminus", TypePrefix::Terminus),
        ("DerivativeMolecule", TypePrefix::DerivativeMolecule),
        ("CopyNumberCount", TypePrefix::CopyNumberCount),
        ("CopyNumberChange", TypePrefix::CopyNumberChange),
        ("SequenceLocation", TypePrefix::SequenceLocation),
        (
            "RelativeSequenceLocation",
            TypePrefix::RelativeSequenceLocation,
        ),
    ] {
        let s = schema("vrs", class);
        assert_eq!(
            s["ga4gh"]["prefix"].as_str().unwrap(),
            prefix.as_str(),
            "{class}"
        );
    }
}

/// Maturity annotations recorded in `kira_vrs::spec` must agree with the schemas.
#[test]
fn maturity_table_matches_upstream_schema() {
    for (class, maturity) in kira_vrs::spec::CLASS_MATURITY {
        let s = schema("vrs", class);
        assert_eq!(
            s["maturity"].as_str().unwrap(),
            maturity.as_str(),
            "{class}"
        );
    }
    let s = schema("vrs", "Allele");
    assert_eq!(
        s["$id"].as_str().unwrap(),
        format!("{}Allele", kira_vrs::spec::VRS_SCHEMA_BASE)
    );
}
