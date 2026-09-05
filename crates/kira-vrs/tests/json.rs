//! JSON interchange: round trips, required/optional properties, and rejection of invalid
//! documents with meaningful errors.

use kira_vrs::digest::Identifiable;
use kira_vrs::json::{JsonErrorKind, from_str, to_string, to_string_pretty};
use kira_vrs::model::*;

const RS7412: &str = r#"{
  "type": "Allele",
  "location": {
    "type": "SequenceLocation",
    "sequenceReference": {
      "type": "SequenceReference",
      "refgetAccession": "SQ.IIB53T8CNeJJdUqzn9V_JnRtQadwWCbl"
    },
    "start": 44908821,
    "end": 44908822
  },
  "state": { "type": "LiteralSequenceExpression", "sequence": "T" }
}"#;

fn err<T: serde::de::DeserializeOwned>(json: &str) -> String {
    let Err(e) = from_str::<T>(json) else {
        panic!("expected a deserialization error for {json}");
    };
    assert_eq!(e.kind(), JsonErrorKind::Data, "{e}");
    e.to_string()
}

#[test]
fn spec_example_round_trip_and_compact_output() {
    let allele: Allele = from_str(RS7412).unwrap();
    assert_eq!(
        allele.identifier().to_string(),
        "ga4gh:VA.0AePZIWZUNsUlQTamyLrjm2HWUw2opLt"
    );
    let json = to_string(&allele).unwrap();
    assert_eq!(
        json,
        r#"{"type":"Allele","location":{"type":"SequenceLocation","sequenceReference":{"type":"SequenceReference","refgetAccession":"SQ.IIB53T8CNeJJdUqzn9V_JnRtQadwWCbl"},"start":44908821,"end":44908822},"state":{"type":"LiteralSequenceExpression","sequence":"T"}}"#
    );
    let back: Allele = from_str(&json).unwrap();
    assert_eq!(back, allele);
    assert!(
        to_string_pretty(&allele)
            .unwrap()
            .contains("\n  \"location\"")
    );
}

#[test]
fn type_property_is_checked_when_present_and_defaulted_when_absent() {
    // A concretely typed object may omit `type` (as the official validation vectors do for
    // nested objects); output always carries it.
    let allele: Allele = from_str(r#"{"location":"ga4gh:SL.4t6JnYWqHwYw9WzBT_lmWBb3tLQNalkT","state":{"type":"LiteralSequenceExpression","sequence":"T"}}"#).unwrap();
    assert!(
        to_string(&allele)
            .unwrap()
            .starts_with(r#"{"type":"Allele""#)
    );
    // But a polymorphic property still needs it to dispatch.
    let msg = err::<Allele>(
        r#"{"location":"ga4gh:SL.4t6JnYWqHwYw9WzBT_lmWBb3tLQNalkT","state":{"sequence":"T"}}"#,
    );
    assert!(msg.contains("missing field `type`"), "{msg}");
    let msg = err::<SequenceLocation>(r#"{"type":"Allele","start":1,"end":2}"#);
    assert!(msg.contains("expected type \"SequenceLocation\""), "{msg}");
    let msg = err::<Variation>(r#"{"type":"Haplotype","members":[]}"#);
    assert!(
        msg.contains("unknown variant") && msg.contains("Haplotype"),
        "{msg}"
    );
}

#[test]
fn required_properties() {
    let msg = err::<Allele>(
        r#"{"type":"Allele","state":{"type":"LiteralSequenceExpression","sequence":"T"}}"#,
    );
    assert!(msg.contains("missing field `location`"), "{msg}");
    let msg = err::<SequenceReference>(r#"{"type":"SequenceReference"}"#);
    assert!(msg.contains("missing field `refgetAccession`"), "{msg}");
    let msg =
        err::<ReferenceLengthExpression>(r#"{"type":"ReferenceLengthExpression","length":3}"#);
    assert!(msg.contains("missing field `repeatSubunitLength`"), "{msg}");
    let msg = err::<SequenceLocation>(r#"{"type":"SequenceLocation"}"#);
    assert!(msg.contains("at least one of start and end"), "{msg}");
}

#[test]
fn unknown_properties_are_rejected() {
    let msg = err::<SequenceReference>(
        r#"{"type":"SequenceReference","refgetAccession":"SQ.IIB53T8CNeJJdUqzn9V_JnRtQadwWCbl","assembly":"GRCh38"}"#,
    );
    assert!(msg.contains("unknown field `assembly`"), "{msg}");
    // `digest` belongs to identifiable objects only; `expressions` to variation only.
    let msg = err::<SequenceReference>(
        r#"{"type":"SequenceReference","refgetAccession":"SQ.IIB53T8CNeJJdUqzn9V_JnRtQadwWCbl","digest":"IIB53T8CNeJJdUqzn9V_JnRtQadwWCbl"}"#,
    );
    assert!(msg.contains("unknown field `digest`"), "{msg}");
    let msg = err::<SequenceLocation>(r#"{"type":"SequenceLocation","start":1,"expressions":[]}"#);
    assert!(msg.contains("unknown field `expressions`"), "{msg}");
}

#[test]
fn invalid_enum_values() {
    let msg = err::<CopyNumberChange>(
        r#"{"type":"CopyNumberChange","location":"ga4gh:SL.4t6JnYWqHwYw9WzBT_lmWBb3tLQNalkT","copyChange":"test"}"#,
    );
    assert!(msg.contains("\"test\" is not a valid copyChange"), "{msg}");
    let msg = err::<Expression>(r#"{"syntax":"hgvs.x","value":"NC_000001.11:g.1A>T"}"#);
    assert!(msg.contains("hgvs.x"), "{msg}");
    let msg = err::<SequenceReference>(
        r#"{"type":"SequenceReference","refgetAccession":"SQ.IIB53T8CNeJJdUqzn9V_JnRtQadwWCbl","residueAlphabet":"dna"}"#,
    );
    assert!(msg.contains("residueAlphabet"), "{msg}");
}

#[test]
fn invalid_intervals_and_ranges() {
    let msg = err::<SequenceLocation>(r#"{"type":"SequenceLocation","start":10,"end":5}"#);
    assert!(msg.contains("start (10) is greater than end (5)"), "{msg}");
    let msg = err::<SequenceLocation>(r#"{"type":"SequenceLocation","start":-1,"end":5}"#);
    assert!(msg.contains("negative"), "{msg}");
    let msg = err::<SequenceLocation>(r#"{"type":"SequenceLocation","start":[null,null],"end":5}"#);
    assert!(msg.contains("at least one integer bound"), "{msg}");
    let msg = err::<SequenceLocation>(r#"{"type":"SequenceLocation","start":[7,3],"end":9}"#);
    assert!(msg.contains("minimum 7 is greater than maximum 3"), "{msg}");
    let msg = err::<SequenceLocation>(r#"{"type":"SequenceLocation","start":[1,2,3],"end":9}"#);
    assert!(msg.contains("invalid length"), "{msg}");
    let msg = err::<CopyNumberCount>(
        r#"{"type":"CopyNumberCount","location":"ga4gh:SL.4t6JnYWqHwYw9WzBT_lmWBb3tLQNalkT","copies":-2}"#,
    );
    assert!(msg.contains("copies must be non-negative"), "{msg}");
    // Circular references may wrap around the origin.
    let circular: SequenceLocation = from_str(
        r#"{"type":"SequenceLocation","sequenceReference":{"type":"SequenceReference","refgetAccession":"SQ.IIB53T8CNeJJdUqzn9V_JnRtQadwWCbl","circular":true},"start":16000,"end":100}"#,
    )
    .unwrap();
    assert_eq!(circular.exact_interval(), Some((16000, 100)));
}

#[test]
fn malformed_identifiers_and_sequences() {
    let msg = err::<SequenceReference>(
        r#"{"type":"SequenceReference","refgetAccession":"ga4gh:SQ.IIB53T8CNeJJdUqzn9V_JnRtQadwWCbl"}"#,
    );
    assert!(msg.contains("invalid RefGet accession"), "{msg}");
    let msg = err::<Allele>(
        r#"{"type":"Allele","digest":"ga4gh:734G5mtNwe40do8F6GKuqQP4QxyjBqVp","location":"x","state":{"type":"LiteralSequenceExpression","sequence":"T"}}"#,
    );
    assert!(msg.contains("invalid sha512t24u digest"), "{msg}");
    let msg = err::<LiteralSequenceExpression>(
        r#"{"type":"LiteralSequenceExpression","sequence":"acgt"}"#,
    );
    assert!(
        msg.contains("invalid residue byte 0x61 at offset 0"),
        "{msg}"
    );
    let msg = err::<LiteralSequenceExpression>(
        r#"{"type":"LiteralSequenceExpression","sequence":"ACG T"}"#,
    );
    assert!(msg.contains("offset 3"), "{msg}");
    let msg = err::<Adjacency>(
        r#"{"type":"Adjacency","adjoinedSequences":[{"type":"SequenceLocation","start":1,"end":2},{"type":"SequenceLocation","start":5}]}"#,
    );
    assert!(msg.contains("must not define both start and end"), "{msg}");
    let msg = err::<CisPhasedBlock>(
        r#"{"type":"CisPhasedBlock","members":["ga4gh:VA.0AePZIWZUNsUlQTamyLrjm2HWUw2opLt"]}"#,
    );
    assert!(msg.contains("requires at least 2 items, got 1"), "{msg}");
}

#[test]
fn syntax_errors_are_classified() {
    let e = from_str::<Allele>("{\"type\": ").unwrap_err();
    assert_eq!(e.kind(), JsonErrorKind::Eof);
    let e = from_str::<Allele>("{\"type\": Allele}").unwrap_err();
    assert_eq!(e.kind(), JsonErrorKind::Syntax);
}

#[test]
fn optional_properties_and_metadata_round_trip() {
    let json = r#"{
      "id": "ga4gh:VA.0AePZIWZUNsUlQTamyLrjm2HWUw2opLt",
      "type": "Allele",
      "digest": "0AePZIWZUNsUlQTamyLrjm2HWUw2opLt",
      "name": "rs7412",
      "description": "APOE e2",
      "aliases": ["rs7412-T"],
      "extensions": [{"name": "clinvar", "value": {"vcv": 17864, "stars": [1, 2]}, "description": "ClinVar"}],
      "expressions": [{"syntax": "hgvs.g", "value": "NC_000019.10:g.44908822C>T", "syntax_version": "21.0"}],
      "location": {
        "type": "SequenceLocation",
        "id": "NC_000019.10:44908821-44908822",
        "sequenceReference": {
          "type": "SequenceReference",
          "id": "NC_000019.10",
          "refgetAccession": "SQ.IIB53T8CNeJJdUqzn9V_JnRtQadwWCbl",
          "residueAlphabet": "na",
          "moleculeType": "genomic",
          "circular": false,
          "sequence": "ACGT"
        },
        "start": 44908821,
        "end": 44908822,
        "sequence": "C"
      },
      "state": {"type": "ReferenceLengthExpression", "length": [2, null], "repeatSubunitLength": 1, "sequence": "TT"}
    }"#;
    let allele: Allele = from_str(json).unwrap();
    assert_eq!(
        allele.id(),
        Some("ga4gh:VA.0AePZIWZUNsUlQTamyLrjm2HWUw2opLt")
    );
    assert_eq!(allele.name(), Some("rs7412"));
    assert_eq!(allele.aliases(), ["rs7412-T"]);
    assert_eq!(allele.extensions()[0].name, "clinvar");
    assert_eq!(allele.expressions()[0].syntax(), Syntax::HgvsG);
    assert_eq!(allele.expressions()[0].syntax_version(), Some("21.0"));
    let loc = allele.sequence_location().unwrap();
    assert_eq!(loc.id(), Some("NC_000019.10:44908821-44908822"));
    assert_eq!(loc.sequence().unwrap().as_str(), "C");
    let reference = loc.inline_sequence_reference().unwrap();
    assert_eq!(reference.molecule_type(), Some(MoleculeType::Genomic));
    assert_eq!(reference.circular(), Some(false));
    let rle = allele.state().as_reference_length().unwrap();
    assert_eq!(rle.length(), IntOrRange::Range(Range::at_least(2).unwrap()));

    let out = to_string(&allele).unwrap();
    let back: Allele = from_str(&out).unwrap();
    assert_eq!(back, allele);
    // Every property survived.
    for key in [
        "\"digest\"",
        "\"aliases\"",
        "\"extensions\"",
        "\"expressions\"",
        "\"syntax_version\"",
        "\"moleculeType\"",
        "\"circular\":false",
        "[2,null]",
    ] {
        assert!(out.contains(key), "{key} missing from {out}");
    }
    // The carried `digest` is decorative and not trusted: the computed digest reflects the
    // actual content (an RLE state here, so it differs from the SNV digest carried in `digest`).
    assert_ne!(allele.digest().as_str(), "0AePZIWZUNsUlQTamyLrjm2HWUw2opLt");
    assert_eq!(allele.digest(), back.digest());
    assert_eq!(
        allele.meta().unwrap().digest.unwrap().as_str(),
        "0AePZIWZUNsUlQTamyLrjm2HWUw2opLt"
    );
}

#[test]
fn iri_references_in_place_of_objects() {
    let allele: Allele = from_str(
        r#"{"type":"Allele","location":"ga4gh:SL.wIlaGykfwHIpPY2Fcxtbx4TINbbODFVz","state":{"type":"LiteralSequenceExpression","sequence":"T"}}"#,
    )
    .unwrap();
    assert!(allele.sequence_location().is_none());
    assert_eq!(
        allele.location().as_iri().unwrap().as_str(),
        "ga4gh:SL.wIlaGykfwHIpPY2Fcxtbx4TINbbODFVz"
    );
    assert_eq!(
        allele.identifier().to_string(),
        "ga4gh:VA.0AePZIWZUNsUlQTamyLrjm2HWUw2opLt"
    );
    let json = to_string(&allele).unwrap();
    assert!(json.contains(r#""location":"ga4gh:SL.wIlaGykfwHIpPY2Fcxtbx4TINbbODFVz""#));
    let block: CisPhasedBlock = from_str(
        r#"{"type":"CisPhasedBlock","members":["ga4gh:VA.0AePZIWZUNsUlQTamyLrjm2HWUw2opLt","locations.json#/1"]}"#,
    )
    .unwrap();
    assert_eq!(block.members().len(), 2);
}

#[test]
fn polymorphic_parsing_with_type_in_any_position() {
    for json in [
        r#"{"type":"Allele","id":"a","location":"ga4gh:SL.wIlaGykfwHIpPY2Fcxtbx4TINbbODFVz","state":{"sequence":"T","type":"LiteralSequenceExpression"}}"#,
        r#"{"id":"a","location":"ga4gh:SL.wIlaGykfwHIpPY2Fcxtbx4TINbbODFVz","state":{"sequence":"T","type":"LiteralSequenceExpression"},"type":"Allele"}"#,
        r#"{"location":{"start":44908821,"end":44908822,"sequenceReference":{"refgetAccession":"SQ.IIB53T8CNeJJdUqzn9V_JnRtQadwWCbl","type":"SequenceReference"},"type":"SequenceLocation"},"id":"a","state":{"sequence":"T","type":"LiteralSequenceExpression"},"type":"Allele"}"#,
    ] {
        let v: Variation = from_str(json).unwrap();
        assert_eq!(v.type_name(), "Allele");
        assert_eq!(
            v.identifier().to_string(),
            "ga4gh:VA.0AePZIWZUNsUlQTamyLrjm2HWUw2opLt"
        );
        assert_eq!(v.as_allele().unwrap().id(), Some("a"));
    }
    let msg = err::<Variation>(r#"{"location":"x"}"#);
    assert!(msg.contains("missing field `type`"), "{msg}");
    let loc: Location =
        from_str(r#"{"type":"SequenceLocation","start":[null,5],"end":9}"#).unwrap();
    assert_eq!(loc.type_name(), "SequenceLocation");
}

#[test]
fn every_variation_class_round_trips() {
    let reference = || SequenceReference::parse("SQ.IIB53T8CNeJJdUqzn9V_JnRtQadwWCbl").unwrap();
    let loc = SequenceLocation::new(reference(), 10, 20).unwrap();
    let snv = Allele::new(loc.clone(), SequenceExpression::literal("A").unwrap());
    let offset =
        SequenceOffsetLocation::new(reference(), 100, AnchorOrientation::Left, 12, 13).unwrap();
    let relative = RelativeAllele::new(
        RelativeSequenceLocation::new(loc.clone(), offset),
        SequenceExpression::literal("A").unwrap(),
        SequenceExpression::literal("T").unwrap(),
    );
    let block = CisPhasedBlock::new(vec![
        snv.clone().into(),
        Iri::new("ga4gh:VA.0AePZIWZUNsUlQTamyLrjm2HWUw2opLt").into(),
    ])
    .unwrap();
    let adjacency = Adjacency::new(
        SequenceLocation::ending_at(reference(), 123).unwrap(),
        SequenceLocation::starting_at(reference(), 456).unwrap(),
    )
    .unwrap()
    .with_linker(LengthExpression::new(20_000).unwrap())
    .with_homology(false);
    let terminus = Terminus::new(SequenceLocation::ending_at(reference(), 123).unwrap());
    let molecule = DerivativeMolecule::new(vec![
        terminus.clone().into(),
        TraversalBlock::new(adjacency.clone(), Orientation::ReverseComplement).into(),
        snv.clone().into(),
        Iri::new("components.json#/3").into(),
    ])
    .unwrap()
    .with_circular(false);
    let count = CopyNumberCount::new(loc.clone(), Range::bounded(3, 5).unwrap()).unwrap();
    let change = CopyNumberChange::new(loc, CopyChange::HighLevelGain);

    let all: Vec<Variation> = vec![
        snv.into(),
        relative.into(),
        block.into(),
        adjacency.into(),
        terminus.into(),
        molecule.into(),
        count.into(),
        change.into(),
    ];
    for v in &all {
        let json = to_string(v).unwrap();
        let back: Variation = from_str(&json).unwrap();
        assert_eq!(&back, v, "{json}");
        assert_eq!(back.identifier(), v.identifier());
        assert!(json.contains(&format!("\"type\":\"{}\"", v.type_name())));
    }
    let json = to_string(&all).unwrap();
    let back: Vec<Variation> = from_str(&json).unwrap();
    assert_eq!(back, all);
}
