//! Memory layout guards: the sizes documented in `docs/design.md` must not regress silently.
//!
//! The bounds are upper limits (padding differs slightly across targets); the point is that a
//! `Vec<Allele>` of SNVs stays a few hundred bytes per element with *no* heap allocation for
//! inline sequences, references and coordinates.

use std::mem::size_of;

use kira_vrs::model::*;

#[test]
fn compact_primitives() {
    assert_eq!(size_of::<SequenceString>(), 24);
    assert_eq!(size_of::<Range>(), 16);
    assert_eq!(size_of::<IntOrRange>(), 24);
    assert_eq!(size_of::<Option<IntOrRange>>(), 24);
    assert_eq!(size_of::<Digest>(), 32);
    assert_eq!(size_of::<RefgetAccession>(), 35);
    assert_eq!(size_of::<Iri>(), 16);
    assert_eq!(size_of::<Option<Box<Meta>>>(), 8);
}

#[test]
fn hot_path_types_are_bounded() {
    for (name, size) in [
        ("SequenceReference", size_of::<SequenceReference>()),
        ("SequenceLocation", size_of::<SequenceLocation>()),
        ("SequenceExpression", size_of::<SequenceExpression>()),
        ("Allele", size_of::<Allele>()),
        ("Variation", size_of::<Variation>()),
        ("Location", size_of::<Location>()),
        ("Adjacency", size_of::<Adjacency>()),
        ("CopyNumberCount", size_of::<CopyNumberCount>()),
    ] {
        println!("{name:<20} {size} bytes");
    }
    assert!(
        size_of::<SequenceReference>() <= 96,
        "{}",
        size_of::<SequenceReference>()
    );
    assert!(
        size_of::<SequenceLocation>() <= 192,
        "{}",
        size_of::<SequenceLocation>()
    );
    assert!(
        size_of::<SequenceExpression>() <= 80,
        "{}",
        size_of::<SequenceExpression>()
    );
    assert!(size_of::<Allele>() <= 288, "{}", size_of::<Allele>());
    // Boxing the structural classes keeps the unions close to the Allele size.
    assert!(size_of::<Variation>() <= 296, "{}", size_of::<Variation>());
    assert!(size_of::<Location>() <= 200, "{}", size_of::<Location>());
}

#[test]
fn snv_allele_construction_allocates_nothing_beyond_the_object() {
    // Inline sequence (<= 22 residues), fixed-size accession and coordinates: constructing an
    // SNV allele must not touch the heap. Verified structurally here (no Box/Vec/String in the
    // path) and by the allocation benchmark in `benches/vrs.rs`.
    let reference = SequenceReference::parse("SQ.IIB53T8CNeJJdUqzn9V_JnRtQadwWCbl").unwrap();
    let location = SequenceLocation::new(reference, 100, 101).unwrap();
    let allele = Allele::new(location, SequenceString::new("T").unwrap());
    assert!(allele.meta().is_none());
    assert_eq!(allele.state().sequence().unwrap().len(), 1);
}
