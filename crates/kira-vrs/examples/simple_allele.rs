//! Build the specification's canonical example allele — `NC_000019.10:g.44908822C>T`
//! (rs7412, APOE) — from typed parts and print its computed identifier.
//!
//! ```text
//! cargo run --example simple_allele
//! ```

use kira_vrs::prelude::*;

fn main() -> Result<(), kira_vrs::Error> {
    // VRS identifies sequences by RefGet accession, never by `NC_000019.10`. The accession is
    // the sha512t24u digest of the chromosome sequence; a sequence service (e.g. SeqRepo)
    // performs the translation.
    let chr19 = SequenceReference::parse("SQ.IIB53T8CNeJJdUqzn9V_JnRtQadwWCbl")?
        .with_id("NC_000019.10")
        .with_residue_alphabet(ResidueAlphabet::NucleicAcid);

    // Inter-residue coordinates: 1-based position 44908822 is the interval [44908821, 44908822).
    let location = SequenceLocation::new(chr19, 44_908_821, 44_908_822)?;

    let allele = Allele::new(location, SequenceExpression::literal("T")?)
        .with_expression(Expression::new(Syntax::HgvsG, "NC_000019.10:g.44908822C>T"));

    println!(
        "class:       {}",
        Variation::from(allele.clone()).type_name()
    );
    println!(
        "location:    {:?}",
        allele.sequence_location().unwrap().exact_interval()
    );
    println!("state:       {}", allele.state().sequence().unwrap());
    println!("identifier:  {}", allele.identifier());
    println!(
        "location id: {}",
        allele.sequence_location().unwrap().identifier()
    );

    assert_eq!(
        allele.identifier().to_string(),
        "ga4gh:VA.0AePZIWZUNsUlQTamyLrjm2HWUw2opLt"
    );
    Ok(())
}
