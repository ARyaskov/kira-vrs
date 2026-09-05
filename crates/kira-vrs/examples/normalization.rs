//! Normalization: the same biological variant written three different ways (VCF-style with an
//! anchor base, left-aligned, right-aligned) collapses to one fully-justified VRS allele with
//! one identifier, and a deletion in a microsatellite becomes a compact reference-length
//! expression.
//!
//! ```text
//! cargo run --example normalization
//! ```

use kira_vrs::json;
use kira_vrs::normalize::{InMemorySequenceProvider, NormalizeOptions, normalize_allele};
use kira_vrs::prelude::*;

fn main() -> Result<(), kira_vrs::Error> {
    // A tiny synthetic reference. Real pipelines implement `SequenceProvider` over SeqRepo,
    // a RefGet server or an indexed FASTA.
    //            0         1         2
    //            0123456789012345678901234
    let genome = b"GGACTTTTTTAGCACACACAGTC";
    let mut provider = InMemorySequenceProvider::new();
    let accession = provider.insert_sequence(genome)?;
    let reference = SequenceReference::new(accession);
    let options = NormalizeOptions::default();

    let allele = |start: i64, end: i64, alt: &str| -> Result<Allele, kira_vrs::Error> {
        Ok(Allele::new(
            SequenceLocation::new(reference.clone(), start, end)?,
            SequenceExpression::literal(alt)?,
        ))
    };

    // One extra T in the T6 homopolymer, expressed three ways.
    let vcf_style = allele(3, 4, "CT")?; // REF "C" ALT "CT" at 1-based position 4
    let left_aligned = allele(4, 4, "T")?; // insertion before the first T
    let right_aligned = allele(10, 10, "T")?; // insertion after the last T

    println!("input                       -> normalized");
    let mut ids = Vec::new();
    for a in [&vcf_style, &left_aligned, &right_aligned] {
        let n = normalize_allele(a, &provider, &options)?;
        println!("{:<27} -> {}", describe(a), describe(&n));
        ids.push(n.identifier());
    }
    assert!(ids.iter().all(|i| i == &ids[0]));
    println!("all three share {}\n", ids[0]);

    // A CA-dinucleotide deletion in the (CA)4 microsatellite: the normalized allele spans the
    // whole repeat and the state is a ReferenceLengthExpression (6 residues left of a
    // 2-residue repeat unit).
    let deletion = allele(12, 14, "")?;
    let n = normalize_allele(&deletion, &provider, &options)?;
    println!("{:<27} -> {}", describe(&deletion), describe(&n));
    println!("{}", json::to_string_pretty(&n)?);
    Ok(())
}

fn describe(a: &Allele) -> String {
    let (s, e) = a.sequence_location().unwrap().exact_interval().unwrap();
    match a.state() {
        SequenceExpression::Literal(l) => format!("[{s}, {e}) {:?}", l.sequence().as_str()),
        SequenceExpression::ReferenceLength(r) => format!(
            "[{s}, {e}) RLE(length={}, unit={})",
            r.length(),
            r.repeat_subunit_length()
        ),
        SequenceExpression::Length(l) => format!("[{s}, {e}) length {:?}", l.length()),
    }
}
