//! Fuzz allele normalization over random reference sequences and alleles: no panics,
//! idempotence, and the normalized allele must describe the same molecule as the input.
#![no_main]

use arbitrary::Arbitrary;
use kira_vrs::model::*;
use kira_vrs::normalize::{
    InMemorySequenceProvider, NormalizeOptions, expand_reference_length_expression,
    normalize_allele,
};
use libfuzzer_sys::fuzz_target;

#[derive(Arbitrary, Debug)]
struct Input {
    reference: Vec<u8>,
    start: u16,
    len: u8,
    alt: Vec<u8>,
}

fn base(b: u8) -> u8 {
    b"ACGT"[(b & 3) as usize]
}

fuzz_target!(|input: Input| {
    if input.reference.is_empty() || input.reference.len() > 4096 {
        return;
    }
    let reference: Vec<u8> = input.reference.iter().map(|b| base(*b)).collect();
    let alt: Vec<u8> = input.alt.iter().take(64).map(|b| base(*b)).collect();
    let start = usize::from(input.start) % reference.len();
    let end = (start + usize::from(input.len)).min(reference.len());

    let mut provider = InMemorySequenceProvider::new();
    let acc = provider.insert_sequence(&reference).unwrap();
    let location = SequenceLocation::new(SequenceReference::new(acc), start as i64, end as i64).unwrap();
    let allele = Allele::new(location, SequenceString::from_bytes(&alt).unwrap());
    let opts = NormalizeOptions::default();

    let n = normalize_allele(&allele, &provider, &opts).unwrap();
    assert_eq!(normalize_allele(&n, &provider, &opts).unwrap(), n, "idempotent");

    // Apply both alleles to the reference: the results must be identical molecules.
    let (ns, ne) = n.sequence_location().unwrap().exact_interval().unwrap();
    let (ns, ne) = (ns as usize, ne as usize);
    let n_alt: Vec<u8> = match n.state() {
        SequenceExpression::Literal(l) => l.sequence().as_bytes().to_vec(),
        SequenceExpression::ReferenceLength(r) => {
            expand_reference_length_expression(&reference[ns..ne], r).unwrap().as_bytes().to_vec()
        }
        SequenceExpression::Length(_) => unreachable!(),
    };
    let mut applied_in = reference[..start].to_vec();
    applied_in.extend_from_slice(&alt);
    applied_in.extend_from_slice(&reference[end..]);
    let mut applied_norm = reference[..ns].to_vec();
    applied_norm.extend_from_slice(&n_alt);
    applied_norm.extend_from_slice(&reference[ne..]);
    assert_eq!(applied_in, applied_norm, "normalized allele must describe the same molecule");
});
