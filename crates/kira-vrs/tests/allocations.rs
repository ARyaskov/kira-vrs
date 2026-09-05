//! Allocation counts on the hot paths, measured with a counting global allocator.
//!
//! These are the numbers behind the "zero allocations for an SNV" claims in the README and
//! `docs/design.md`; the assertions are upper bounds so that the test also documents what a
//! regression looks like.
//!
//! The counting allocator is the one place in the workspace that needs `unsafe`
//! (`GlobalAlloc` is an unsafe trait); it is a thin, well-known wrapper over `System`.
//! The file is a custom (`harness = false`), single-threaded test binary so that no other
//! test's allocations are counted.
#![allow(unsafe_code)]

use std::alloc::{GlobalAlloc, Layout, System};
use std::sync::atomic::{AtomicUsize, Ordering};

use kira_vrs::digest::Identifiable;
use kira_vrs::model::*;
use kira_vrs::normalize::{InMemorySequenceProvider, NormalizeOptions, normalize_allele};

struct Counting;

static ALLOCATIONS: AtomicUsize = AtomicUsize::new(0);

// SAFETY: every call is forwarded unchanged to `System`, which upholds the `GlobalAlloc`
// contract; the counter is a relaxed atomic increment with no other side effects.
unsafe impl GlobalAlloc for Counting {
    unsafe fn alloc(&self, layout: Layout) -> *mut u8 {
        ALLOCATIONS.fetch_add(1, Ordering::Relaxed);
        unsafe { System.alloc(layout) }
    }
    unsafe fn dealloc(&self, ptr: *mut u8, layout: Layout) {
        unsafe { System.dealloc(ptr, layout) }
    }
    unsafe fn realloc(&self, ptr: *mut u8, layout: Layout, new_size: usize) -> *mut u8 {
        ALLOCATIONS.fetch_add(1, Ordering::Relaxed);
        unsafe { System.realloc(ptr, layout, new_size) }
    }
}

#[global_allocator]
static GLOBAL: Counting = Counting;

/// Run `f` and return how many allocations (including reallocations) it performed.
fn count<T>(f: impl FnOnce() -> T) -> (usize, T) {
    let before = ALLOCATIONS.load(Ordering::Relaxed);
    let out = f();
    let after = ALLOCATIONS.load(Ordering::Relaxed);
    (after - before, out)
}

/// Measure and report one operation.
fn report<T>(label: &str, f: impl FnOnce() -> T) -> (usize, T) {
    let (n, out) = count(f);
    println!("{label:<48} {n} allocation(s)");
    (n, out)
}

const CHR19: &str = "SQ.IIB53T8CNeJJdUqzn9V_JnRtQadwWCbl";

fn snv() -> Allele {
    let reference = SequenceReference::parse(CHR19).unwrap();
    Allele::new(
        SequenceLocation::new(reference, 44_908_821, 44_908_822).unwrap(),
        SequenceString::new("T").unwrap(),
    )
}

fn constructing_an_snv_allele_does_not_allocate() {
    let (n, allele) = report("construct SNV allele", snv);
    assert_eq!(n, 0, "SNV construction allocated {n} times");
    assert_eq!(allele.state().sequence().unwrap().len(), 1);

    // Up to 22 residues are inline; longer sequences take exactly one allocation.
    let (n, _) = count(|| SequenceString::new("ACGTACGTACGTACGTACGTAC").unwrap());
    assert_eq!(n, 0);
    let (n, _) = count(|| SequenceString::new("ACGTACGTACGTACGTACGTACG").unwrap());
    assert_eq!(n, 1);
}

fn identifier_allocates_once_or_never_with_a_reused_buffer() {
    let allele = snv();
    let (n, id) = report("identifier()", || allele.identifier());
    assert_eq!(
        n, 1,
        "identifier() allocated {n} times (expected the one serialization buffer)"
    );
    assert_eq!(id.to_string(), "ga4gh:VA.0AePZIWZUNsUlQTamyLrjm2HWUw2opLt");

    let mut scratch = Vec::with_capacity(256);
    let (n, id2) = report("identifier_with(scratch)", || {
        allele.identifier_with(&mut scratch)
    });
    assert_eq!(n, 0, "identifier_with() allocated {n} times");
    assert_eq!(id2, id);
}

fn json_deserialization_of_an_snv_is_allocation_free_when_type_comes_first() {
    let allele = snv();
    let json = kira_vrs::json::to_string(&allele).unwrap();
    let (n, back) = report("json::from_str::<Allele> (SNV)", || {
        kira_vrs::json::from_str::<Allele>(&json).unwrap()
    });
    assert_eq!(back, allele);
    assert_eq!(n, 0, "deserializing an SNV allocated {n} times");

    let (n, v) = report("json::from_str::<Variation> (SNV)", || {
        kira_vrs::json::from_str::<Variation>(&json).unwrap()
    });
    assert_eq!(v.identifier(), allele.identifier());
    assert_eq!(n, 0, "polymorphic deserialization allocated {n} times");

    // Serialization grows one output buffer.
    let (n, _) = report("json::to_string (SNV)", || {
        kira_vrs::json::to_string(&allele).unwrap()
    });
    assert!(n <= 3, "serialization allocated {n} times");
}

fn normalization_of_small_variants_allocates_a_bounded_number_of_times() {
    let mut provider = InMemorySequenceProvider::new();
    let seq: Vec<u8> = (0..4096usize)
        .map(|i| b"ACGTTGCAAGCTTAGC"[i % 16])
        .collect();
    let acc = provider.insert_sequence(&seq).unwrap();
    let opts = NormalizeOptions::default();
    let mk = |s: i64, e: i64, alt: &str| {
        Allele::new(
            SequenceLocation::new(SequenceReference::new(acc), s, e).unwrap(),
            SequenceExpression::literal(alt).unwrap(),
        )
    };
    for (name, allele, bound) in [
        ("snv", mk(1000, 1001, "T"), 1),
        ("insertion", mk(1000, 1000, "TTTT"), 3),
        ("deletion", mk(1000, 1004, ""), 3),
    ] {
        let (n, _) = report(&format!("normalize_allele ({name})"), || {
            normalize_allele(&allele, &provider, &opts).unwrap()
        });
        assert!(
            n <= bound,
            "normalizing {name} allocated {n} times (bound {bound})"
        );
    }
}

fn main() {
    // Warm up: first use of a code path may lazily allocate (e.g. thread-local buffers).
    let _ = snv().identifier();
    constructing_an_snv_allele_does_not_allocate();
    identifier_allocates_once_or_never_with_a_reused_buffer();
    json_deserialization_of_an_snv_is_allocation_free_when_type_comes_first();
    normalization_of_small_variants_allocates_a_bounded_number_of_times();
    println!("allocation counts within bounds");
}
