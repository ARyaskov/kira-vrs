//! Criterion benchmarks for the hot paths: construction, JSON (de)serialization, digest
//! serialization, identifiers and normalization.
//!
//! Run with `cargo bench -p kira-vrs`; HTML reports land in `target/criterion`.

use std::hint::black_box;

use criterion::{BenchmarkId, Criterion, Throughput, criterion_group, criterion_main};
use kira_vrs::digest::{DigestSerialize, Identifiable};
use kira_vrs::model::*;
use kira_vrs::normalize::{InMemorySequenceProvider, NormalizeOptions, normalize_allele};

const CHR19: &str = "SQ.IIB53T8CNeJJdUqzn9V_JnRtQadwWCbl";

fn reference() -> SequenceReference {
    SequenceReference::parse(CHR19).unwrap()
}

fn snv() -> Allele {
    Allele::new(
        SequenceLocation::new(reference(), 44_908_821, 44_908_822).unwrap(),
        SequenceExpression::literal("T").unwrap(),
    )
}

fn deletion(len: usize) -> Allele {
    Allele::new(
        SequenceLocation::new(reference(), 44_908_821, 44_908_821 + len as i64).unwrap(),
        SequenceExpression::literal("").unwrap(),
    )
}

fn insertion(len: usize) -> Allele {
    let seq: String = (0..len).map(|i| ['A', 'C', 'G', 'T'][i % 4]).collect();
    Allele::new(
        SequenceLocation::new(reference(), 44_908_821, 44_908_821).unwrap(),
        SequenceExpression::literal(&seq).unwrap(),
    )
}

fn bench_construct(c: &mut Criterion) {
    let acc = RefgetAccession::parse(CHR19).unwrap();
    c.bench_function("construct/snv_allele", |b| {
        b.iter(|| {
            let loc = SequenceLocation::new(
                SequenceReference::new(black_box(acc)),
                black_box(44_908_821),
                black_box(44_908_822),
            )
            .unwrap();
            Allele::new(loc, SequenceString::new(black_box("T")).unwrap())
        });
    });
    c.bench_function("construct/refget_accession_parse", |b| {
        b.iter(|| RefgetAccession::parse(black_box(CHR19)).unwrap());
    });
}

fn bench_json(c: &mut Criterion) {
    let mut group = c.benchmark_group("json");
    let allele = snv();
    let json = kira_vrs::json::to_string(&allele).unwrap();
    group.throughput(Throughput::Bytes(json.len() as u64));
    group.bench_function("serialize/snv_allele", |b| {
        b.iter(|| kira_vrs::json::to_string(black_box(&allele)).unwrap());
    });
    group.bench_function("serialize_to_vec/snv_allele", |b| {
        b.iter(|| kira_vrs::json::to_vec(black_box(&allele)).unwrap());
    });
    group.bench_function("deserialize/snv_allele", |b| {
        b.iter(|| kira_vrs::json::from_str::<Allele>(black_box(&json)).unwrap());
    });
    group.bench_function("deserialize/snv_allele_as_variation", |b| {
        b.iter(|| kira_vrs::json::from_str::<Variation>(black_box(&json)).unwrap());
    });
    // `type` last forces the streaming dispatcher to buffer the preceding properties.
    let type_last = r#"{"location":{"start":44908821,"end":44908822,"sequenceReference":{"refgetAccession":"SQ.IIB53T8CNeJJdUqzn9V_JnRtQadwWCbl","type":"SequenceReference"},"type":"SequenceLocation"},"state":{"sequence":"T","type":"LiteralSequenceExpression"},"type":"Allele"}"#;
    group.bench_function("deserialize/snv_allele_as_variation_type_last", |b| {
        b.iter(|| kira_vrs::json::from_str::<Variation>(black_box(type_last)).unwrap());
    });
    group.finish();

    let mut group = c.benchmark_group("json/long_insertion");
    for len in [100usize, 1_000, 10_000] {
        let allele = insertion(len).with_expression(Expression::new(Syntax::Spdi, "x"));
        let json = kira_vrs::json::to_string(&allele).unwrap();
        group.throughput(Throughput::Bytes(json.len() as u64));
        group.bench_with_input(BenchmarkId::new("serialize", len), &allele, |b, a| {
            b.iter(|| kira_vrs::json::to_string(black_box(a)).unwrap());
        });
        group.bench_with_input(BenchmarkId::new("deserialize", len), &json, |b, json| {
            b.iter(|| kira_vrs::json::from_str::<Allele>(black_box(json)).unwrap());
        });
    }
    group.finish();
}

fn bench_digest(c: &mut Criterion) {
    let mut group = c.benchmark_group("digest");
    let allele = snv();
    group.bench_function("canonicalize/snv_allele", |b| {
        let mut buf = Vec::with_capacity(256);
        b.iter(|| {
            buf.clear();
            black_box(&allele).write_digest_serialization(&mut buf);
            black_box(buf.len())
        });
    });
    group.bench_function("digest/snv_allele", |b| {
        b.iter(|| black_box(&allele).digest());
    });
    group.bench_function("identifier/snv_allele", |b| {
        b.iter(|| black_box(&allele).identifier());
    });
    group.bench_function("identifier/snv_allele_reused_buffer", |b| {
        let mut scratch = Vec::with_capacity(256);
        b.iter(|| black_box(&allele).identifier_with(&mut scratch));
    });
    group.bench_function("identifier/deletion_1kb", |b| {
        let allele = deletion(1000);
        b.iter(|| black_box(&allele).identifier());
    });
    group.bench_function("identifier/sequence_location", |b| {
        let loc = allele.sequence_location().unwrap();
        b.iter(|| black_box(loc).identifier());
    });
    group.bench_function("sha512t24u/128_bytes", |b| {
        let data = [b'A'; 128];
        b.iter(|| kira_vrs::digest::sha512t24u(black_box(&data)));
    });
    for len in [1_000usize, 100_000] {
        let allele = insertion(len);
        group.bench_with_input(
            BenchmarkId::new("identifier/insertion", len),
            &allele,
            |b, a| {
                b.iter(|| black_box(a).identifier());
            },
        );
    }
    let block = CisPhasedBlock::new(
        (0..8)
            .map(|i| {
                IriOr::Object(Allele::new(
                    SequenceLocation::new(reference(), 100 + i, 101 + i).unwrap(),
                    SequenceExpression::literal("A").unwrap(),
                ))
            })
            .collect(),
    )
    .unwrap();
    group.bench_function("identifier/cis_phased_block_8", |b| {
        b.iter(|| black_box(&block).identifier());
    });
    group.finish();
}

fn bench_normalize(c: &mut Criterion) {
    let mut group = c.benchmark_group("normalize");
    let mut provider = InMemorySequenceProvider::new();
    // A synthetic 1 Mb chromosome with a CAG microsatellite and a homopolymer run.
    let mut seq: Vec<u8> = (0..1_000_000usize)
        .map(|i| b"ACGTTGCAAGCTTAGC"[i % 16])
        .collect();
    seq[500_000..500_060].copy_from_slice(&b"CAG".repeat(20));
    seq[600_000..600_030].fill(b'A');
    let acc = provider.insert_sequence(&seq).unwrap();
    let opts = NormalizeOptions::default();
    let mk = |s: i64, e: i64, alt: &str| {
        Allele::new(
            SequenceLocation::new(SequenceReference::new(acc), s, e).unwrap(),
            SequenceExpression::literal(alt).unwrap(),
        )
    };
    let cases = [
        ("snv", mk(1000, 1001, "T")),
        ("insertion_unambiguous", mk(1000, 1000, "TTTT")),
        ("deletion_unambiguous", mk(1000, 1004, "")),
        ("insertion_microsatellite", mk(500_030, 500_030, "CAG")),
        ("deletion_microsatellite", mk(500_030, 500_033, "")),
        ("insertion_homopolymer", mk(600_010, 600_010, "A")),
        ("deletion_1kb", mk(2000, 3000, "")),
        ("vcf_style_indel", mk(999, 1004, "A")),
    ];
    for (name, allele) in &cases {
        group.bench_with_input(BenchmarkId::new("allele", name), allele, |b, a| {
            b.iter(|| normalize_allele(black_box(a), &provider, &opts).unwrap());
        });
    }
    group.finish();
}

criterion_group!(
    benches,
    bench_construct,
    bench_json,
    bench_digest,
    bench_normalize
);
criterion_main!(benches);
