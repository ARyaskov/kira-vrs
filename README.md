# kira-vrs

A native Rust implementation of the GA4GH **Variation Representation Specification (VRS) 2.1**:
a strongly typed domain model, standards-compatible JSON, RFC 8785 digest serialization with
`sha512t24u` computed identifiers, and the VRS normalization algorithms — built as the
foundation for high-throughput genomic variation tooling (VCF/BCF conversion, cohort indexes,
variation graphs) rather than as a JSON wrapper.

```text
genomic variation ─▶ typed VRS model ─▶ canonical (normalized) form ─▶ stable identifier
```

| | |
|---|---|
| Specification | VRS **2.1.0** (`ga4gh/vrs` tag `2.1.0`, commit `cf33bfa7…`, 2026-09-01), gkm-core 1.2.0 |
| Rust | 1.95.0 (pinned), edition 2024, stable only, `#![forbid(unsafe_code)]` in the library |
| Conformance | all official validation vectors, examples and JSON Schemas pass (`cargo test -p kira-vrs-validation`) |
| Licence | MIT (vendored GA4GH test material: Apache-2.0) |

## Why VRS, in one paragraph

VCF is a file format; VRS is a **semantic model** of what a variant *is*: the state of a
sequence at a location on an identified reference. Its conventions — RefGet accessions
instead of contig names, inter-residue coordinates, fully-justified normalization — make two
independent parties describe the same variant with the same object, and its computed
identifier (`ga4gh:VA.…`, a digest of the canonical object) gives that variant one globally
stable key without any registry. JSON is only VRS's interchange encoding; kira-vrs keeps a
compact typed model in memory and treats JSON as an import/export format. See
[`docs/vcf-to-vrs.md`](docs/vcf-to-vrs.md) if you come from VCF.

## Quick start

```toml
[dependencies]
kira-vrs = "0.1"
```

```rust
use kira_vrs::prelude::*;

fn main() -> Result<(), kira_vrs::Error> {
    // NC_000019.10:g.44908822C>T (rs7412). VRS identifies the chromosome by its RefGet
    // accession — the sha512t24u digest of the sequence — not by name.
    let chr19 = SequenceReference::parse("SQ.IIB53T8CNeJJdUqzn9V_JnRtQadwWCbl")?
        .with_id("NC_000019.10");

    // Inter-residue coordinates: 1-based position 44908822 is [44908821, 44908822).
    let location = SequenceLocation::new(chr19, 44_908_821, 44_908_822)?;
    let allele = Allele::new(location, SequenceExpression::literal("T")?);

    println!("{}", allele.identifier()); // ga4gh:VA.0AePZIWZUNsUlQTamyLrjm2HWUw2opLt

    let json = kira_vrs::json::to_string_pretty(&allele)?;
    let back: Allele = kira_vrs::json::from_str(&json)?;
    assert_eq!(back, allele);
    Ok(())
}
```

The JSON above is the specification's own example:

```json
{
  "type": "Allele",
  "location": {
    "type": "SequenceLocation",
    "sequenceReference": {
      "type": "SequenceReference",
      "id": "NC_000019.10",
      "refgetAccession": "SQ.IIB53T8CNeJJdUqzn9V_JnRtQadwWCbl"
    },
    "start": 44908821,
    "end": 44908822
  },
  "state": { "type": "LiteralSequenceExpression", "sequence": "T" }
}
```

Polymorphic documents parse into the `Variation` / `Location` unions; every class (`Allele`,
`RelativeAllele`, `CisPhasedBlock`, `Adjacency`, `Terminus`, `DerivativeMolecule`,
`CopyNumberCount`, `CopyNumberChange`, and the locations, references and sequence
expressions) has a constructor that enforces its invariants, so invalid objects cannot be
built or parsed. More in [`crates/kira-vrs/examples/`](crates/kira-vrs/examples/):
`simple_allele`, `serialization`, `computed_identifier`, `normalization`, `validation`.

## Computed identifiers

```rust
use kira_vrs::digest::{DigestSerialize, Identifiable, sha512t24u};

let blob = allele.digest_serialization();
// {"location":"wIlaGykfwHIpPY2Fcxtbx4TINbbODFVz","state":{"sequence":"T","type":"LiteralSequenceExpression"},"type":"Allele"}
assert_eq!(sha512t24u(&blob), allele.digest());
assert_eq!(allele.identifier().to_string(), "ga4gh:VA.0AePZIWZUNsUlQTamyLrjm2HWUw2opLt");
```

Only the class's *inherent* properties take part; `id`, `name`, `expressions`, extensions and
the carried `digest` never change an identifier. Nested identifiable objects are replaced by
their digests, unordered lists are sorted, and the result is canonical JSON per RFC 8785 —
written straight into one byte buffer, no `serde_json::Value` involved. `Identifiable::
identifier_with(&mut scratch)` reuses a caller buffer for allocation-free loops. VRS 1.3
identifiers are available in `kira_vrs::digest::legacy` for migrating older data.

## Normalization

```rust
use kira_vrs::normalize::{InMemorySequenceProvider, NormalizeOptions, normalize_allele};

let mut reference = InMemorySequenceProvider::new();
let acc = reference.insert_sequence(b"GGTATATATACC")?;      // any SequenceProvider works
let r = SequenceReference::new(acc);
let a = Allele::new(SequenceLocation::new(r.clone(), 2, 2)?, SequenceExpression::literal("TA")?);
let b = Allele::new(SequenceLocation::new(r.clone(), 1, 2)?, SequenceExpression::literal("GTA")?); // VCF style
let opts = NormalizeOptions::default();
assert_eq!(
    normalize_allele(&a, &reference, &opts)?.identifier(),
    normalize_allele(&b, &reference, &opts)?.identifier(),
);
```

The same insertion written before the repeat, after it, or with a VCF anchor base yields one
fully-justified allele spanning the whole `(TA)4` run, with a compact
`ReferenceLengthExpression` state. `normalize` dispatches on any `Variation`; alleles,
cis-phased blocks, adjacencies and relative alleles have rules, everything else passes
through unchanged as the specification requires. Algorithms, complexity and edge cases:
[`docs/normalization.md`](docs/normalization.md).

## Architecture

```text
             External formats
          ┌───────┬───────┬───────┐
          │ VCF   │ HGVS  │ JSON  │      (VCF/HGVS: future crates; JSON: kira_vrs::json)
          └───┬───┴───┬───┴───┬───┘
              │       │       │
              └───────┼───────┘
                      ▼
               kira_vrs::model            typed, invariant-checked, compact
                      │
          ┌───────────┼────────────┐
          ▼           ▼            ▼
   normalize      digest         validation      kira_vrs::normalize / ::digest /
   (semantic)   (RFC 8785 +      (constructors,   kira-vrs-validation (official vectors,
                 sha512t24u)      JSON Schema)     examples, schemas)
                      │
                      ▼
               ga4gh:VA.<digest>
```

* `crates/kira-vrs` — the library: `model`, `json`, `digest`, `normalize`, `spec`.
* `crates/kira-vrs-validation` — the official GA4GH validation suite, vendored at the pinned
  revision and executed as tests (vectors, examples, JSON-Schema conformance, and
  cross-checks of the inherent-property tables, prefixes and maturity levels against the
  schemas).
* `fuzz/` — cargo-fuzz targets for JSON parsing, round trips, digests and normalization.
* `docs/` — [design](docs/design.md), [normalization](docs/normalization.md),
  [VCF → VRS](docs/vcf-to-vrs.md), [pinned revision](docs/spec-revision.md).

Design in brief ([`docs/design.md`](docs/design.md)): the JSON schema is not the internal
representation. Sequences up to 22 residues are stored inline, RefGet accessions and digests
are fixed-size arrays, ranges are 16 bytes, decorative metadata sits behind one nullable
pointer, and the rarely used structural classes are boxed inside the `Variation` union. An
SNV `Allele` is 224 bytes with **zero heap allocations**; ten million of them fit in a
`Vec` in 2.2 GB. Polymorphic JSON is dispatched on `type` without buffering the object.
Every deviation from the specification prose (notably that `null` inherent properties are
serialized, as the normative vectors require) is documented there.

## Performance

Criterion, single thread, x86-64 desktop, Rust 1.95 release profile
(`cargo bench -p kira-vrs`; `target/criterion` has the reports):

| Operation | Time |
|---|---|
| construct SNV allele (parse accession, location, state) | 27 ns |
| JSON serialize SNV allele (~190 B) | 209 ns (1.1 GiB/s) |
| JSON deserialize SNV allele | 500 ns (`type` first) / 1.18 µs (`type` last, buffered) |
| JSON serialize / deserialize allele with a 10 kb insertion | 3.4 µs (2.8 GiB/s) / 3.8 µs (2.5 GiB/s) |
| canonical digest serialization (incl. nested location digest) | 400 ns |
| `identifier()` — two SHA-512 rounds | 785 ns (754 ns with a reused buffer) |
| `sha512t24u` of 128 bytes | 344 ns |
| normalize SNV / VCF-style indel | 216 ns / 218 ns |
| normalize microsatellite insertion / deletion | 432 ns / 336 ns |
| normalize 1 kb deletion | 1.96 µs |

Allocation counts (`cargo test -p kira-vrs --test allocations`, a counting global
allocator): SNV construction 0, `identifier()` 1 (the serialization buffer; 0 with
`identifier_with`), JSON deserialization of an SNV 0 (as `Allele` or as `Variation`), JSON
serialization 2 (the output string), normalization of an SNV 1 and of a small indel 3. Identifier generation is dominated by
SHA-512 itself; no SIMD or platform-specific code is used, and none is needed for the byte
loops in trimming and rolling, which are bounded by the width of the repeat region.

## Development

```bash
cargo test --workspace            # unit, integration, property, validation and doc tests
cargo clippy --workspace --all-targets -- -D warnings
cargo doc --workspace --no-deps
cargo bench -p kira-vrs
cargo run -p kira-vrs --example normalization
```

CI runs formatting, clippy (pedantic), tests, docs and the examples on Linux, macOS (ARM64)
and Windows with the pinned toolchain, plus a stable-toolchain build and a nightly fuzz
build. See [`CONTRIBUTING.md`](CONTRIBUTING.md).

## Status and scope

Version 0.1 covers the complete VRS 2.1.0 model. Trial-use classes (`Allele`,
`CisPhasedBlock`, `Adjacency`, `CopyNumberCount`, `SequenceLocation`, …) are stable;
draft classes (`RelativeAllele`, `Terminus`, `DerivativeMolecule`, `CopyNumberChange`,
`LengthExpression`, relative locations) are implemented but follow upstream draft semantics.
Not in scope for this version, by design: VCF/BCF/GVCF readers, HGVS parsing, sequence
services (implement `SequenceProvider`), and adjacency homology expansion (not yet specified
upstream). Planned crates: `kira-vcf`, `kira-gvcf`, `kira-bcf`, `kira-gfa`, `kira-gbz`, and an
indexed variation store, all building on this model.
