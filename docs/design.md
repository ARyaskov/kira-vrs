# kira-vrs design

This document records the engineering decisions behind kira-vrs: which specification revision
is implemented, how the crate is structured, how the domain model differs from the JSON
schema, how digests and normalization are implemented, and where the implementation
deliberately deviates from (or extends) the upstream text. It is the place to look before
changing anything that affects computed identifiers.

## 1. Specification revision

| Item | Value |
|---|---|
| Specification | GA4GH Variation Representation Specification (VRS) **2.1.0** |
| Repository | <https://github.com/ga4gh/vrs> |
| Git tag / commit | `2.1.0` / `cf33bfa7618011087655d5a5898e518c9d96bcdb` (2026-09-01) |
| Core dependency | gkm-core **1.2.0**, `91abbb7d0f8f05a183303853c121abd76b8b765a` |
| JSON Schema `$id` base | `https://w3id.org/ga4gh/schema/vrs/2.1.0/json/` |
| Reference implementation consulted | vrs-python `main` (2.4.0-a4 / `87ffb9d`, branch `vrs/2.1.0`) |

2.1.0 is a released (non-ballot) version: the tag exists, `vrs-source.yaml` declares
`$id …/vrs/2.1.0/…`, and the release notes describe it as the 2.1.0 release. The constants are
recorded in machine-readable form in `crates/kira-vrs/src/spec.rs` and
`crates/kira-vrs/Cargo.toml` (`[package.metadata.vrs]`), and the upstream test material is
vendored at that exact commit in `crates/kira-vrs-validation/upstream/` (see `REVISION.md`
there). Nothing is fetched at build time.

### What is 2.0, what is 2.1, what is draft, what is ours

VRS annotates every class with a GA4GH maturity level. kira-vrs implements the complete 2.1.0
schema and labels each class (`kira_vrs::spec::CLASS_MATURITY`, cross-checked against the
schemas by a test):

| Level | Classes |
|---|---|
| Stable VRS 2.0 functionality (trial use, unchanged in 2.1) | `Allele`, `CisPhasedBlock`, `Adjacency`, `CopyNumberCount`, `SequenceLocation`, `SequenceReference`, `LiteralSequenceExpression`, `ReferenceLengthExpression`, `Expression`, `Range`, `sequenceString` |
| Changed in VRS 2.1 (minor-version change) | Allele normalization now selects the **smallest** repeat-subunit factor for reference-derived ambiguous insertions (PR #700). This changes the normalized form, hence the digest, of affected alleles relative to 2.0. |
| New in 2.1 at **draft** maturity | `RelativeAllele`, `RelativeSequenceLocation`, `SequenceOffsetLocation`, `AnchorOrientation`, the relative-allele normalization rule |
| Pre-existing **draft** classes | `Terminus`, `DerivativeMolecule`, `TraversalBlock`, `CopyNumberChange`, `LengthExpression`, `Adjacency.homology` |
| Implementation-specific extensions (not VRS) | `digest::legacy` (VRS 1.3 digests for migration; matches the `ga4gh_1_3_*` validation vectors), `IriOr` resolution of `ga4gh:SQ.…` IRIs during normalization, `NormalizeOptions::rle_sequence_limit`, the in-memory `SequenceProvider`, memory-layout choices |

Draft classes are fully implemented (types, JSON, digests, validation vectors) but carry no
stability promise across VRS patch releases; their computed identifiers may change upstream.

## 2. Crate layout

```
kira-vrs/
├── crates/kira-vrs/             the library (model, digest, json, normalize, spec)
│   ├── benches/vrs.rs           Criterion benchmarks
│   ├── examples/                simple_allele, serialization, computed_identifier,
│   │                            normalization, validation
│   └── tests/                   json, normalize, properties (proptest), layout, allocations
├── crates/kira-vrs-validation/  upstream compatibility harness (not published)
│   ├── upstream/                vendored VRS 2.1.0 vectors, examples, schemas
│   └── tests/                   vectors, examples, JSON-Schema conformance
├── fuzz/                        cargo-fuzz targets (excluded from the workspace)
├── docs/                        this file, normalization.md, vcf-to-vrs.md, spec-revision.md
└── scripts/sync-upstream.sh     re-vendor upstream at an exact ref
```

The prompt proposed separate `core` / `serde` / `digest` / `normalize` crates. They were not
created: the digest layer needs the model, normalization needs the digest layer (member
ordering), and serde is the interchange format of the model itself — there is no real
dependency boundary to place a crate at, and splitting would only add version coordination.
The one boundary that exists — *library* versus *conformance harness with vendored fixtures
and a heavy JSON-Schema dev-dependency* — is a crate boundary. Fuzzing lives outside the
workspace because it needs nightly.

## 3. Domain model versus JSON

```
JSON text ──serde wire structs──▶ validating constructors ──▶ typed model
typed model ──hand-written Serialize──▶ JSON text
typed model ──DigestSerialize (RFC 8785 writer)──▶ bytes ──sha512t24u──▶ identifier
```

Principles:

* **Invariants at construction.** `SequenceLocation::new` rejects `start > end` on linear
  references and negative coordinates; `Range::new` rejects `[null, null]` and inverted
  bounds; `CisPhasedBlock::new` needs two members; `Adjacency::new` rejects adjoined
  sequences with both `start` and `end`; value sets are enums. JSON input goes through the
  same constructors, so a value of a model type is schema-valid by construction.
* **Enums for `oneOf`.** `SequenceExpression`, `Location`, `Variation`,
  `MolecularVariation`, `SystemicVariation`, `DerivativeComponent`; `IriOr<T>` for
  "IRI or object" properties; `IntOrRange` for "integer or Range".
* **Decorative metadata behind one pointer.** Every entity holds `Option<Box<Meta>>` with
  `id`, `name`, `description`, `aliases`, `extensions`, the carried `digest` and
  `expressions`. Objects without metadata cost 8 bytes for it.
* **No `serde_json::Value` on the hot path.** The only dynamic JSON in the model is
  `Extension.value`, which the schema defines as arbitrary JSON.
* **`type` is required on output, tolerated on input.** The JSON Schema marks `type`
  required on every class, but the official validation vectors themselves omit it on nested
  objects whose class is fixed by the containing property (the `RelativeAllele` vector's
  `relativeLocation`), because the reference implementation defaults it. kira-vrs therefore
  accepts a missing `type` where the class is already known, rejects a wrong one, and always
  emits it. Polymorphic properties still require it for dispatch.

### Memory layout

Measured on x86-64 (`tests/layout.rs` guards the bounds):

| Type | Size | Notes |
|---|---|---|
| `SequenceString` | 24 B | ≤ 22 residues inline, longer boxed; every SNV/most indels allocate nothing |
| `RefgetAccession` | 35 B | fixed array incl. `SQ.`; comparable and copyable |
| `Digest` | 32 B | base64url characters, so ordering = the order VRS sorts by |
| `Range` | 16 B | unbounded sides encoded with reserved `i64` extremes |
| `Option<IntOrRange>` | 24 B | niche-optimised |
| `SequenceReference` | 72 B | accession inline, annotations as small enums |
| `SequenceLocation` / `Location` | 152 B | reference inline; relative variant boxed |
| `SequenceExpression` | 64 B | |
| `Allele` | 224 B | zero heap allocations for an SNV with inline reference |
| `Variation` | 224 B | rare structural classes boxed, so the union costs no more than an `Allele` |
| `Adjacency` | 384 B | two half-open locations |

Allocation counts on the hot paths (`tests/allocations.rs`, counting global allocator): SNV
construction 0; `identifier()` 1 (the serialization buffer) and `identifier_with()` 0; JSON
deserialization of an SNV 0, as `Allele` or as `Variation`; JSON serialization 2;
normalization 1 (SNV) to 3 (small indel: reference window, seed, expanded alternate).

A cohort of ten million SNV alleles held as `Vec<Allele>` is therefore 2.2 GB with no
per-element heap traffic; a future columnar store (see §8) will do far better, and the model
is deliberately free of interior pointers that would prevent it.

## 4. JSON layer

* Serialization is hand-written per class (`json/ser.rs`): property order is `type`, `id`,
  `digest`, class properties, then decorative properties; absent optionals are omitted.
* Deserialization uses private wire structs generated by a macro (`json/de.rs`) with
  `deny_unknown_fields` (= `additionalProperties: false`) and per-class metadata field sets
  (`digest` only on identifiable classes, `expressions` only on variation classes).
* Polymorphic dispatch (`json/tagged.rs`) reads keys until `type`, buffers only the entries
  before it, and streams the remainder to the concrete deserializer — so canonical VRS JSON
  (`type` first or second) is parsed without buffering the object, unlike serde's
  internally-tagged enums.
* Conformance is proven by re-serializing every validation vector and example and validating
  the output against the official JSON Schemas with the `jsonschema` crate
  (`kira-vrs-validation/tests/json_schema.rs`).

## 5. Computed identifiers

`digest/serialize.rs` writes, per class, exactly the `ga4gh.inherent` properties in Unicode
code-point order (fixed at compile time; a test compares the tables with the schemas). Nested
identifiable objects are serialized into the tail of the same buffer, hashed and replaced by
their digest — one allocation per identifier. `sha512t24u` is SHA-512 (the `sha2` crate),
24 bytes, hand-rolled base64url (24 bytes → 32 chars is trivial and allocation-free).

### Deviation from the prose: `null` inherent properties

The *Digest Serialization* prose says implementations "filter out fields with null values".
The normative validation vectors do the opposite: an `Adjacency` without a linker serializes
as `{"adjoinedSequences":[…],"linker":null,"type":"Adjacency"}`, a `SequenceLocation`
without `sequenceReference` includes `"sequenceReference":null`, and an adjoined sequence
without `end` includes `"end":null` — this was verified numerically against the expected
digests (`elmvUghL59i1XrD-Y7cwS__tBR6EEA98`, `VJIUKfuj7QCxPI-bplNjh5bv2Y8nkvW7`, …) before
writing any code, and it is what vrs-python's `_ValueObject.ga4gh_serialize` produces.
kira-vrs follows the vectors and the reference implementation: **every inherent property is
always emitted, `null` when absent.** Interoperable identifiers matter more than the prose.

Other details worth knowing:

* IRIs of the form `ga4gh:<prefix>.<digest>` are serialized as the bare digest; other IRIs
  verbatim (gkm-core `iriReference.ga4gh_serialize`).
* `CisPhasedBlock.members` (`ordered: false`) is sorted by digest; `adjoinedSequences` and
  `components` (`ordered: true`) are not.
* `TraversalBlock` and `SequenceOffsetLocation` are value objects (inherent keys but no
  prefix): they are inlined in their parent's serialization, with their own nested
  identifiable objects (the adjacency, the base location) replaced by digests.
* Integers are written as exact decimals. RFC 8785 mandates ECMAScript formatting, which is
  identical for every integer below 2^53; the reference implementation also writes exact
  integers.
* The carried `digest`/`id` properties are never trusted: identifiers are always recomputed.

## 6. Normalization

See `docs/normalization.md` for the algorithms, complexity, edge cases and the exact
relationship to the specification text. Summary of decisions:

* Allele normalization follows VRS 2.1.0 step by step, including the **smallest-factor**
  rule. vrs-python (through 2.4.0-a4, including its `vrs/2.1.0` branch) still picks the
  greatest factor (`_factor_gen` yields descending); kira-vrs follows the 2.1.0 text and tests
  the difference explicitly (`smallest_repeat_subunit_factor_per_vrs_2_1`). All other
  vrs-python normalization test cases are reproduced with their recorded reference bases.
* The reference-derived check is the semantic one: the alternate must equal the cyclic
  extension of the first `d` reference residues truncated to its length — exactly what
  `ReferenceLengthExpression` denotes and what `expand_reference_length_expression` computes.
* Reference bases are fetched as a window (interval + context) and extended lazily during
  rolling; providers with partial segments are supported through an exact-range fallback.
* Adjacency normalization implements the ordering convention (forward first, accessions
  ascending by RefGet accession, coordinates ascending) as a choice between the two orderings
  of the adjoined sequences; orientation is not changed when a literal linker is present,
  because reversing it would require reverse complementation that the specification does
  not define. Homology expansion is "to be described" upstream and is not implemented.
* Relative-allele normalization normalizes the base representation and applies state-only
  changes; when the base location moves, the mapped offsets would have to be re-derived from
  the transcript alignment the object does not carry, so the input is returned. The anchor
  selection rule is exposed as `preferred_anchor` for callers that have both candidates.

## 7. Errors

`error.rs` defines one enum per failure domain — `CoordinateError`, `SequenceStringError`,
`IdentifierError`, `ModelError`, `SequenceError` (provider), `NormalizeError`,
`LegacyDigestError`, `JsonError` (with a `kind()` classification) and `UnsupportedError` —
plus the umbrella `Error`. Digests and identifiers are infallible by construction. No user
input is ever validated with a panic.

## 8. Performance and future integration

The design targets the pipeline `VCF/BCF record → typed Allele → normalized Allele →
identifier` with no JSON in between:

* `SequenceString::from_bytes`, `RefgetAccession` values and integer coordinates let a VCF
  reader build alleles from its own byte slices;
* `SequenceProvider` is the only external dependency of normalization and can wrap an
  indexed FASTA, SeqRepo or RefGet with borrowed (`Cow`) returns;
* identifiers are computed into a caller-reusable buffer with one allocation;
* JSON is produced directly from the model (no `Value`).

Benchmarks (`cargo bench -p kira-vrs`) cover construction, JSON in both directions, canonical
serialization, digest/identifier and the normalization cases; numbers are in the README.
SIMD/`memchr` were not introduced: the byte loops in trimming and rolling are bounded by the
width of the ambiguity region and the profile is dominated by SHA-512, so there is nothing
measurable to gain yet.

Future crates (`kira-vcf`, `kira-gvcf`, `kira-bcf`, `kira-gfa`, `kira-gbz`, an indexed
variation store) will depend on this crate's model and traits only; nothing in the public API
requires JSON, owned strings or heap-allocated sequences at the boundary.
