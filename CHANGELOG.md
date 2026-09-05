# Changelog

All notable changes to this project are documented here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/) and the project uses
[Semantic Versioning](https://semver.org/).

## [Unreleased]

## [0.1.0] — 2026-09-05

First release. Implements GA4GH VRS **2.1.0** (`ga4gh/vrs` tag `2.1.0`, commit
`cf33bfa7618011087655d5a5898e518c9d96bcdb`) with gkm-core 1.2.0.

### Added

- Typed model of every VRS 2.1.0 class: `Allele`, `RelativeAllele` (draft), `CisPhasedBlock`,
  `Adjacency`, `Terminus` (draft), `DerivativeMolecule` / `TraversalBlock` (draft),
  `CopyNumberCount`, `CopyNumberChange` (draft), `SequenceLocation`,
  `RelativeSequenceLocation` / `SequenceOffsetLocation` (draft), `SequenceReference`,
  `LiteralSequenceExpression`, `ReferenceLengthExpression`, `LengthExpression` (draft),
  `Expression`, `Extension`, `Range`, `sequenceString`, value sets, `IriOr` references and the
  `Variation` / `MolecularVariation` / `SystemicVariation` / `Location` unions, with
  invariants enforced at construction.
- JSON serialization and deserialization conforming to the official JSON Schemas (verified
  with the `jsonschema` crate), with a non-buffering `type` dispatcher for polymorphic
  properties.
- RFC 8785 digest serialization, `sha512t24u`, GA4GH computed identifiers for all identifiable
  classes, plus VRS 1.3 legacy identifiers for migration.
- Normalization: fully-justified allele normalization with reference-length encoding (VRS 2.1
  smallest-factor rule), `CisPhasedBlock` member normalization and ordering, `Adjacency`
  orientation convention, `RelativeAllele` base normalization and anchor selection;
  `SequenceProvider` trait with an in-memory implementation.
- `kira-vrs-validation`: the official validation vectors (24 vectors, 15 classes, 70+
  expectations including VRS 1.3), the specification examples, and JSON-Schema conformance,
  all passing.
- Unit, integration, property (proptest) and layout tests; cargo-fuzz targets; Criterion
  benchmarks; GitHub Actions CI on Linux, macOS (ARM64) and Windows with the pinned 1.95.0
  toolchain.
- Documentation: README, `docs/design.md`, `docs/normalization.md`, `docs/vcf-to-vrs.md`,
  `docs/spec-revision.md`.

### Notes on identifier stability

- Identifiers of trial-use classes match the upstream vectors byte for byte.
- Alleles whose normalization involves a repeat unit that is itself periodic get the VRS 2.1
  smallest-factor encoding; vrs-python ≤ 2.4.0-a4 still uses the VRS 2.0 greatest-factor
  encoding for those (see `docs/normalization.md`).
- Draft classes carry no stability promise across VRS patch releases.
