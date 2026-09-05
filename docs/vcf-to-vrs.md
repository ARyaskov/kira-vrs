# From VCF to VRS — a conceptual mapping

This page is for people who know VCF/BCF and are new to VRS.

## What VRS is (and is not)

* **VCF ≠ VRS.** VCF is a *file format* for calls: one line per site, alleles as REF/ALT
  strings relative to a named contig, 1-based coordinates, with an anchor base for indels.
* **VRS ≠ a JSON file format.** VRS has a JSON Schema, but the schema is only the interchange
  encoding. VRS is a *semantic model*: what a variant *is* (a state at a location on an
  identified sequence), with conventions that make two independent parties produce the same
  object and the same identifier for the same variant.
* **VRS = domain model + conventions.** Classes (`Allele`, `SequenceLocation`,
  `CisPhasedBlock`, `Adjacency`, `CopyNumberCount`, …), an identifier algorithm, and
  normalization rules.

## The mapping

```
VCF record                          VRS
──────────────────────────────────  ──────────────────────────────────────────────
CHROM  "chr19" / "NC_000019.10"  →  SequenceReference { refgetAccession: "SQ.IIB53T8…" }
                                    (RefGet accession = sha512t24u of the sequence; translate
                                     with a sequence service, e.g. SeqRepo; keep the name in `id`)
POS    44908822 (1-based)        →  inter-residue interval [POS-1, POS-1+len(REF))
REF    C                         →  reference sequence at the location (optionally carried in
                                     SequenceLocation.sequence; not part of the identity)
ALT    T                         →  LiteralSequenceExpression { sequence: "T" }
(record)                         →  Allele { location, state }
normalization (left-align etc.)  →  VRS fully-justified normalization → canonical Allele
ID / stable key                  →  ga4gh:VA.<digest>, computed from the canonical Allele
```

Concretely, `NC_000019.10:g.44908822C>T` (rs7412) becomes:

```rust
use kira_vrs::prelude::*;

let chr19 = SequenceReference::parse("SQ.IIB53T8CNeJJdUqzn9V_JnRtQadwWCbl")?.with_id("NC_000019.10");
let location = SequenceLocation::new(chr19, 44_908_821, 44_908_822)?;   // [POS-1, POS-1+|REF|)
let allele = Allele::new(location, SequenceExpression::literal("T")?);  // ALT
assert_eq!(allele.identifier().to_string(), "ga4gh:VA.0AePZIWZUNsUlQTamyLrjm2HWUw2opLt");
```

and, as JSON:

```json
{
  "type": "Allele",
  "location": {
    "type": "SequenceLocation",
    "sequenceReference": { "type": "SequenceReference", "refgetAccession": "SQ.IIB53T8CNeJJdUqzn9V_JnRtQadwWCbl" },
    "start": 44908821,
    "end": 44908822
  },
  "state": { "type": "LiteralSequenceExpression", "sequence": "T" }
}
```

## Details that trip people up

**Coordinates.** VRS uses inter-residue coordinates: positions *between* residues, starting at
0. The residue at 1-based position `p` occupies `[p-1, p)`. An insertion between residues `p`
and `p+1` is the empty interval `[p, p)`. This is why VRS needs no anchor base.

**Indels and the anchor base.** VCF writes a deletion of `G` after `A` as `REF=AG ALT=A`. In
VRS the same variant is `[pos, pos+2)` with state `A`, and after normalization it is the
interval covering just the `G` (or the whole repeat run it lies in) with an empty or
reference-length state. Both inputs normalize to the same allele; you do not have to strip the
anchor yourself, but you do need the reference sequence.

**Multiallelic sites.** One VRS `Allele` per ALT. Phased genotypes across records are a
`CisPhasedBlock` of alleles; VRS itself does not model samples or genotype calls.

**Symbolic alleles.** `<DEL>`, `<DUP>`, `<CNV>` with `END`/`SVLEN` map to `CopyNumberCount` /
`CopyNumberChange` over a `SequenceLocation`; breakends (`BND`) map to `Adjacency` with two
half-open locations (defined by `end` = sequence extends left, by `start` = extends right)
and an optional linker; assembled rearrangements to `DerivativeMolecule`.

**Contig names.** VRS refuses conventional accessions in `refgetAccession`. The accession is
the digest of the sequence, so it is the same for `chr19`, `19`, `NC_000019.10` and
`CM000681.2` — and different for every patch release that changes the sequence. Translation
requires a sequence service; the `SequenceProvider` trait is where kira-vrs meets it, and the
same service supplies the bases normalization needs.

**Reference allele.** A record with `ALT` equal to `REF` (or a `.` ALT you wish to represent)
normalizes to a `ReferenceLengthExpression` whose length equals the location length.

## Why computed identifiers

`ga4gh:VA.0AePZIWZUNsUlQTamyLrjm2HWUw2opLt` is derived from the normalized allele's content.
Any lab, any pipeline, any language produces the same string for the same variant without a
registry, so identifiers can be joined across resources, generated offline in a clinical
setting, and used as database keys. They are only as stable as the classes they are built
from: identifiers of trial-use classes are stable within a major version; draft classes may
change.

## Where a `kira-vcf` crate will plug in

```
VCF/BCF record ─▶ (contig → RefgetAccession, POS/REF/ALT → typed Allele)
              ─▶ normalize_allele(&allele, &sequence_provider, &opts)
              ─▶ allele.identifier()  /  kira_vrs::json::to_writer(...)
```

No JSON or string allocation is needed on that path: `SequenceString::from_bytes` takes the
record's bytes, `RefgetAccession` is a copyable 35-byte value, and coordinates are integers.
