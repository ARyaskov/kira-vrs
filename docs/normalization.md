# Normalization

*Specification reference:* VRS 2.1.0, *Conventions → Normalization*
(`docs/source/conventions/normalization.rst` at commit `cf33bfa`), release note for PR #700.

Normalization rewrites a variation into the canonical form the specification prescribes, so
that every representation of the same biological variant produces the same computed
identifier. The rules are semantic — they need the reference sequence — and apply to the
typed model, never to JSON text.

## Why the same variant has many textual forms

Take a reference `…GG TATATATA CC…` and one extra `TA`. It can be written as an insertion of
`TA` before the first `TA`, after any repeat unit, as `AT` one base to the right, or in VCF
style as `G → GTA` with an anchor base. Every one is correct and every one is a different
JSON object with a different digest. Fully-justified normalization expands the allele over
the *entire* region of ambiguity — the whole `(TA)4` run — so there is only one
representation left.

## Allele normalization (`normalize_allele`)

Inputs: an `Allele` with a `LiteralSequenceExpression` state and an inline `SequenceLocation`
whose reference is resolvable (inline `SequenceReference`, or a `ga4gh:SQ.…` IRI), plus a
`SequenceProvider`.

| Step (spec numbering) | Implementation |
|---|---|
| 0. reference allele sequence `ref = reference[start, end)`, alternate `alt = state.sequence` | fetched as one window `[start − c, end + c)` with `c = 2·max(|alt|, |ref|, 16)`; extended lazily if a roll reaches its edge |
| 1. trim common suffix, then common prefix | byte comparisons, O(min(|ref|, |alt|)) |
| 2a. both empty → reference allele | `ReferenceLengthExpression { length = repeatSubunitLength = |location| }` at the **original** location |
| 2b. both non-empty → substitution | `LiteralSequenceExpression(trimmed alt)` at the trimmed location |
| 2c. one empty → seed = the non-empty sequence | continue |
| 3a. left roll | while the base before the window equals the seed's last base (under cyclic rotation): move left |
| 3b. right roll | while the base after the window equals the seed's first base (under cyclic rotation): move right |
| 4. expand | `ref' = reference[left, right)`, `alt' = reference[left, trimmed start) + alt + reference[trimmed end, right)` |
| 5a. `ref'` empty → unambiguous insertion | `LiteralSequenceExpression(alt')` |
| 5b. deletion | `ReferenceLengthExpression { length = |alt'|, repeatSubunitLength = |seed| }` |
| 5c/5d. insertion | smallest factor `d` of `|seed|`, `d ≤ |ref'|`, such that `alt'` equals the first `d` residues of `ref'` repeated cyclically and truncated to `|alt'|` → `ReferenceLengthExpression { |alt'|, d }`; otherwise `LiteralSequenceExpression(alt')` |

Rolling is implemented with a rotation index rather than by rotating a buffer, so each step
is O(1); the whole algorithm is O(|ref| + |alt| + w) where `w` is the width of the region of
ambiguity (the length of the repeat run). Provider calls: one `sequence_length`, one
`sequence` for the window, plus one per window extension (geometric growth, so O(log w)).

### Edge cases

* **Coordinates as indefinite ranges** (`[null, x]`, `[x, null]`): normalized on the defined
  bound; the output keeps the same form. Definite ranges are returned unchanged (as the
  reference implementation does).
* **Non-literal states, IRI locations, unknown references:** returned unchanged, as the
  specification requires for objects without applicable rules. Reference-length states are
  the *output* of normalization and are already canonical; to re-normalize one, expand it with
  `expand_reference_length_expression` first.
* **Circular references:** the algorithm is undefined for wrap-around intervals; an
  `UnsupportedError` is returned rather than a wrong answer.
* **`start > end`:** `UnsupportedError`.
* **Sequence ends:** rolls stop at coordinate 0 and at the sequence length.
* **Empty reference allele** (`start == end`, empty alt): a reference allele of length 0
  (`ReferenceLengthExpression { 0, 0 }`), per step 2a.
* **Literal RLE sequence:** `NormalizeOptions::rle_sequence_limit` (default 50, as in
  vrs-python) attaches the decorative literal sequence to reference-length states up to that
  length; it never affects digests.

### The VRS 2.1 smallest-factor rule

VRS 2.0 selected the *greatest* factor of the seed length; 2.1.0 selects the *smallest* one
(release notes, PR #700). The two differ when a repeat unit is itself a repeat: inserting
`CACA` into `…CACACA…` is reference-derived with period 2 and with period 4; 2.1 encodes
`repeatSubunitLength: 2`. vrs-python (through 2.4.0-a4, `_factor_gen`) still yields the
greatest factor, so its identifiers for such alleles differ from VRS 2.1.0 identifiers.
kira-vrs follows the specification; the case is pinned by
`tests/normalize.rs::smallest_repeat_subunit_factor_per_vrs_2_1`.

### Tests and benchmarks

* `tests/normalize.rs`: every vrs-python normalization case (with the reference bases from
  its recorded SeqRepo responses), the 2.1 factor rule, shift invariance, idempotence, IRI
  references, error paths, `CisPhasedBlock`.
* `tests/properties.rs`: proptest — shift invariance for random insertions, deletions
  preserving the described molecule, idempotence.
* `fuzz/fuzz_targets/normalize_allele.rs`: the same invariants under libFuzzer.
* `benches/vrs.rs` (`normalize/*`): SNV, unambiguous indels, microsatellite indels,
  homopolymer insertion, 1 kb deletion, VCF-style indel.

## CisPhasedBlock (`normalize_cis_phased_block`)

Member alleles are normalized (members that omit `sequenceReference` borrow the block's, and
keep omitting it afterwards), then ordered by digest. The digest serialization sorts members
regardless; ordering the JSON form as well makes the block canonical everywhere.

## Adjacency (`normalize_adjacency`)

*Spec:* "1. The first of the adjoined sequences MUST have a forward orientation (location
defined by `end`). 2. The adjoined sequence accessions are equal or in ascending
lexicographical order. 3. The defined adjoined sequence coordinates are in ascending
numerical order."

The two orientations of an adjacency are the two orders of its adjoined sequences. Each order
is scored by the three criteria in sequence and the better order is returned; ties keep the
input order. Accessions are compared by RefGet accession, the only identifier that takes part
in digests. The orientation is left unchanged when the linker is a literal or
reference-length sequence (reversing the traversal would reverse-complement it, which the
specification does not describe) or when an adjoined sequence is an IRI or a relative
location. Ambiguity expansion for homologous breakpoints is marked "to be described" upstream
and is not implemented. O(1).

## RelativeAllele (`normalize_relative_allele`, `preferred_anchor`) — draft

*Spec:* normalize the base representation like an `Allele`, then choose between the
left-anchor and right-anchor representations of the mapped location: the one whose largest
offset magnitude is smaller wins, ties go to the left anchor.

The base representation is normalized with `normalize_allele`. State-only results (identity
→ reference-length encoding, trimmed substitutions at the same span) are applied. When the
base location moves, the mapped offsets depend on the transcript alignment (which side of the
gap, strand), information a `RelativeAllele` does not carry; the input is then returned
unchanged rather than guessed. `preferred_anchor(left, right)` applies the selection rule for
callers that can compute both candidates from an alignment.
