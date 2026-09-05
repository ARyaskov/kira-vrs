# Fuzz targets

Fuzzing uses [cargo-fuzz](https://github.com/rust-fuzz/cargo-fuzz) (libFuzzer), which needs a
nightly toolchain to *run*; the library itself is stable-only.

```bash
cargo install cargo-fuzz
cargo +nightly fuzz run json_variation        # arbitrary bytes → JSON parser
cargo +nightly fuzz run json_round_trip       # structured alleles → JSON → identifier stability
cargo +nightly fuzz run normalize_allele      # random references/alleles → normalization
cargo +nightly fuzz run digest_serialization  # canonical JSON invariants
```

Targets assert crash-freedom plus the crate's core invariants: JSON round trips are lossless,
identifiers are deterministic, digest serializations are canonical JSON, and normalization is
idempotent and preserves the molecule it describes.

This directory is excluded from the workspace so that `cargo test --workspace` stays
stable-only; it depends on `kira-vrs` by path.
