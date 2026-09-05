# Contributing

Thank you for considering a contribution. kira-vrs aims to be a *standards-compliant*
implementation first and a fast one second; both are checked mechanically.

## Setup

The repository pins Rust 1.95.0 through `rust-toolchain.toml`; `rustup` installs it on first
use. No other tools are needed for the library, tests, examples or benchmarks. Fuzzing needs a
nightly toolchain and `cargo-fuzz` (see `fuzz/README.md`); re-vendoring the specification
needs `python3` with `pyyaml` (see `scripts/sync-upstream.sh`).

## Before opening a pull request

```bash
cargo fmt --all --check
cargo clippy --workspace --all-targets --all-features -- -D warnings
cargo test --workspace --all-targets
cargo test --workspace --doc
cargo doc --workspace --no-deps          # with RUSTDOCFLAGS="-D warnings"
cargo bench -p kira-vrs --no-run         # benchmarks must still compile
```

CI runs exactly these on Linux, macOS and Windows. Clippy runs at the `pedantic` level with a
short, documented list of opt-outs in the workspace `Cargo.toml`; please do not add
`#[allow]` attributes without a comment explaining why.

## Rules of the road

* **Anything that changes a computed identifier is a compatibility change.** Digest
  serialization (`crates/kira-vrs/src/digest/serialize.rs`) and normalization
  (`crates/kira-vrs/src/normalize/`) must match the pinned specification revision and the
  vendored validation vectors; a change there needs a specification reference in the PR and
  a `CHANGELOG.md` entry stating which identifiers change.
* **The specification revision is pinned.** Do not update schemas, vectors or examples by
  hand; use `scripts/sync-upstream.sh` and follow `docs/spec-revision.md`.
* **No `unsafe`.** The library crate forbids it (`#![forbid(unsafe_code)]`) and the workspace
  denies it; the single exception is the counting allocator in the allocation test. If a
  benchmark shows a compelling need, open an issue first.
* **No new dependencies without a reason** stated in the PR (what it does that the standard
  library or an existing dependency cannot, licence, maintenance status, MSRV).
* **Every public item is documented**, and every algorithm notes its specification section,
  complexity and edge cases (`docs/normalization.md` is the model).
* **Tests accompany behaviour.** Prefer adding a case to the existing integration or
  property tests over a one-off unit test; use upstream test material when it exists.

## Reporting a compatibility problem

If kira-vrs produces an identifier different from another implementation for the same input,
please include the input JSON, both identifiers, both digest serializations (kira-vrs:
`DigestSerialize::digest_serialization`) and the other implementation's version. Differences
in normalization of ambiguous insertions between VRS 2.0 and 2.1 implementations are expected
(`docs/normalization.md`).
