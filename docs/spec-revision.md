# Pinned specification revision

kira-vrs implements a single, exact revision of the GA4GH VRS specification. The revision is
recorded in three places that must agree (a test in `kira-vrs-validation` checks the schema
`$id` and maturity table against the vendored schemas):

| Where | What |
|---|---|
| `crates/kira-vrs/src/spec.rs` | `VRS_VERSION`, `VRS_TAG`, `VRS_REVISION`, `VRS_REVISION_DATE`, `VRS_SCHEMA_BASE`, `GKM_CORE_VERSION`, `GKM_CORE_REVISION`, `CLASS_MATURITY` |
| `crates/kira-vrs/Cargo.toml` → `[package.metadata.vrs]` | the same values in machine-readable Cargo metadata |
| `crates/kira-vrs-validation/upstream/REVISION.md` | the revision the vendored fixtures were copied from |

Current values:

```toml
[package.metadata.vrs]
spec_version      = "2.1.0"
spec_tag          = "2.1.0"
spec_repository   = "https://github.com/ga4gh/vrs"
spec_revision     = "cf33bfa7618011087655d5a5898e518c9d96bcdb"
spec_date         = "2026-09-01"
gkm_core_version  = "1.2.0"
gkm_core_repository = "https://github.com/ga4gh/gkm-core"
gkm_core_revision = "91abbb7d0f8f05a183303853c121abd76b8b765a"
```

## Updating to a new upstream revision

1. `scripts/sync-upstream.sh <tag-or-commit>` re-vendors the validation vectors, examples and
   schemas and rewrites `upstream/REVISION.md`.
2. Update the constants in `spec.rs` and `Cargo.toml`.
3. `cargo test --workspace`. Three tests are the early warning system:
   * `inherent_property_tables_match_upstream_schema` — a class gained or lost an inherent
     property (digests change: update `digest/serialize.rs`);
   * `type_prefixes_match_upstream_schema` / `maturity_table_matches_upstream_schema`;
   * `reserialized_vectors_validate_against_upstream_schema` — a property was added, removed
     or renamed (update the model, `json/ser.rs`, `json/de.rs`).
4. Read the upstream release notes for normalization changes; they are not detectable from
   the schemas. Add a test pinning any new rule.
5. Record the change in `CHANGELOG.md`, including whether any computed identifiers change.

Ballot and snapshot pre-releases are pinned by commit hash, never by branch.
