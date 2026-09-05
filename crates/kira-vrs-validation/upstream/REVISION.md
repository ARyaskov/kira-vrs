# Vendored upstream revision

| Repository | Ref | Commit | Date |
|---|---|---|---|
| https://github.com/ga4gh/vrs | `2.1.0` | `cf33bfa7618011087655d5a5898e518c9d96bcdb` | 2026-09-01 |
| https://github.com/ga4gh/gkm-core (submodule) | `1.2.0` | `91abbb7d0f8f05a183303853c121abd76b8b765a` | — |

Files under this directory are copied verbatim from those revisions and are licensed under
the Apache License 2.0 (see `LICENSE-APACHE-2.0`):

| Path | Upstream source |
|---|---|
| `validation/models.json`, `validation/functions.json` | `validation/models.yaml`, `validation/functions.yaml`, converted to JSON (the pre-generated `validation/json/*.json` upstream was stale and lacked the 2.1 relative-location classes) |
| `examples/*.json` | `examples/json/*.json` |
| `examples/test_definitions.yaml` | `tests/test_definitions.yaml` |
| `schema/vrs/*` | `schema/vrs/json/*` |
| `schema/gkm-core/*` | `submodules/gkm-core/schema/gkm-core/json/*` |

Regenerate with `scripts/sync-upstream.sh 2.1.0`.
