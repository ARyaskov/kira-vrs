//! Compatibility harness for the official GA4GH VRS validation suite.
//!
//! This crate vendors, under `upstream/`, the language-neutral test material of the pinned
//! VRS release (see `upstream/REVISION.md`):
//!
//! * `validation/models.json` — per-class vectors giving the expected digest serialization,
//!   digest and identifier (`ga4gh_serialize` / `ga4gh_digest` / `ga4gh_identify`, plus the
//!   VRS 1.3 variants);
//! * `validation/functions.json` — `sha512t24u` vectors;
//! * `examples/*.json` — the example objects from the specification, with
//!   `test_definitions.yaml` transcribed into [`EXAMPLES`];
//! * `schema/vrs/*` and `schema/gkm-core/*` — the generated JSON Schemas.
//!
//! The library part exposes typed loaders; the actual assertions live in `tests/`.

use std::path::{Path, PathBuf};

use serde::Deserialize;
use serde_json::Value;

/// Root of the vendored upstream material.
pub fn upstream_dir() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("upstream")
}

/// One `models.yaml` vector.
#[derive(Debug, Clone, Deserialize)]
pub struct ModelVector {
    /// Optional human-readable name.
    #[serde(default)]
    pub name: Option<String>,
    /// The input object (JSON, as written upstream).
    #[serde(rename = "in")]
    pub input: Value,
    /// Expected outputs keyed by function name; `null` means "not applicable / not
    /// identifiable".
    pub out: serde_json::Map<String, Value>,
}

/// Load `validation/models.json` grouped by class name, in file order.
pub fn model_vectors() -> Vec<(String, Vec<ModelVector>)> {
    let path = upstream_dir().join("validation/models.json");
    let text = std::fs::read_to_string(&path).expect("read models.json");
    let map: serde_json::Map<String, Value> =
        serde_json::from_str(&text).expect("parse models.json");
    map.into_iter()
        .map(|(class, vectors)| {
            let vectors: Vec<ModelVector> = serde_json::from_value(vectors).expect("vector shape");
            (class, vectors)
        })
        .collect()
}

/// One `functions.yaml` vector for `sha512t24u`.
#[derive(Debug, Clone, Deserialize)]
pub struct Sha512t24uVector {
    /// Input blob (UTF-8 text).
    #[serde(rename = "in")]
    pub input: Sha512t24uInput,
    /// Expected digest.
    pub out: String,
}

/// Input of a `sha512t24u` vector.
#[derive(Debug, Clone, Deserialize)]
pub struct Sha512t24uInput {
    /// The bytes to hash, as text.
    pub blob: String,
}

/// Load the `sha512t24u` vectors.
pub fn sha512t24u_vectors() -> Vec<Sha512t24uVector> {
    let path = upstream_dir().join("validation/functions.json");
    let text = std::fs::read_to_string(&path).expect("read functions.json");
    let map: serde_json::Map<String, Value> = serde_json::from_str(&text).expect("parse");
    serde_json::from_value(map["sha512t24u"].clone()).expect("vector shape")
}

/// An upstream example (`tests/test_definitions.yaml`, transcribed).
#[derive(Debug, Clone, Copy)]
pub struct Example {
    /// File name under `upstream/examples/`.
    pub file: &'static str,
    /// The VRS class the example is validated against.
    pub class: &'static str,
    /// Whether upstream expects schema validation to fail.
    pub should_fail: bool,
}

/// The upstream example definitions.
pub const EXAMPLES: &[Example] = &[
    Example {
        file: "simple_breakpoint.json",
        class: "Adjacency",
        should_fail: false,
    },
    Example {
        file: "revcomp_breakpoint.json",
        class: "Adjacency",
        should_fail: false,
    },
    Example {
        file: "terminal_breakend.json",
        class: "Terminus",
        should_fail: false,
    },
    Example {
        file: "sequence_homology.json",
        class: "Adjacency",
        should_fail: false,
    },
    Example {
        file: "precise_linker.json",
        class: "Adjacency",
        should_fail: false,
    },
    Example {
        file: "ambiguous_linker.json",
        class: "Adjacency",
        should_fail: false,
    },
    Example {
        file: "sv_derivative_molecule.json",
        class: "DerivativeMolecule",
        should_fail: false,
    },
    Example {
        file: "simple_haplotype.json",
        class: "CisPhasedBlock",
        should_fail: false,
    },
    Example {
        file: "SPDI_contraction.json",
        class: "Allele",
        should_fail: false,
    },
    Example {
        file: "SPDI_expansion.json",
        class: "Allele",
        should_fail: false,
    },
    Example {
        file: "invalid_adjacency.json",
        class: "Adjacency",
        should_fail: true,
    },
];

/// Read an example file as text.
pub fn example_text(file: &str) -> String {
    std::fs::read_to_string(upstream_dir().join("examples").join(file)).expect("read example")
}

/// Read a JSON Schema by class name (`vrs` or `gkm-core`).
pub fn schema(module: &str, class: &str) -> Value {
    let path = upstream_dir().join("schema").join(module).join(class);
    let text =
        std::fs::read_to_string(&path).unwrap_or_else(|e| panic!("read {}: {e}", path.display()));
    serde_json::from_str(&text).expect("parse schema")
}

/// All schemas of a module as `(class, schema)`.
pub fn schemas(module: &str) -> Vec<(String, Value)> {
    let dir = upstream_dir().join("schema").join(module);
    let mut out: Vec<(String, Value)> = std::fs::read_dir(&dir)
        .expect("schema dir")
        .map(|e| {
            let e = e.expect("dir entry");
            let name = e.file_name().to_string_lossy().into_owned();
            let text = std::fs::read_to_string(e.path()).expect("read schema");
            (name, serde_json::from_str(&text).expect("parse schema"))
        })
        .collect();
    out.sort_by(|a, b| a.0.cmp(&b.0));
    out
}
