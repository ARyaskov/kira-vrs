//! JSON interchange conforming to the official VRS 2.1 JSON Schema.
//!
//! Every model type implements [`serde::Serialize`] and [`serde::Deserialize`], so the
//! ordinary `serde_json` entry points work. The functions in this module are thin wrappers
//! that return the crate's [`JsonError`].
//!
//! Design:
//!
//! * **Serialization** is implemented by hand on the domain types and writes straight to the
//!   serializer — no intermediate `serde_json::Value`. Absent optional properties are omitted;
//!   `type` is always emitted.
//! * **Deserialization** goes through private *wire* structs that mirror the schema exactly
//!   (`additionalProperties: false` ⇒ `deny_unknown_fields`, `required` ⇒ non-`Option`) and
//!   are then converted through the same validating constructors used by Rust callers, so
//!   parsed data satisfies every model invariant.
//! * **Polymorphic properties** (`oneOf` discriminated by `type`) are dispatched by a custom
//!   visitor that reads keys until it finds `type` and then streams the *rest* of the object to
//!   the concrete class's deserializer. Only properties that precede `type` (usually just
//!   `id`) are buffered; nothing is buffered when `type` comes first.
//! * `IriOr<T>` accepts a JSON string (IRI) or an object; `IntOrRange` accepts an integer or a
//!   two-element array.

mod de;
mod ser;
mod tagged;

pub use crate::error::{JsonError, JsonErrorKind};

/// Serialize to a compact JSON string.
///
/// # Errors
/// Only if the underlying writer fails (never for in-memory strings in practice).
pub fn to_string<T: serde::Serialize + ?Sized>(value: &T) -> Result<String, JsonError> {
    serde_json::to_string(value).map_err(JsonError)
}

/// Serialize to an indented JSON string.
///
/// # Errors
/// See [`to_string`].
pub fn to_string_pretty<T: serde::Serialize + ?Sized>(value: &T) -> Result<String, JsonError> {
    serde_json::to_string_pretty(value).map_err(JsonError)
}

/// Serialize to compact JSON bytes.
///
/// # Errors
/// See [`to_string`].
pub fn to_vec<T: serde::Serialize + ?Sized>(value: &T) -> Result<Vec<u8>, JsonError> {
    serde_json::to_vec(value).map_err(JsonError)
}

/// Serialize compact JSON into a writer.
///
/// # Errors
/// I/O errors from the writer.
pub fn to_writer<W: std::io::Write, T: serde::Serialize + ?Sized>(
    writer: W,
    value: &T,
) -> Result<(), JsonError> {
    serde_json::to_writer(writer, value).map_err(JsonError)
}

/// Deserialize from a JSON string (borrowing where possible).
///
/// # Errors
/// Syntax errors, schema violations (unknown properties, missing required properties, bad
/// `type`) and model invariant violations (negative coordinates, `start > end`, ...).
pub fn from_str<'a, T: serde::Deserialize<'a>>(s: &'a str) -> Result<T, JsonError> {
    serde_json::from_str(s).map_err(JsonError)
}

/// Deserialize from JSON bytes (borrowing where possible).
///
/// # Errors
/// See [`from_str`].
pub fn from_slice<'a, T: serde::Deserialize<'a>>(bytes: &'a [u8]) -> Result<T, JsonError> {
    serde_json::from_slice(bytes).map_err(JsonError)
}

/// Deserialize from a reader.
///
/// # Errors
/// See [`from_str`], plus I/O errors.
pub fn from_reader<R: std::io::Read, T: serde::de::DeserializeOwned>(
    reader: R,
) -> Result<T, JsonError> {
    serde_json::from_reader(reader).map_err(JsonError)
}
