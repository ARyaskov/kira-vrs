//! Minimal RFC 8785 (JSON Canonicalization Scheme) primitives for writing digest
//! serializations straight into a byte buffer.
//!
//! Only the subset VRS needs is implemented: strings, integers, `null` and the structural
//! characters. Key ordering is the caller's responsibility (keys are emitted in compile-time
//! sorted order by the per-class serializers).

use crate::model::{IntOrRange, Range};

const HEX: &[u8; 16] = b"0123456789abcdef";

/// Write a JSON string literal with RFC 8785 §3.2.2.2 escaping.
pub(crate) fn write_str(out: &mut Vec<u8>, s: &str) {
    out.push(b'"');
    let bytes = s.as_bytes();
    let mut start = 0;
    for (i, &b) in bytes.iter().enumerate() {
        let escape: &[u8] = match b {
            b'"' => b"\\\"",
            b'\\' => b"\\\\",
            0x08 => b"\\b",
            0x0C => b"\\f",
            b'\n' => b"\\n",
            b'\r' => b"\\r",
            b'\t' => b"\\t",
            0x00..=0x1F => {
                out.extend_from_slice(&bytes[start..i]);
                out.extend_from_slice(b"\\u00");
                out.push(HEX[usize::from(b >> 4)]);
                out.push(HEX[usize::from(b & 0xF)]);
                start = i + 1;
                continue;
            }
            _ => continue,
        };
        out.extend_from_slice(&bytes[start..i]);
        out.extend_from_slice(escape);
        start = i + 1;
    }
    out.extend_from_slice(&bytes[start..]);
    out.push(b'"');
}

/// Write a quoted string known to contain only ASCII characters that never need escaping
/// (digests, RefGet accessions, residues).
#[inline]
pub(crate) fn write_quoted_ascii(out: &mut Vec<u8>, s: &[u8]) {
    debug_assert!(
        s.iter()
            .all(|b| b.is_ascii_graphic() && *b != b'"' && *b != b'\\')
    );
    out.reserve(s.len() + 2);
    out.push(b'"');
    out.extend_from_slice(s);
    out.push(b'"');
}

/// Write an integer in exact decimal form.
pub(crate) fn write_i64(out: &mut Vec<u8>, v: i64) {
    let mut buf = [0u8; 20];
    out.extend_from_slice(itoa(&mut buf, v));
}

/// Write `null`.
#[inline]
pub(crate) fn write_null(out: &mut Vec<u8>) {
    out.extend_from_slice(b"null");
}

/// Write an optional integer-or-range (`null` when absent).
pub(crate) fn write_opt_int_or_range(out: &mut Vec<u8>, v: Option<IntOrRange>) {
    match v {
        None => write_null(out),
        Some(v) => write_int_or_range(out, v),
    }
}

/// Write an integer or a `[min, max]` range (with `null` for unbounded sides).
pub(crate) fn write_int_or_range(out: &mut Vec<u8>, v: IntOrRange) {
    match v {
        IntOrRange::Int(i) => write_i64(out, i),
        IntOrRange::Range(r) => write_range(out, r),
    }
}

/// Write a `[min, max]` range.
pub(crate) fn write_range(out: &mut Vec<u8>, r: Range) {
    out.push(b'[');
    match r.min() {
        Some(v) => write_i64(out, v),
        None => write_null(out),
    }
    out.push(b',');
    match r.max() {
        Some(v) => write_i64(out, v),
        None => write_null(out),
    }
    out.push(b']');
}

/// Write `"key":` (key must be a plain ASCII identifier, which all VRS property names are).
#[inline]
pub(crate) fn write_key(out: &mut Vec<u8>, key: &str) {
    out.push(b'"');
    out.extend_from_slice(key.as_bytes());
    out.extend_from_slice(b"\":");
}

/// Format an `i64` into `buf`, returning the used slice.
fn itoa(buf: &mut [u8; 20], v: i64) -> &[u8] {
    let negative = v < 0;
    let mut n = v.unsigned_abs();
    let mut i = buf.len();
    loop {
        i -= 1;
        buf[i] = b'0' + (n % 10) as u8;
        n /= 10;
        if n == 0 {
            break;
        }
    }
    if negative {
        i -= 1;
        buf[i] = b'-';
    }
    &buf[i..]
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn escapes_per_rfc8785() {
        // a "b \c <LF> <SOH> €
        let input = String::from_utf8(vec![
            b'a', b'"', b'b', b'\\', b'c', b'\n', 0x01, 0xE2, 0x82, 0xAC,
        ])
        .unwrap();
        let mut out = Vec::new();
        write_str(&mut out, &input);
        let mut expected = Vec::new();
        expected.extend_from_slice(b"\"a\\\"b\\\\c\\n\\u0001");
        expected.extend_from_slice(&[0xE2, 0x82, 0xAC, b'"']);
        assert_eq!(out, expected);
    }

    #[test]
    fn integers() {
        for v in [0i64, 7, -7, 44_908_822, i64::MAX - 1, i64::MIN + 1] {
            let mut out = Vec::new();
            write_i64(&mut out, v);
            assert_eq!(out, v.to_string().into_bytes());
        }
    }

    #[test]
    fn ranges() {
        let mut out = Vec::new();
        write_range(&mut out, Range::at_least(3).unwrap());
        assert_eq!(out, b"[3,null]");
    }
}
