//! Reference sequence access.

use std::borrow::Cow;
use std::collections::HashMap;
use std::sync::Arc;

use crate::error::SequenceError;
use crate::model::{RefgetAccession, SequenceString};

/// Source of reference sequence data, keyed by RefGet accession.
///
/// This is the *data proxy* the VRS *Required External Data* section describes. Implementations
/// may wrap an in-memory FASTA, a SeqRepo directory, a RefGet server or a database; only
/// sub-sequence retrieval and total length are needed by normalization.
///
/// Coordinates are inter-residue (0-based, half-open). Returned bytes must be valid
/// `sequenceString` residues (`[A-Z*-]`), i.e. upper case.
pub trait SequenceProvider {
    /// The residues in `[start, end)` of the sequence `accession`.
    ///
    /// # Errors
    /// [`SequenceError::UnknownSequence`] if the accession is unknown;
    /// [`SequenceError::OutOfBounds`] if the interval is not available.
    fn sequence(
        &self,
        accession: &RefgetAccession,
        start: u64,
        end: u64,
    ) -> Result<Cow<'_, [u8]>, SequenceError>;

    /// Total length of the sequence `accession`.
    ///
    /// # Errors
    /// [`SequenceError::UnknownSequence`] if the accession is unknown.
    fn sequence_length(&self, accession: &RefgetAccession) -> Result<u64, SequenceError>;
}

impl<P: SequenceProvider + ?Sized> SequenceProvider for &P {
    fn sequence(
        &self,
        accession: &RefgetAccession,
        start: u64,
        end: u64,
    ) -> Result<Cow<'_, [u8]>, SequenceError> {
        (**self).sequence(accession, start, end)
    }
    fn sequence_length(&self, accession: &RefgetAccession) -> Result<u64, SequenceError> {
        (**self).sequence_length(accession)
    }
}

impl<P: SequenceProvider + ?Sized> SequenceProvider for Box<P> {
    fn sequence(
        &self,
        accession: &RefgetAccession,
        start: u64,
        end: u64,
    ) -> Result<Cow<'_, [u8]>, SequenceError> {
        (**self).sequence(accession, start, end)
    }
    fn sequence_length(&self, accession: &RefgetAccession) -> Result<u64, SequenceError> {
        (**self).sequence_length(accession)
    }
}

impl<P: SequenceProvider + ?Sized> SequenceProvider for Arc<P> {
    fn sequence(
        &self,
        accession: &RefgetAccession,
        start: u64,
        end: u64,
    ) -> Result<Cow<'_, [u8]>, SequenceError> {
        (**self).sequence(accession, start, end)
    }
    fn sequence_length(&self, accession: &RefgetAccession) -> Result<u64, SequenceError> {
        (**self).sequence_length(accession)
    }
}

struct Stored {
    /// Coordinate of `bytes[0]` on the full sequence.
    offset: u64,
    bytes: Box<[u8]>,
    /// Length of the full sequence (may exceed `offset + bytes.len()` for partial segments).
    total_length: u64,
}

/// An in-memory [`SequenceProvider`]: whole sequences or partial segments with a known offset.
///
/// Partial segments make it practical to test normalization against a handful of bases from a
/// real chromosome without shipping the chromosome.
#[derive(Default)]
pub struct InMemorySequenceProvider {
    sequences: HashMap<RefgetAccession, Stored>,
}

impl InMemorySequenceProvider {
    /// An empty provider.
    pub fn new() -> Self {
        Self::default()
    }

    /// Store a complete sequence under `accession`.
    ///
    /// # Errors
    /// [`SequenceError::InvalidSequence`] if the bytes are not valid residues.
    pub fn insert(
        &mut self,
        accession: RefgetAccession,
        sequence: impl AsRef<[u8]>,
    ) -> Result<(), SequenceError> {
        let bytes = sequence.as_ref();
        crate::model::primitives::validate_residues(bytes)?;
        let total_length = bytes.len() as u64;
        self.sequences.insert(
            accession,
            Stored {
                offset: 0,
                bytes: bytes.into(),
                total_length,
            },
        );
        Ok(())
    }

    /// Store a complete sequence, computing its RefGet accession from the bytes.
    ///
    /// # Errors
    /// [`SequenceError::InvalidSequence`] if the bytes are not valid residues.
    pub fn insert_sequence(
        &mut self,
        sequence: impl AsRef<[u8]>,
    ) -> Result<RefgetAccession, SequenceError> {
        let accession = RefgetAccession::from_sequence(sequence.as_ref());
        self.insert(accession, sequence)?;
        Ok(accession)
    }

    /// Store the segment `[offset, offset + segment.len())` of a sequence whose total length
    /// is `total_length`.
    ///
    /// # Errors
    /// [`SequenceError::InvalidSequence`] if the bytes are not valid residues;
    /// [`SequenceError::OutOfBounds`] if the segment does not fit in `total_length`.
    pub fn insert_segment(
        &mut self,
        accession: RefgetAccession,
        offset: u64,
        segment: impl AsRef<[u8]>,
        total_length: u64,
    ) -> Result<(), SequenceError> {
        let bytes = segment.as_ref();
        crate::model::primitives::validate_residues(bytes)?;
        let end = offset + bytes.len() as u64;
        if end > total_length {
            return Err(SequenceError::OutOfBounds {
                accession: accession.to_string(),
                start: offset,
                end,
            });
        }
        self.sequences.insert(
            accession,
            Stored {
                offset,
                bytes: bytes.into(),
                total_length,
            },
        );
        Ok(())
    }

    /// Number of stored sequences.
    pub fn len(&self) -> usize {
        self.sequences.len()
    }

    /// `true` if nothing is stored.
    pub fn is_empty(&self) -> bool {
        self.sequences.is_empty()
    }

    /// The residues of a stored sequence as a [`SequenceString`] (whole sequences only).
    pub fn get(&self, accession: &RefgetAccession) -> Option<SequenceString> {
        let s = self.sequences.get(accession)?;
        (s.offset == 0 && s.bytes.len() as u64 == s.total_length)
            .then(|| SequenceString::from_valid_bytes(&s.bytes))
    }
}

impl SequenceProvider for InMemorySequenceProvider {
    fn sequence(
        &self,
        accession: &RefgetAccession,
        start: u64,
        end: u64,
    ) -> Result<Cow<'_, [u8]>, SequenceError> {
        let s = self
            .sequences
            .get(accession)
            .ok_or_else(|| SequenceError::UnknownSequence(accession.to_string()))?;
        let available_end = s.offset + s.bytes.len() as u64;
        if start > end || start < s.offset || end > available_end {
            return Err(SequenceError::OutOfBounds {
                accession: accession.to_string(),
                start,
                end,
            });
        }
        let lo = (start - s.offset) as usize;
        let hi = (end - s.offset) as usize;
        Ok(Cow::Borrowed(&s.bytes[lo..hi]))
    }

    fn sequence_length(&self, accession: &RefgetAccession) -> Result<u64, SequenceError> {
        self.sequences
            .get(accession)
            .map(|s| s.total_length)
            .ok_or_else(|| SequenceError::UnknownSequence(accession.to_string()))
    }
}
