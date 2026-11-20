//! Reads from BAM files (aligned or unaligned) through pasta.

use std::io::{self, Read};

use pasta::bam;

use crate::io::format::BAM_MAGIC;

/// The BAM 4-bit base codes, indexed by code (SAM specification §4.2).
const BASES: &[u8; 16] = b"=ACMGRSVTWYHKDBN";

/// Yields the primary records of a BAM stream as reads in their original orientation.
pub(crate) struct BamReads<R> {
    reader: bam::io::Reader<R>,
    record: bam::Record,
    records_read: u64,
}

impl<R: Read> BamReads<R> {
    /// Checks the magic number and skips the header of the decompressed BAM stream `inner`.
    ///
    /// The SAM header text is discarded unparsed, so unusual headers cannot fail the read.
    pub(crate) fn new(inner: R) -> io::Result<Self> {
        let mut reader = bam::io::Reader::from(inner);
        let mut header = reader.header_reader();
        if header.read_magic_number()? != BAM_MAGIC {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "invalid BAM magic number",
            ));
        }
        header.raw_sam_header_reader()?.discard_to_end()?;
        header.read_reference_sequences()?;
        Ok(Self {
            reader,
            record: bam::Record::default(),
            records_read: 0,
        })
    }

    /// Returns the number of BAM records read so far, skipped ones included.
    pub(crate) fn records_read(&self) -> u64 {
        self.records_read
    }

    /// Reads the next primary record into `name` and `seq`, returning `Ok(false)` at the end.
    ///
    /// Secondary and supplementary records are skipped. A record flagged as reverse
    /// complemented is reverse complemented back, as `samtools fastq` does.
    pub(crate) fn read(&mut self, name: &mut Vec<u8>, seq: &mut Vec<u8>) -> io::Result<bool> {
        loop {
            if self.reader.read_record(&mut self.record)? == 0 {
                return Ok(false);
            }
            self.records_read += 1;
            let flags = self.record.flags();
            if flags.is_secondary() || flags.is_supplementary() {
                continue;
            }
            name.clear();
            if let Some(record_name) = self.record.name() {
                name.extend_from_slice(record_name);
            }
            let sequence = self.record.sequence();
            decode_sequence(sequence.as_bytes(), sequence.len(), seq)?;
            if flags.is_reverse_complemented() {
                reverse_complement_in_place(seq);
            }
            return Ok(true);
        }
    }
}

/// Replaces `out` with the `len` bases packed two per byte in `packed`.
fn decode_sequence(packed: &[u8], len: usize, out: &mut Vec<u8>) -> io::Result<()> {
    if packed.len() < len.div_ceil(2) {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!(
                "BAM sequence of {len} bases is stored in only {} bytes",
                packed.len()
            ),
        ));
    }
    out.clear();
    out.reserve(len);
    for &pair in &packed[..len.div_ceil(2)] {
        out.push(BASES[usize::from(pair >> 4)]);
        out.push(BASES[usize::from(pair & 0x0f)]);
    }
    out.truncate(len);
    Ok(())
}

fn reverse_complement_in_place(seq: &mut [u8]) {
    seq.reverse();
    for base in seq {
        *base = crate::seq::complement(*base);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn decodes_packed_bases() {
        let mut out = Vec::new();
        // = A C M G R S V T W Y H K D B N, then an odd-length tail.
        decode_sequence(
            &[0x01, 0x23, 0x45, 0x67, 0x89, 0xab, 0xcd, 0xef],
            16,
            &mut out,
        )
        .unwrap();
        assert_eq!(out, BASES);
        decode_sequence(&[0x12, 0x80], 3, &mut out).unwrap();
        assert_eq!(out, b"ACT");
        assert!(decode_sequence(&[0x12], 3, &mut out).is_err());
    }

    #[test]
    fn reverse_complements_in_place() {
        let mut seq = b"AACGTN".to_vec();
        reverse_complement_in_place(&mut seq);
        assert_eq!(seq, b"NACGTT");
    }
}
