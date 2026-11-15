//! Reading every supported input format from tiny files.

use std::{
    error::Error,
    fs::{self, File},
    io::Write,
    num::NonZero,
    path::{Path, PathBuf},
};

use flate2::{Compression as GzLevel, write::GzEncoder};
use pasta::{
    bam, bgzf,
    core::Position,
    sam::{
        self,
        alignment::{
            RecordBuf,
            io::Write as _,
            record::{
                Flags,
                cigar::{Op, op::Kind},
            },
            record_buf::{Cigar, Sequence},
        },
        header::record::value::{Map, map::ReferenceSequence},
    },
};
use picoasm_io::{
    Compression, InputFormat, InputOptions, ReadBatch, ReadError, ReadRecord, ReadSource,
    RecordFormat, detect_format, fastq,
};
use picoasm_seq::{ReadId, reverse_complement};
use tempfile::TempDir;

type TestResult = Result<(), Box<dyn Error>>;

const FASTA: &[u8] = b">f1 first\nACGT\nAC\r\n>f2\n>f3 third read\nNNacgtRY\n";
const FASTQ: &[u8] = b"@q1 desc\nACGTT\n+\nIIIII\n@q2\nGG\n+q2\n!!\n";

fn write_file(dir: &TempDir, name: &str, data: &[u8]) -> Result<PathBuf, Box<dyn Error>> {
    let path = dir.path().join(name);
    fs::write(&path, data)?;
    Ok(path)
}

fn gzip(data: &[u8]) -> Result<Vec<u8>, Box<dyn Error>> {
    let mut encoder = GzEncoder::new(Vec::new(), GzLevel::fast());
    encoder.write_all(data)?;
    Ok(encoder.finish()?)
}

fn bgzip(data: &[u8]) -> Result<Vec<u8>, Box<dyn Error>> {
    let mut writer = bgzf::io::Writer::new(Vec::new());
    writer.write_all(data)?;
    Ok(writer.finish()?)
}

/// Returns `(name, seq)` of every read, checking that ids count up from 0.
fn read_all(source: &ReadSource) -> Result<Vec<(String, Vec<u8>)>, ReadError> {
    let mut reads = Vec::new();
    for (i, record) in source.open()?.enumerate() {
        let record = record?;
        assert_eq!(
            record.id,
            ReadId(u32::try_from(i).expect("test ids fit u32"))
        );
        reads.push((
            String::from_utf8(record.name).expect("test names are UTF-8"),
            record.seq,
        ));
    }
    Ok(reads)
}

fn names(reads: &[(String, Vec<u8>)]) -> Vec<&str> {
    reads.iter().map(|(name, _)| name.as_str()).collect()
}

fn expected_fasta() -> Vec<(String, Vec<u8>)> {
    vec![
        ("f1".to_owned(), b"ACGTAC".to_vec()),
        ("f2".to_owned(), Vec::new()),
        ("f3".to_owned(), b"NNacgtRY".to_vec()),
    ]
}

fn expected_fastq() -> Vec<(String, Vec<u8>)> {
    vec![
        ("q1".to_owned(), b"ACGTT".to_vec()),
        ("q2".to_owned(), b"GG".to_vec()),
    ]
}

fn unmapped(name: &str, seq: &[u8], flags: Flags) -> RecordBuf {
    RecordBuf::builder()
        .set_name(name)
        .set_flags(flags | Flags::UNMAPPED)
        .set_sequence(Sequence::from(seq.to_vec()))
        .build()
}

fn mapped(name: &str, stored_seq: &[u8], flags: Flags) -> Result<RecordBuf, Box<dyn Error>> {
    Ok(RecordBuf::builder()
        .set_name(name)
        .set_flags(flags)
        .set_reference_sequence_id(0)
        .set_alignment_start(Position::try_from(1)?)
        .set_cigar(Cigar::from(vec![Op::new(Kind::Match, stored_seq.len())]))
        .set_sequence(Sequence::from(stored_seq.to_vec()))
        .build())
}

fn write_bam(path: &Path, records: &[RecordBuf]) -> TestResult {
    let header = sam::Header::builder()
        .add_reference_sequence(
            "chr1",
            Map::<ReferenceSequence>::new(NonZero::new(1000).ok_or("zero length")?),
        )
        .add_comment("picoasm-io test")
        .build();
    let mut writer = bam::io::Writer::new(File::create(path)?);
    writer.write_header(&header)?;
    for record in records {
        writer.write_alignment_record(&header, record)?;
    }
    writer.try_finish()?;
    Ok(())
}

#[test]
fn reads_plain_fasta_and_fastq() -> TestResult {
    let dir = TempDir::new()?;
    let fasta = write_file(&dir, "reads.txt", FASTA)?;
    let fastq = write_file(&dir, "reads.fasta", FASTQ)?; // Extension deliberately misleading.
    assert_eq!(read_all(&ReadSource::new([&fasta]))?, expected_fasta());
    assert_eq!(read_all(&ReadSource::new([&fastq]))?, expected_fastq());
    assert_eq!(
        detect_format(&fastq)?,
        InputFormat::new(Compression::None, RecordFormat::Fastq)
    );
    Ok(())
}

#[test]
fn reads_gzip_including_multiple_members() -> TestResult {
    let dir = TempDir::new()?;
    let mut data = gzip(&FASTQ[..20])?;
    data.extend(gzip(&FASTQ[20..])?);
    let path = write_file(&dir, "reads.fq.gz", &data)?;
    assert_eq!(read_all(&ReadSource::new([&path]))?, expected_fastq());
    assert_eq!(detect_format(&path)?.compression, Compression::Gzip);
    Ok(())
}

#[test]
fn reads_bgzf_with_several_workers() -> TestResult {
    let dir = TempDir::new()?;
    let path = write_file(&dir, "reads.fa.bgz", &bgzip(FASTA)?)?;
    let options = InputOptions::new().with_bgzf_workers(NonZero::new(2).ok_or("zero")?);
    let source = ReadSource::new([&path]).with_options(options);
    assert_eq!(read_all(&source)?, expected_fasta());
    assert_eq!(
        source.detect_formats()?,
        [InputFormat::new(Compression::Bgzf, RecordFormat::Fasta)]
    );
    Ok(())
}

#[test]
fn reads_unaligned_bam() -> TestResult {
    let dir = TempDir::new()?;
    let path = dir.path().join("movie.subreads");
    write_bam(
        &path,
        &[
            unmapped("m/1/ccs", b"ACGTNACGT", Flags::empty()),
            unmapped("m/2/ccs", b"", Flags::empty()),
            unmapped("m/3/ccs", b"TTTAC", Flags::empty()),
        ],
    )?;
    assert_eq!(
        detect_format(&path)?,
        InputFormat::new(Compression::Bgzf, RecordFormat::Bam)
    );
    let reads = read_all(&ReadSource::new([&path]))?;
    assert_eq!(names(&reads), ["m/1/ccs", "m/2/ccs", "m/3/ccs"]);
    assert_eq!(reads[0].1, b"ACGTNACGT");
    assert!(reads[1].1.is_empty());
    assert_eq!(reads[2].1, b"TTTAC");
    Ok(())
}

#[test]
fn aligned_bam_skips_non_primary_and_restores_orientation() -> TestResult {
    let dir = TempDir::new()?;
    let path = dir.path().join("aligned.bam");
    let original = b"AAACGTCAGG";
    let stored = reverse_complement(original);
    write_bam(
        &path,
        &[
            mapped("fwd", b"ACGTA", Flags::empty())?,
            mapped("rev", &stored, Flags::REVERSE_COMPLEMENTED)?,
            mapped(
                "rev",
                &stored[..4],
                Flags::SECONDARY | Flags::REVERSE_COMPLEMENTED,
            )?,
            mapped("fwd", b"ACG", Flags::SUPPLEMENTARY)?,
            unmapped("unm", b"GGT", Flags::empty()),
        ],
    )?;
    let reads = read_all(&ReadSource::new([&path]))?;
    assert_eq!(names(&reads), ["fwd", "rev", "unm"]);
    assert_eq!(reads[0].1, b"ACGTA");
    assert_eq!(reads[1].1, original);
    assert_eq!(reads[2].1, b"GGT");
    Ok(())
}

#[test]
fn ids_continue_across_files_and_passes() -> TestResult {
    let dir = TempDir::new()?;
    let fasta = write_file(&dir, "a.fa", FASTA)?;
    let fastq = write_file(&dir, "b.fq.gz", &gzip(FASTQ)?)?;
    let empty = write_file(&dir, "c.fa", b"")?;
    let bam_path = dir.path().join("d.bam");
    write_bam(&bam_path, &[unmapped("b1", b"ACGT", Flags::empty())])?;
    let source = ReadSource::new([&fasta, &fastq, &empty, &bam_path, &fasta]);

    let first = read_all(&source)?;
    assert_eq!(
        names(&first),
        ["f1", "f2", "f3", "q1", "q2", "b1", "f1", "f2", "f3"]
    );
    // A second pass sees exactly the same reads with the same ids.
    assert_eq!(read_all(&source)?, first);
    assert_eq!(source.detect_formats()?[2].records, RecordFormat::Empty);
    Ok(())
}

#[test]
fn batches_are_contiguous_and_respect_the_base_budget() -> TestResult {
    let dir = TempDir::new()?;
    let fasta = write_file(&dir, "a.fa", FASTA)?;
    let fastq = write_file(&dir, "b.fq", FASTQ)?;
    let source = ReadSource::new([&fasta, &fastq, &fasta]);
    let expected: Vec<ReadRecord> = source.open()?.collect::<Result<_, _>>()?;

    for max_bases in [0, 1, 6, 7, 1000] {
        let mut stream = source.open()?;
        let mut batch = ReadBatch::with_max_bases(max_bases);
        let mut seen = Vec::new();
        while stream.next_batch(&mut batch)? {
            let first = batch.first_id().ok_or("empty batch")?;
            assert_eq!(first.index(), seen.len(), "max_bases {max_bases}");
            assert_eq!(
                batch.bases(),
                batch.iter().map(|r| r.seq.len()).sum::<usize>()
            );
            // Only the last read may take the batch over its budget.
            let before_last = batch.bases() - batch.records()[batch.len() - 1].seq.len();
            assert!(batch.len() == 1 || before_last < max_bases, "{batch:?}");
            seen.extend(batch.records().iter().cloned());
        }
        assert_eq!(seen, expected, "max_bases {max_bases}");
        assert!(batch.is_empty());
        assert_eq!(stream.reads_read(), 8);
    }
    Ok(())
}

#[test]
fn malformed_fastq_reports_file_record_and_line() -> TestResult {
    let dir = TempDir::new()?;
    let path = write_file(
        &dir,
        "bad.fq.gz",
        &gzip(b"@q1\nACGT\n+\nIIII\n@q2\nACGT\n+\nII")?,
    )?;
    let mut stream = ReadSource::new([&path]).open()?;
    assert!(stream.next().ok_or("missing first read")?.is_ok());
    let error = stream.next().ok_or("missing error")?.unwrap_err();
    match &error {
        ReadError::Fastq {
            path: error_path,
            record,
            source,
        } => {
            assert_eq!(error_path, &path);
            assert_eq!(*record, 2);
            assert!(matches!(
                source,
                fastq::ParseError::QualityLength { line: 8, .. }
            ));
        }
        other => panic!("unexpected error {other:?}"),
    }
    assert!(error.to_string().contains("bad.fq.gz"), "{error}");
    // The stream is finished after an error.
    assert!(stream.next().is_none());

    let truncated = write_file(&dir, "truncated.fq", b"@q1\nACGT\n+\nIIII\n@q2\nAC")?;
    let error = read_all(&ReadSource::new([&truncated])).unwrap_err();
    assert!(
        matches!(
            error,
            ReadError::Fastq {
                record: 2,
                source: fastq::ParseError::Truncated { line: 7, .. },
                ..
            }
        ),
        "{error:?}"
    );
    Ok(())
}

#[test]
fn rejects_unknown_formats_missing_files_and_bad_bam() -> TestResult {
    let dir = TempDir::new()?;
    let text = write_file(&dir, "notes.fa", b"ACGT\n")?;
    let error = ReadSource::new([&text]).open().unwrap_err();
    assert!(
        matches!(
            error,
            ReadError::UnrecognizedFormat {
                first_byte: b'A',
                ..
            }
        ),
        "{error:?}"
    );
    assert_eq!(error.path(), text);

    let missing = dir.path().join("missing.fq");
    let good = write_file(&dir, "good.fq", FASTQ)?;
    let source = ReadSource::new([&good, &missing]);
    assert!(matches!(
        source.detect_formats(),
        Err(ReadError::Open { .. })
    ));
    // The missing second file is only noticed when the stream reaches it.
    let mut stream = source.open()?;
    assert!(stream.next().ok_or("missing read")?.is_ok());
    assert!(stream.next().ok_or("missing read")?.is_ok());
    assert!(matches!(stream.next(), Some(Err(ReadError::Open { .. }))));

    let not_bam = write_file(&dir, "fake.bam", &bgzip(b"BAD!")?)?;
    let error = ReadSource::new([&not_bam]).open().unwrap_err();
    assert!(
        matches!(
            error,
            ReadError::UnrecognizedFormat {
                first_byte: b'B',
                ..
            }
        ),
        "{error:?}"
    );

    // The magic number is right but the header is cut short.
    let cut_header = write_file(&dir, "cut_header.bam", &bgzip(b"BAM\x01\x10\x00")?)?;
    let error = ReadSource::new([&cut_header]).open().unwrap_err();
    assert!(matches!(error, ReadError::BamHeader { .. }), "{error:?}");

    let bam_path = dir.path().join("cut.bam");
    write_bam(&bam_path, &[unmapped("r", b"ACGTACGT", Flags::empty())])?;
    let mut raw = Vec::new();
    std::io::Read::read_to_end(&mut bgzf::io::Reader::new(File::open(&bam_path)?), &mut raw)?;
    raw.truncate(raw.len() - 3);
    let cut = write_file(&dir, "cut2.bam", &bgzip(&raw)?)?;
    let error = read_all(&ReadSource::new([&cut])).unwrap_err();
    assert!(
        matches!(error, ReadError::Bam { record: 1, .. }),
        "{error:?}"
    );
    Ok(())
}

#[test]
fn empty_source_yields_nothing() -> TestResult {
    let source = ReadSource::new(Vec::<PathBuf>::new());
    let mut stream = source.open()?;
    let mut batch = ReadBatch::new();
    assert!(!stream.next_batch(&mut batch)?);
    assert!(stream.next().is_none());
    Ok(())
}

#[test]
fn detect_formats_rejects_text_starting_with_b_and_sam() -> TestResult {
    let dir = TempDir::new()?;
    let text = write_file(&dir, "notes.bam", b"Bonjour, this is not a BAM file\n")?;
    let error = ReadSource::new([&text]).detect_formats().unwrap_err();
    assert!(
        matches!(
            error,
            ReadError::UnrecognizedFormat {
                first_byte: b'B',
                ..
            }
        ),
        "{error:?}"
    );

    let sam = write_file(
        &dir,
        "reads.sam.gz",
        &gzip(b"@HD\tVN:1.6\tSO:unknown\n@PG\tID:x\nr1\t4\t*\t0\t0\t*\t*\t0\t0\tACGT\t*\n")?,
    )?;
    let error = ReadSource::new([&sam]).detect_formats().unwrap_err();
    assert!(
        matches!(error, ReadError::UnsupportedFormat { format: "SAM", .. }),
        "{error:?}"
    );
    assert!(error.to_string().contains("SAM"), "{error}");

    let cram = write_file(&dir, "reads.cram", b"CRAM\x03\x00rest")?;
    let error = ReadSource::new([&cram]).open().unwrap_err();
    assert!(
        matches!(error, ReadError::UnsupportedFormat { format: "CRAM", .. }),
        "{error:?}"
    );
    Ok(())
}

/// Returns a deterministic FASTA of `reads` reads of `len` pseudo-random bases each.
fn toy_fasta(reads: usize, len: usize) -> Vec<u8> {
    let mut state: u64 = 0x9e37_79b9_7f4a_7c15;
    let mut fasta = Vec::with_capacity(reads * (len + 16));
    for i in 0..reads {
        fasta.extend_from_slice(format!(">r{i}\n").as_bytes());
        for _ in 0..len {
            state = state
                .wrapping_mul(6_364_136_223_846_793_005)
                .wrapping_add(1_442_695_040_888_963_407);
            let index = usize::try_from(state >> 62).expect("two bits fit usize");
            fasta.push(b"ACGT"[index]);
        }
        fasta.push(b'\n');
    }
    fasta
}

#[test]
fn bgzf_with_many_blocks_reads_the_same_with_any_worker_count() -> TestResult {
    let dir = TempDir::new()?;
    // About 200 KB: at least four 64 KiB BGZF blocks.
    let fasta = toy_fasta(50, 4000);
    let compressed = bgzip(&fasta)?;
    let path = write_file(&dir, "big.fa.bgz", &compressed)?;
    let plain = write_file(&dir, "big.fa", &fasta)?;
    let expected = read_all(&ReadSource::new([&plain]))?;
    assert_eq!(expected.len(), 50);
    assert_eq!(expected[49].1.len(), 4000);
    for workers in [1, 3] {
        let options = InputOptions::new().with_bgzf_workers(NonZero::new(workers).ok_or("zero")?);
        let source = ReadSource::new([&path]).with_options(options);
        assert_eq!(read_all(&source)?, expected, "{workers} workers");
    }
    Ok(())
}

#[test]
fn reads_fastq_in_bgzf_and_bam_uncompressed_or_in_gzip() -> TestResult {
    let dir = TempDir::new()?;
    let fastq = write_file(&dir, "reads.fq.bgz", &bgzip(FASTQ)?)?;
    assert_eq!(read_all(&ReadSource::new([&fastq]))?, expected_fastq());
    assert_eq!(
        detect_format(&fastq)?,
        InputFormat::new(Compression::Bgzf, RecordFormat::Fastq)
    );

    let bam = dir.path().join("reads.bam");
    write_bam(
        &bam,
        &[
            unmapped("u1", b"ACGTT", Flags::empty()),
            // Unmapped but flagged reverse: still reverse-complemented back.
            unmapped("u2", b"AACCG", Flags::REVERSE_COMPLEMENTED),
        ],
    )?;
    let mut raw = Vec::new();
    std::io::Read::read_to_end(&mut bgzf::io::Reader::new(File::open(&bam)?), &mut raw)?;
    let uncompressed = write_file(&dir, "plain.bam", &raw)?;
    let gzipped = write_file(&dir, "gzipped.bam", &gzip(&raw)?)?;
    let expected = vec![
        ("u1".to_owned(), b"ACGTT".to_vec()),
        ("u2".to_owned(), b"CGGTT".to_vec()),
    ];
    for (path, compression) in [
        (&bam, Compression::Bgzf),
        (&uncompressed, Compression::None),
        (&gzipped, Compression::Gzip),
    ] {
        assert_eq!(
            detect_format(path)?,
            InputFormat::new(compression, RecordFormat::Bam)
        );
        assert_eq!(
            read_all(&ReadSource::new([path]))?,
            expected,
            "{compression}"
        );
    }
    Ok(())
}

#[test]
fn malformed_fasta_reports_file_and_record() -> TestResult {
    let dir = TempDir::new()?;
    let path = write_file(&dir, "bad.fa", b">r1\nACGT\n>\nACGT\n")?;
    let mut stream = ReadSource::new([&path]).open()?;
    assert!(stream.next().ok_or("missing first read")?.is_ok());
    let error = stream.next().ok_or("missing error")?.unwrap_err();
    assert!(
        matches!(&error, ReadError::Fasta { path: p, record: 2, .. } if p == &path),
        "{error:?}"
    );
    assert!(error.to_string().contains("record 2"), "{error}");
    assert!(stream.next().is_none());
    Ok(())
}

#[test]
fn next_batch_keeps_the_reads_before_an_error() -> TestResult {
    let dir = TempDir::new()?;
    let path = write_file(
        &dir,
        "bad.fq",
        b"@a\nAC\n+\nII\n@b\nGG\n+\nII\n@c\nTT\n+\nI\n",
    )?;
    let mut stream = ReadSource::new([&path]).open()?;
    let mut batch = ReadBatch::with_max_bases(1000);
    let error = stream.next_batch(&mut batch).unwrap_err();
    assert!(
        matches!(error, ReadError::Fastq { record: 3, .. }),
        "{error:?}"
    );
    let names: Vec<&[u8]> = batch.iter().map(|r| r.name.as_slice()).collect();
    assert_eq!(names, [&b"a"[..], b"b"]);
    assert_eq!(batch.bases(), 4);
    // The stream is finished: no more batches.
    assert!(!stream.next_batch(&mut batch)?);
    assert!(batch.is_empty());
    Ok(())
}

#[test]
fn current_path_and_format_follow_the_files() -> TestResult {
    let dir = TempDir::new()?;
    let fasta = write_file(&dir, "a.fa", b">r1\nAC\n")?;
    let fastq = write_file(&dir, "b.fq.gz", &gzip(b"@q1\nGG\n+\nII\n")?)?;
    let mut stream = ReadSource::new([&fasta, &fastq]).open()?;
    // The first file is open before any read.
    assert_eq!(stream.current_path(), Some(fasta.as_path()));
    let plain_fasta = InputFormat::new(Compression::None, RecordFormat::Fasta);
    assert_eq!(stream.current_format(), Some(plain_fasta));

    let mut record = ReadRecord::default();
    assert!(stream.read_record(&mut record)?);
    assert_eq!(stream.current_path(), Some(fasta.as_path()));
    assert!(stream.read_record(&mut record)?);
    assert_eq!(record.name, b"q1");
    assert_eq!(stream.current_path(), Some(fastq.as_path()));
    assert_eq!(
        stream.current_format(),
        Some(InputFormat::new(Compression::Gzip, RecordFormat::Fastq))
    );
    assert!(!stream.read_record(&mut record)?);
    assert_eq!(stream.current_path(), None);
    assert_eq!(stream.current_format(), None);
    Ok(())
}

#[test]
fn batch_memory_follows_its_contents() -> TestResult {
    let dir = TempDir::new()?;
    // Five 2-base reads fill the first batch, then long and short reads alternate.
    let mut fastq = Vec::new();
    for i in 0..5 {
        fastq.extend_from_slice(format!("@s{i}\nAC\n+\nII\n").as_bytes());
    }
    for (i, len) in [100_000, 100, 100_000, 100].into_iter().enumerate() {
        let seq = "A".repeat(len);
        fastq.extend_from_slice(format!("@l{i}\n{seq}\n+\n{seq}\n").as_bytes());
    }
    let path = write_file(&dir, "mixed.fq", &fastq)?;
    let mut stream = ReadSource::new([&path]).open()?;
    let mut batch = ReadBatch::with_max_bases(10);

    assert!(stream.next_batch(&mut batch)?);
    assert_eq!(batch.len(), 5);
    let mut sizes = Vec::new();
    while stream.next_batch(&mut batch)? {
        assert_eq!(batch.len(), 1);
        sizes.push((batch.bases(), batch.heap_bytes()));
        // Copies hold only the current read.
        let copy = batch.clone();
        assert_eq!(copy.records(), batch.records());
        assert!(copy.heap_bytes() <= batch.heap_bytes(), "{sizes:?}");
        assert_eq!(copy.max_bases(), 10);
    }
    assert_eq!(sizes.len(), 4);
    for &(bases, heap) in &sizes {
        // One read's buffers plus a record table of at most 1024 slots: far less than a
        // 100 kb buffer left over from an earlier read.
        let table = 1024 * size_of::<ReadRecord>();
        assert!(heap <= 4 * bases.max(256) + 1024 + table, "{sizes:?}");
    }

    // Clearing keeps the buffers, but a clone of the empty batch holds none of them.
    batch.clear();
    let empty = batch.clone();
    assert!(empty.is_empty());
    assert_eq!(empty.heap_bytes(), 0);

    // `clone_from` reuses buffers and copies only the current reads.
    let mut stream = ReadSource::new([&path]).open()?;
    let mut first = ReadBatch::with_max_bases(10);
    assert!(stream.next_batch(&mut first)?);
    let mut target = ReadBatch::new();
    target.clone_from(&first);
    assert_eq!(target.records(), first.records());
    assert_eq!(target.first_id(), Some(ReadId(0)));
    assert_eq!((target.bases(), target.max_bases()), (10, 10));
    Ok(())
}
