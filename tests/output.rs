//! Writer round trips: FASTA back through the reader, GFA line by line.

use std::{error::Error, fs};

use picoasm::io::{
    FastaWriter, GfaLink, GfaPath, GfaSegment, GfaTag, GfaWriter, OrientedSegment, Overlap,
    ReadSource, TagValue, WriteError,
};
use picoasm::seq::Strand;
use tempfile::TempDir;

type TestResult = Result<(), Box<dyn Error>>;

#[test]
fn fasta_round_trips_through_the_reader() -> TestResult {
    let dir = TempDir::new()?;
    let records: Vec<(String, Vec<u8>)> = vec![
        ("ctg1".to_owned(), b"ACGT".repeat(50)),
        ("ctg2".to_owned(), b"acgtNNRY".to_vec()),
        ("ctg3".to_owned(), Vec::new()),
        ("ctg4".to_owned(), b"G".repeat(161)),
    ];
    for width in [0, 1, 7, 80] {
        let path = dir.path().join(format!("contigs_{width}.fa"));
        let mut writer = FastaWriter::create(&path)?;
        if width != 80 {
            writer = FastaWriter::with_line_width(writer.finish()?, width);
        }
        for (name, seq) in &records {
            writer.write_record(format!("{name} len={}", seq.len()), seq)?;
        }
        writer.finish()?;

        let read: Vec<(String, Vec<u8>)> = ReadSource::new([&path])
            .open()?
            .map(|record| record.map(|r| (String::from_utf8_lossy(&r.name).into_owned(), r.seq)))
            .collect::<Result<_, _>>()?;
        assert_eq!(read, records, "width {width}");

        let text = fs::read_to_string(&path)?;
        if width > 0 {
            let longest = text
                .lines()
                .filter(|l| !l.starts_with('>'))
                .map(str::len)
                .max();
            assert_eq!(longest, Some(width), "width {width}");
        }
    }
    Ok(())
}

#[test]
fn gfa_writes_every_line_type() -> TestResult {
    let dir = TempDir::new()?;
    let path = dir.path().join("graph.gfa");
    let mut gfa = GfaWriter::create(&path)?;
    gfa.write_header(&[GfaTag::new(*b"PG", TagValue::String("picoasm test"))])?;
    let tags = [
        GfaTag::length(4),
        GfaTag::kmer_count(30),
        GfaTag::depth(7.5),
    ];
    gfa.write_segment(
        &GfaSegment::new("u1")
            .with_sequence(b"ACGT")
            .with_tags(&tags),
    )?;
    gfa.write_segment(&GfaSegment::new("u2").with_tags(&[GfaTag::length(1_000_000)]))?;
    gfa.write_segment(&GfaSegment::new("u3").with_sequence(b""))?;
    let u1 = OrientedSegment::new("u1", Strand::Forward);
    let u2 = OrientedSegment::new("u2", Strand::Reverse);
    gfa.write_link(&GfaLink::new(u1, u2, Overlap::Cigar("2M1I1M")))?;
    gfa.write_link(&GfaLink::new(u2, u1, Overlap::Unknown).with_tags(&[GfaTag::read_count(3)]))?;
    let steps = [u1, u2, OrientedSegment::new("u3", Strand::Forward)];
    gfa.write_path(&GfaPath::new("p1", &steps))?;
    let overlaps = [Overlap::Matches(2), Overlap::Cigar("0M")];
    gfa.write_path(&GfaPath::new("p2", &steps).with_overlaps(&overlaps))?;
    gfa.finish()?;

    let text = fs::read_to_string(&path)?;
    let expected = [
        "H\tVN:Z:1.0\tPG:Z:picoasm test",
        "S\tu1\tACGT\tLN:i:4\tKC:i:30\tdp:f:7.5",
        "S\tu2\t*\tLN:i:1000000",
        "S\tu3\t*",
        "L\tu1\t+\tu2\t-\t2M1I1M",
        "L\tu2\t-\tu1\t+\t*\tRC:i:3",
        "P\tp1\tu1+,u2-,u3+\t*",
        "P\tp2\tu1+,u2-,u3+\t2M,0M",
    ];
    assert_eq!(text.lines().collect::<Vec<_>>(), expected);
    assert!(text.ends_with('\n'));
    Ok(())
}

#[test]
fn gfa_rejects_invalid_lines_without_writing() {
    let mut gfa = GfaWriter::new(Vec::new());
    let ok = OrientedSegment::new("a", Strand::Forward);
    let bad = OrientedSegment::new("a b", Strand::Forward);
    let steps = [ok, ok];

    let results = [
        gfa.write_segment(&GfaSegment::new("*")),
        gfa.write_segment(&GfaSegment::new("a").with_sequence(b"AC\tGT")),
        gfa.write_segment(&GfaSegment::new("a").with_tags(&[GfaTag::depth(f64::INFINITY)])),
        gfa.write_link(&GfaLink::new(ok, bad, Overlap::Unknown)),
        gfa.write_link(&GfaLink::new(ok, ok, Overlap::Cigar("5"))),
        gfa.write_path(&GfaPath::new("p", &[])),
        gfa.write_path(&GfaPath::new("p", &steps).with_overlaps(&[Overlap::Matches(1); 2])),
        gfa.write_path(&GfaPath::new("p", &steps).with_overlaps(&[Overlap::Unknown])),
        gfa.write_path(&GfaPath::new("p", &[ok, bad])),
    ];
    for (i, result) in results.iter().enumerate() {
        assert!(result.is_err(), "case {i} should fail");
    }
    assert!(matches!(results[0], Err(WriteError::InvalidName { .. })));
    assert!(matches!(results[5], Err(WriteError::EmptyPath { .. })));
    assert!(matches!(
        results[6],
        Err(WriteError::PathOverlapCount {
            steps: 2,
            overlaps: 2,
            ..
        })
    ));
    assert!(
        gfa.get_ref().is_empty(),
        "nothing is written for invalid lines"
    );
}

/// A writer that always fails, standing in for a full disk.
struct FullDisk;

impl std::io::Write for FullDisk {
    fn write(&mut self, _: &[u8]) -> std::io::Result<usize> {
        Err(std::io::Error::other("disk full"))
    }

    fn flush(&mut self) -> std::io::Result<()> {
        Ok(())
    }
}

#[test]
fn write_errors_name_the_output_file() -> TestResult {
    let error = FastaWriter::new(FullDisk)
        .write_record("ctg", "ACGT")
        .unwrap_err();
    assert!(
        matches!(error, WriteError::Io { path: None, .. }),
        "{error:?}"
    );
    let error = GfaWriter::new(FullDisk).write_header(&[]).unwrap_err();
    assert!(
        matches!(error, WriteError::Io { path: None, .. }),
        "{error:?}"
    );

    // On Linux, /dev/full fails every write with "no space left on device", like a full disk.
    let full = std::path::Path::new("/dev/full");
    if cfg!(target_os = "linux") && full.exists() {
        let mut fasta = FastaWriter::create(full)?;
        fasta.write_record("ctg", "ACGT")?; // Buffered, so this succeeds.
        let error = fasta.finish().unwrap_err();
        assert!(
            matches!(&error, WriteError::Io { path: Some(p), .. } if p == full),
            "{error:?}"
        );
        assert!(error.to_string().contains("/dev/full"), "{error}");

        let mut gfa = GfaWriter::create(full)?;
        gfa.write_header(&[])?;
        let error = gfa.finish().unwrap_err();
        assert!(
            matches!(&error, WriteError::Io { path: Some(p), .. } if p == full),
            "{error:?}"
        );
    }
    Ok(())
}
