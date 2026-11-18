# picoASM

[![CI](https://github.com/DylanCosto/picoASM/actions/workflows/ci.yml/badge.svg)](https://github.com/DylanCosto/picoASM/actions/workflows/ci.yml)
[![Release](https://img.shields.io/github/v/release/DylanCosto/picoASM)](https://github.com/DylanCosto/picoASM/releases/latest)

picoASM is a de novo assembler for bacterial isolate genomes sequenced with Oxford Nanopore
(R10.4.1) or PacBio HiFi reads. It aims to give every replicon (the chromosome and each plasmid,
small high-copy plasmids included) as one circular contig, without a separate polishing step,
in a few CPU-minutes.

On the 30 R10.4.1 read sets of the Autocycler benchmark (Wick, Howden & Stinear 2025),
picoASM completed every replicon in all 30 assemblies. Its accuracy (median 2 errors per genome,
86 in all) was close to that of Autocycler followed by Medaka (2 and 78; 25 of 30 complete).
It took a median of 56 s wall and 6.1 CPU-minutes per genome at 8 threads, with at most 1.5 GB
of memory. On a validation set of 20 isolates from genera not used in development (10 ONT,
10 HiFi), it was ahead of Flye (+ Medaka for ONT) on 6, level on 8 and behind on 6, using 4 to 5
times less CPU time.

- [Installation](#installation)
- [Usage](#usage)
- [Outputs](#outputs)
- [Benchmark](#benchmark)
- [Results](#results)
- [Limitations](#limitations)
- [Reproducing the benchmark](#reproducing-the-benchmark)
- [Citation](#citation)

## Installation

**Prebuilt binaries.** Each [release](https://github.com/DylanCosto/picoASM/releases) has binaries
for Linux (x86-64) and macOS (Apple silicon):

```bash
curl -LO https://github.com/DylanCosto/picoASM/releases/download/v0.5.0/picoasm-v0.5.0-x86_64-unknown-linux-gnu.tar.gz
tar -xzf picoasm-v0.5.0-x86_64-unknown-linux-gnu.tar.gz
./picoasm-v0.5.0-x86_64-unknown-linux-gnu/picoasm --version
```

**From source.** picoASM is written in Rust and needs a Rust toolchain, 1.91 or later
(https://rustup.rs):

```bash
git clone https://github.com/DylanCosto/picoASM.git
cd picoASM
cargo build --release
```

The binary is `target/release/picoasm`. It has no run-time dependencies.

## Usage

```bash
picoasm assemble --reads ont -t 8 -o contigs.fa reads.fastq.gz
```

The input can be FASTQ or FASTA, gzipped or not, or unaligned BAM. Several files may be given.
`-t` sets the number of threads (default: all cores). `--reads` sets the read type:

| `--reads` | Reads |
|---|---|
| `ont` (default) | Oxford Nanopore R10.4.1, Dorado sup or hac. Noisier reads (R9.4.1, or fast basecalling: 5-18% divergence) are detected from their overlaps and corrected first. |
| `ont-noisy` | Nanopore reads, always corrected first. |
| `ont-clean` | Nanopore reads, never corrected. |
| `hifi` | PacBio HiFi (CCS). |

`-o` writes the contigs; [Outputs](#outputs) describes them and the other files picoASM can
write.

**Illumina polishing.** `--short-reads R1 --short-reads R2` polishes the contigs at the end with
short reads of the same isolate. `picoasm polish-short ASSEMBLY R1 R2` does the same for an
existing assembly.

**Optional GPU polishing.** On heavily methylated genomes, `dorado polish --bacteria` after
picoASM corrects some remaining errors at modified motifs. On genomes without strong
methylation it changes little. It needs reads basecalled by Dorado:

```bash
dorado aligner -t 8 contigs.fa reads.fastq.gz | samtools sort -o aligned.bam -
samtools index aligned.bam
dorado polish --bacteria -t 8 aligned.bam contigs.fa > polished.fa
```

## Outputs

```bash
picoasm assemble --reads ont -t 8 -o contigs.fa --suspect suspect.bed --gfa contigs.gfa \
    --shallow shallow.fa reads.fastq.gz
```

| Option | Format | Contents |
|---|---|---|
| `-o` | FASTA | The contigs, longest first, named `contig_1`, `contig_2` and on. |
| `--suspect` | BED | Where the reads show a contig may be wrong. |
| `--gfa` | GFA 1 | The contigs as a graph, for Bandage and other graph tools. |
| `--shallow` | FASTA | Contigs left out for low depth. |

Only `-o` is required. The contigs are named the same way in every file.

### Contigs (`-o`)

Each header gives the contig's length, whether it is a closed circular replicon
(`circular=true`), any hairpin telomeres (`telomeres=hairpin`, or `hairpin=start` or
`hairpin=end` for one), its read depth, and a check of the contig against the reads:

```
>contig_1 length=5153435 circular=true depth=40.8 check=stops:5
```

- `check=ok`: the reads show nothing wrong;
- or what is wrong: `linear` (not closed into a whole replicon), `shallow:N` (N 1 kb bins
  under a quarter of the median depth), `stops:N` (N bins where 20% or more of the reads stop
  aligning).

Any problem is also logged as a warning. `--suspect` gives where each one is.

### Suspect regions (`--suspect`)

A BED file of the places the reads show a contig may be wrong, to open beside the assembly in
IGV or intersect with an annotation (`bedtools intersect`):

```
#contig    start    end      reason
contig_1   5094000  5095000  stops:reads=14,depth=27.2
contig_1   5096000  5097000  stops:reads=18,depth=28.4
contig_2   0        403967   linear
contig_2   222130   224130   shallow:depth=24.3,median=151.8
contig_2   395203   403967   shallow:depth=5.4,median=151.8
```

| Reason | Meaning | Usual causes |
|---|---|---|
| `shallow` | Under a quarter of the contig's median depth; gives the lowest depth and the median. | A gap bridged on thin evidence, another replicon run in, part of the population lacking the region. |
| `stops` | 20% or more of the reads stop aligning; gives the bin where most stop and its depth. | A join that is not in the genome, a locus present in two versions (an invertible segment, a subclone), a collapsed repeat. |
| `linear` | The whole contig: it is not closed into a whole replicon. | Repeats longer than the reads, or too few reads across the circle's join. |

Regions are 1 kb bins (the last takes the rest of the contig, 1-2 kb), merged where adjacent,
in the coordinates of the final contigs. A contig with `check=ok` has no line, and an assembly
the reads fully support gives a file with only its header. In the example, the first two lines
(*E. coli* H10407, HiFi) bracket a 1.9 kb segment the reference has the other way round: an
invertible segment the reads hold both ways. The rest are from *A. tenax* TSD-302 (ONT).

### Graph (`--gfa`)

The contigs as GFA 1 segments with their sequences, lengths (`LN:i`) and depths (`dp:f`). A
circular contig links its end to its start (`L contig_1 + contig_1 + 0M`); a hairpin end links
the contig to its own reverse complement. Contigs are not linked to each other, since picoASM
resolves repeats itself.

### Left-out contigs (`--shallow`)

Contigs below 10% of the chromosome's depth (25% if linear) are left out of the assembly:
typically reads of another sample through a barcode, or chimeras. `--shallow FILE` writes them,
named `shallow_1` and on, with the same headers as the contigs.

## Benchmark

`scripts/bact_bench.py` runs an assembler on each dataset of a list (`scripts/bact_datasets*.tsv`)
under `/usr/bin/time -v`, and scores the contigs against the dataset's reference with
`scripts/bact_score.py`. Per reference replicon:

- **Complete**: one contig contains the replicon whole, within 1% of its length.
- **Errors**:
  - Of a complete replicon: every mismatch, inserted base and deleted base in an end-to-end
    alignment of the rotated contig to the replicon, so a faulty circle junction counts.
  - Of an incomplete replicon: the errors within its contigs' alignments. Uncovered bases are
    counted separately as **missing**.
- **Misassemblies**: consecutive alignments of 1 kb or more from one contig that fall on
  different replicons or strands, or more than 1 kb out of place.
- **Extra contigs**: contigs that align nowhere, or duplicate a replicon held by another contig.
- **Perfect**: every replicon complete, with no error, extra contig or misassembly.

On the Autocycler paper's own assemblies these error counts agree with the paper's
`assess_assembly.py` to within a few errors. Errors shared by every tool at the same position
are usually differences between the sequenced culture and the reference.

**Datasets.** Development data and test data are kept apart. Lists marked *development* were
used to design and tune picoASM. Every other list was run after the method was fixed, except
where noted. Per-dataset results are in `benchmarks/results.tsv`.

| List | Datasets | Reads | Use |
|---|---|---|---|
| `subsets` | 5 genomes x 6 read subsets (~50x) | R10.4.1, Dorado sup v5 (Wick, Howden & Stinear 2025) | development |
| `unseen`, `heldout` | 11 isolates (~60x) | R10.4.1, sup/hac v4.3 (Hall et al. 2024) | test |
| `wide`, `wide2`, `extra` | 24 isolates (~60x) | R10.4.1 rapid and ligation, R9.4.1; public runs | test (failures in `wide2` led to general fixes) |
| `fresh` | 12 ONT, 4 HiFi isolates | R10.4.1 (~60x); PacBio 96-plex HiFi (~35x) | test (4 later used for fixes) |
| `more` | 24 ONT, 20 HiFi isolates | recent ONT runs; PacBio 96-plex HiFi | test (3 used in development) |
| `final` | 10 ONT, 10 HiFi isolates (~100x), each of a genus used nowhere else | R10.4.1; HiFi | validation, run once |
| `lowdepth`, `verylow`, `highdepth` | depth series of the `subsets` genomes | as `subsets` | development |

**Compared tools.**
- Flye 2.9.6 (`--nano-hq` or `--pacbio-hifi`).
- Flye followed by Medaka 2.2.2 with the `r1041_e82_400bps_bacterial_methylation` model. This
  is the core of ONT's EPI2ME `wf-bacterial-genomes`.
- hifiasm 0.25.0 (primary contigs).
- On the Autocycler subsets: the paper's own assemblies by every tool it evaluated (figshare
  10.26188/28910099), scored with the same script.

**Hardware.** Every run used 8 threads, on a shared 24-core Intel Xeon Gold 5418Y, with up to
three runs at a time. The exception is Flye and Flye + Medaka on the `subsets`, which used 16
threads. Wall times are therefore indicative; CPU times are comparable.

## Results

Cells give complete replicons / reference replicons, then errors. Where there are any, they
then give missing bases and, after a slash, misassemblies. All picoASM results are from version
0.5.0 with default settings.

### Autocycler benchmark (30 read subsets; development data)

| Tool | All replicons complete | Perfect | Median errors | Total errors | Missing bases | Extra contigs | Misassemblies |
|---|---|---|---|---|---|---|---|
| **picoASM** | **30/30** | **6** | **2** | 86 | **0** | **0** | **0** |
| Autocycler + Medaka | 25/30 | 5 | **2** | **78** | 12,475 | **0** | **0** |
| Autocycler (manual curation) | **30/30** | 2 | 4 | 135 | **0** | **0** | **0** |
| Autocycler | 25/30 | 1 | 4 | 133 | 12,475 | **0** | **0** |
| Flye + Medaka | 12/30 | 0 | 3 | 4,358 | 29,993 | 8 | **0** |
| Flye | 12/30 | 0 | 4.5 | 4,454 | 29,993 | 8 | 1 |
| Hybracter | 28/30 | 1 | 11 | 7,409 | 2,398 | **0** | **0** |
| Canu | 17/30 | 0 | 10 | 105,007 | 36,337 | 9 | 1 |

The other tools in the paper (miniasm, Dragonflye, LJA, MAECI, hifiasm, NECAT, metaMDBG, Raven,
Myloasm, Redbean, NextDenovo) had median errors of 32 or more. They are in
`benchmarks/autocycler_paper.tsv`. Median errors per genome:

| Genome | picoASM | Autocycler + Medaka | Autocycler | Flye + Medaka | Flye |
|---|---|---|---|---|---|
| *Shigella flexneri* | 6 | 4 | 8 | 57 | 60.5 |
| *Klebsiella pneumoniae* | 4 | 6.5 | 8.5 | 9 | 14.5 |
| *Enterobacter hormaechei* | 0.5 | 0 | 1 | 1 | 2.5 |
| *Listeria innocua* | 0.5 | 0 | 1 | 1 | 2 |
| *Providencia rettgeri* | 2 | 2 | 3.5 | 3 | 4.5 |

Three of picoASM's error sites appear to be errors in the references; every tool makes them.
- *P. rettgeri* 1332350: 45 of 47 reads, and the Illumina-polished Hybracter assembly, have
  the other base.
- *P. rettgeri* 3535633: a G6 that all reads on both strands give as G5.
- *S. flexneri* 1178130: most reads give 11 G where the reference has 8.

Together these account for 34 of the 86 errors.

**Resources.** Median 56 s wall time and 6.1 CPU-minutes per genome. Median peak memory was
1.3 GB, at most 1.5 GB. For comparison, Autocycler run end to end on one *S. flexneri* subset
took 7 h wall time and 84 CPU-hours at 16 threads.

### Hall et al. 2024 isolates (test)

| Isolate | picoASM | Flye | Flye + Medaka |
|---|---|---|---|
| *L. ivanovii* | 1/1, 10 | 1/1, 21 | 1/1, 10 |
| *L. welshimeri* | 1/1, 6 | 1/1, 12 | 1/1, 7, 0 kb / 1 |
| *C. lari* | 1/1, 12 | 1/1, 62 | 1/1, 19 |
| *S. dysgalactiae* | 1/1, 16 | 1/1, 27 | 1/1, 9 |
| *S. aureus* | 3/3, 1 | 2/3, 2 | 2/3, 2 |
| *S. pyogenes* | 1/1, 3 | 1/1, 7 | 1/1, 4, 0 kb / 1 |
| *C. jejuni* | 1/1, 19 | 1/1, 74 | 1/1, 17 |
| *L. monocytogenes* | 1/1, 7 | 1/1, 13 | 1/1, 1 |
| *K. variicola* | 1/1, 2 | 1/1, 4 | 1/1, 2 |
| *E. coli* (reads N50 1.2 kb) | 3/5, 5,815, 226.3 kb / 6 | 1/5, 6,057, 805.3 kb / 14 | 1/5, 5,897, 804.7 kb / 15 |
| *S. enterica* (reads N50 1.2 kb) | 1/2, 290, 62.3 kb / 5 | 1/2, 617, 416.4 kb / 17 | 1/2, 484, 422.1 kb / 17 |
| *M. tuberculosis* mc2 6030 (development) | 1/1, 0 | 1/1, 12 | 1/1, 1 |
| *K. pneumoniae* KPC2 (development) | 4/4, 2 | 2/4, 31 | 2/4, 2 |

On the nine test isolates with reads of normal length, picoASM made 76 errors, against 222
for Flye and 71 for Flye + Medaka.

### Final validation set (20 isolates, run once)

Each isolate is of a genus that appears in no other list, and is compared with its NCBI complete
genome. Median wall time:
- ONT: picoASM 66 s, Flye + Medaka 562 s.
- HiFi: picoASM 53 s, Flye 422 s, hifiasm 332 s.

Median CPU time:
- ONT: picoASM 7.5 min, Flye + Medaka 29.5 min.
- HiFi: picoASM 6.1 min, Flye 29.7 min, hifiasm 26.0 min.

| ONT isolate | picoASM | Flye + Medaka |
|---|---|---|
| *Seonamhaeicola* sp. ABD4_340 | 1/1, 3 | 1/1, 3 |
| *Ligilactobacillus salivarius* 11254D021BM | 1/2, 28 | 1/2, 1,586 |
| *Lactococcus lactis* S11-599 | 2/2, 17, 0 kb / 1 | 1/2, 9,595 |
| *Psychrobacillus psychrotolerans* BC2 | 1/1, 6 | 1/1, 1 |
| *Yersinia pestis* 80 | 4/4, 2,330 | 3/4, 3,585 |
| *Citrobacter portucalensis* MEI2O2-23 (ligation) | 3/4, 6, 3.4 kb | 4/4, 2 |
| *Pantoea septica* TB5-P01A | 2/8, 2,815, 41.1 kb / 1 | 2/8, 2,810, 41.1 kb / 1 |
| *Pseudarthrobacter* sp. WHRI 8279 | 3/4, 36, 0 kb / 1 | 4/4, 4 |
| *Cereibacter sphaeroides* DSM 158 | 6/6, 1,635 | 6/6, 45 |
| *Arsenicicoccus dermatophilus* KM18/12 | 3/4, 151, 53.4 kb | 3/4, 5,436 |

| HiFi isolate | picoASM | Flye | hifiasm |
|---|---|---|---|
| *Anaerosalibacter bizertensis* JCM 17239 | 1/1, 1 | 1/1, 0 | 1/1, 8,280 |
| *Priestia megaterium* IMGN6 | 12/12, 0 | 9/12, 0, 24.6 kb | 1/12, 8, 43.6 kb / 2 |
| *Furfurilactobacillus cerevisiae* LTH 5742 | 2/3, 42,122, 131.3 kb | 2/3, 42,107, 108.5 kb | 1/3, 17,487, 42.2 kb / 2 |
| *Aneurinibacillus thermoaerophilus* CCM 8960 | 2/2, 103 | 2/2, 1,643 | 1/2, 103 |
| *Psychrobacter* sp. KFRI-CH2-11 | 4/5, 21 | 1/5, 27 | 1/5, 42 |
| *Kozakia baliensis* LTH 7215 | 4/6, 6,099, 8.4 kb / 1 | 4/6, 6, 6.1 kb | 1/6, 6,098, 8.8 kb / 4 |
| *Trinickia caryophylli* ATCC 25418 | 2/2, 1,904 | 2/2, 2,045 | 2/2, 1,349 |
| *Rhodospirillum rubrum* DSM 467 | 2/2, 3,759 | 2/2, 3,760 | 1/2, 25,354 |
| *Caldimonas thermodepolymerans* DSM 15264 | 1/1, 12 | 1/1, 12 | 0/1, 132,469 |
| *Micromonospora soli* SL3-70 | 1/1, 3,228 | 1/1, 1 | 1/1, 4,323 |

Totals over the 20 isolates:
- ONT: picoASM completed 26 of 36 replicons with 7,027 errors; Flye + Medaka 26 of 36 with
  23,067.
- HiFi: picoASM completed 31 of 35 replicons; Flye 25 and hifiasm 10.

picoASM was ahead on 6 isolates, level on 8 and behind on 6. The six where it was behind:
- A 3.4 kb plasmid from a ligation library is missed. Only 5 reads contain it whole, and small
  circles need whole-read support at a quarter of the chromosome's depth.
- A 212.8 kb plasmid is left linear, its ends overlapping by about 27 kb.
- A 1.6 kb insertion that 5 of 76 reads carry is kept. The other reads align straight across
  it, so no reads stop at the locus.
- A 6.1 kb plasmid is laid out inside the chromosome, at a repeat the two share.
- At a mixed site, picoASM keeps an insertion that 24 reads have and 29 lack.
- Five more errors than Flye + Medaka on *P. psychrotolerans*.

Every assembler misses the same sequence of *P. septica* and *F. cerevisiae*. Errors shared with
Flye at the same positions (*Y. pestis*, *T. caryophylli*, *R. rubrum*) are mostly differences
from the reference.

### Wider validation (test)

| Isolate | picoASM | Flye + Medaka |
|---|---|---|
| *P. aeruginosa* PAO1 (R10, rapid) | 1/1, 142 | 1/1, 138 |
| *P. aeruginosa* PAO1 (R10, ligation) | 1/1, 167 | 1/1, 132 |
| *P. aeruginosa* PAO1 (R9.4.1) | 1/1, 334 | 1/1, 513 |
| *E. coli* CFT073 | 1/1, 53 | 1/1, 52 |
| *K. pneumoniae* MGH 78578 | 6/6, 308 | 2/6, 437, 60.5 kb / 2 |
| *S. aureus* MRSA252 | 1/1, 29 | 1/1, 29 |
| *V. parahaemolyticus* ATCC 17802 (reads N50 1.2 kb) | 2/2, 1,613 | 1/2, 1,659, 2.9 kb |
| *B. arboris* LMG 14939 | 5/5, 0 | 5/5, 5 |
| *B. arboris* R-132 | 3/3, 3 | 3/3, 22 |
| *B. arboris* LMG 24066 | 4/5, 8, 40.7 kb | 4/5, 17,822, 40.7 kb |
| *Burkholderia* sp. LMG 33096 | 4/4, 1 | 4/4, 2 |
| *Burkholderia* sp. LMG 33096 (second run) | 3/3, 0 | 3/3, 1 |
| *Burkholderia* sp. LMG 34232 | 3/3, 17 | 3/3, 8 |
| *B. mallei* ATCC 23344 | 2/2, 352 | 2/2, 385 |
| *E. faecalis* V583 | 4/4, 1,216 | 4/4, 1,143 |
| *E. coli* ATCC 25922 (reads 8.3% divergent) | 5/5, 300 | 3/5, 5,780, 2.2 kb |
| *S. enterica* ATCC 10708 | 2/2, 10 | 2/2, 13 |
| *A. baumannii* ATCC 19606 (reads 10.6% divergent) | 3/3, 20,242 | 1/3, 21,235, 0.3 kb / 2 |
| *S. turgidiscabies* (reads 7.3% divergent, 14x) | 1/4, 12,526, 231.1 kb / 3 | 2/4, 18,080, 18.1 kb / 2 |
| *P. syringae* pv. tomato DC3000 | 2/3, 3,939, 0 kb / 2 | 0/3, 3,923, 48.1 kb / 3 |
| *P. savastanoi* pv. phaseolicola 1448A (structure only) | 2/3, n/a, 51.7 kb / 3 | 2/3, n/a, 51.7 kb / 3 |

The two *Pseudomonas* lab strains are compared with Sanger-finished references, from which the
cultures have drifted. For *P. savastanoi* every tool shows the same ~70,000 differences, so
only structure is scored. Three *Streptomyces* isolates with draft references are left out of
the table: there the missing bases and misassemblies measure the reference as much as the
assembly. They are in `benchmarks/results.tsv`.

### Fresh isolates (test)

| ONT isolate | picoASM | Flye + Medaka |
|---|---|---|
| *B. subtilis* DSM 10 | 2/2, 4,438 | 2/2, 4,438 |
| *C. difficile* NCTC 12727 | 0/1, 222, 18.4 kb | 0/1, 635 |
| *S. capitis* DSM 20326 | 4/4, 24 | 4/4, 9 |
| *S. xylosus* CSKR35 | 6/6, 999 | 4/6, 1,154, 8.7 kb |
| *P. siamensis* LMG 24279 | 3/3, 933 | 3/3, 72 |
| *C. beijerinckii* DSM 105335 | 1/1, 63 | 1/1, 62 |
| *L. plantarum* ATCC 14917 | 2/2, 1,239 | 1/2, 1,603 |
| *Xanthomonas* sp. NCPPB 1131 | 1/1, 15 | 1/1, 2,695 |
| *C. necator* ATCC 25207 | 2/2, 101, 0 kb / 3 | 2/2, 62, 0 kb / 3 |
| *S. thermophilus* DSM 20745 | 2/2, 472 | 2/2, 473 |
| *M. ruber* DSM 1279 | 1/1, 23 | 1/1, 19 |
| *S. oneidensis* MR-1 | 1/2, 6,797, 3.6 kb / 2 | 1/2, 3,848, 8.6 kb / 1 |

| HiFi isolate | picoASM | Flye | hifiasm |
|---|---|---|---|
| *B. multivorans* ATCC 17616 | 4/4, 58 | 4/4, 59 | 4/4, 9,843 |
| *R. palustris* CGA009 | 2/2, 2 | 1/2, 2, 8.4 kb | 1/2, 2 |
| *N. meningitidis* FAM18 | 1/1, 13 | 1/1, 14 | 1/1, 13 |
| *X. oryzae* M10 | 1/1, 54 | 0/1, 71 | 1/1, 54 |

*C. difficile*, *P. siamensis*, *N. meningitidis* and *X. oryzae* M10 were later used to find
fixes, so their rows are not a clean test.

### 44 more isolates (test)

- **ONT** (23 isolates also assembled by Flye + Medaka): every replicon was complete in 16
  isolates with picoASM, against 15 with Flye + Medaka. picoASM completed 35 of 55 replicons
  and made 26,238 errors; Flye + Medaka completed 32 and made 48,976.
- **HiFi** (20 isolates): every replicon was complete in 17 isolates with picoASM, against 16
  with Flye and 12 with hifiasm. Of 39 replicons, picoASM completed 33, Flye 31 and hifiasm 20.

The HiFi references are PacBio's own assemblies from the same reads, so they are not
independent.

### Read depth

| Read set | picoASM | Flye | Flye + Medaka |
|---|---|---|---|
| *S. flexneri* ~25x | 8/8, 30 | 2/8, 3,073, 1.2 kb | 3/8, 2,812, 1.2 kb / 1 |
| *K. pneumoniae* ~25x | 5/5, 6 | 3/5, 27, 1.2 kb | 3/5, 11, 1.2 kb |
| *E. hormaechei* ~25x | 4/4, 6 | 4/4, 4 | 4/4, 2 |
| *L. innocua* ~25x | 2/2, 5 | 2/2, 5 | 2/2, 4 |
| *P. rettgeri* ~25x | 2/2, 11 | 2/2, 10 | 2/2, 6 |
| *K. pneumoniae* 20x | 5/5, 26 | | 3/5, 32, 1.2 kb |
| *S. flexneri* 20x | 7/8, 108 | | 1/8, 229, 3.6 kb / 1 |
| *E. hormaechei* 15x | 3/4, 184, 13.8 kb | | 2/4, 134, 6.7 kb |
| *L. innocua* 12x | 0/2, 409, 5.9 kb | | 1/2, 406, 2.1 kb / 1 |
| *P. rettgeri* 12x | 1/2, 450, 10.1 kb | | 1/2, 324, 2.4 kb |
| *K. pneumoniae* ~200x | 5/5, 3 | | |
| *S. flexneri* ~200x | 8/8, 4 | | |

## Limitations

- **Small plasmids from ligation libraries.** These are under-represented among whole reads,
  and can fall below the support picoASM requires for a small circle (see *C. portucalensis*).
- **Minority versions that reads align across.** A locus is repaired only where reads stop
  aligning. An insertion or deletion that most reads span within one alignment is not
  detected (see *C. sphaeroides*).
- **Near-even mixed populations.** A locus takes the reads' version only when three in five
  reads carry it.
- **Repeats longer than the reads.** Without a read across the repeat, the layout's choice of
  flanks stands. A *C. difficile* chromosome with a 34.5 kb repeat in two copies comes out in
  pieces.
- **Low depth (12-15x).** Chromosomes can break where a single read crosses a dip in
  coverage. Flye leaves fewer bases out at these depths.
- **Long homopolymers (11 bases or more).** These are usually called short. Reads on both
  strands under-call them, so the genome's own reads cannot show by how much.
- **Noisy reads** (R9.4.1 fast, or reads more than 7.5% divergent). Correction adds 10-30
  minutes at 8 threads. Medaka still polishes some of these genomes better.
- **Short reads** (N50 about 1.2 kb). No tool tested completes these chromosomes.

## Reproducing the benchmark

Reads, references and assemblies are not in this repository. The dataset lists in
`scripts/bact_datasets*.tsv` give the expected paths under `data/bacteria/`. Sources:
- The Autocycler paper's reads (figshare 10.26188/28910213), and its references and
  assemblies (10.26188/28910099).
- Hall et al. 2024 (BioProject PRJNA1087001).
- PacBio's 2021-11-Microbial-96plex dataset.
- The ENA/SRA runs and NCBI assemblies of the individual isolates.

```bash
# One tool on one list, 3 runs at a time with 8 threads each:
BACT_DATASETS=scripts/bact_datasets_final.tsv python3 scripts/bact_bench.py --tool picoasm --label picoasm --threads 8 --jobs 3
# The results tables of this README, from benchmarks/results.tsv:
python3 scripts/bact_tables.py
```

The compared tools were installed with micromamba (Flye, Medaka, hifiasm and the Autocycler
assemblers, following the Autocycler environment file); `scripts/bact_bench.py` describes how
each is run. `scripts/bact_core.sh` runs a 17-dataset regression set and reports any assembly
that changed.

## Citation

If you use picoASM, please cite this repository (version 0.5.0). Please also cite the data used
here:
- R. R. Wick, B. P. Howden, T. P. Stinear (2025). Autocycler: long-read consensus assembly for
  bacterial genomes. *Bioinformatics*, btaf474.
- M. B. Hall et al. (2024). Benchmarking reveals superiority of deep learning variant callers on
  bacterial nanopore sequence data. *eLife* 13, RP98300.

## License

MIT; see `LICENSE`.
