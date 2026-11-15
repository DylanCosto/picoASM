"""Prints the README's results tables from benchmarks/results.tsv and benchmarks/autocycler_paper.tsv.

Usage: python3 scripts/bact_tables.py

benchmarks/results.tsv holds one row per dataset, read type and tool (picoasm, Flye, Flye + Medaka,
hifiasm, ...), as scripts/bact_bench.py writes them, with the dataset's list in `set`;
benchmarks/autocycler_paper.tsv holds the Autocycler paper's own assemblies of its 30 read subsets,
scored the same way (scripts/bact_score_paper.py).
"""
import collections
import csv
import pathlib
import statistics

ROOT = pathlib.Path(__file__).resolve().parent.parent
BENCH = ROOT / 'benchmarks'

NAMES = {
    # Hall et al. 2024
    'livanovii': '*L. ivanovii*', 'lwelshimeri': '*L. welshimeri*', 'clari': '*C. lari*',
    'sdysgalactiae': '*S. dysgalactiae*', 'mtb': '*M. tuberculosis* mc2 6030', 'kpc2': '*K. pneumoniae* KPC2',
    'saureus': '*S. aureus*', 'spyogenes': '*S. pyogenes*', 'cjejuni': '*C. jejuni*',
    'lmono': '*L. monocytogenes*', 'kvari': '*K. variicola*', 'ecoli': '*E. coli* (reads N50 1.2 kb)',
    'senterica': '*S. enterica* (reads N50 1.2 kb)',
    # wider validation
    'paeru': '*P. aeruginosa* PAO1 (R10, rapid)', 'paeru_lig': '*P. aeruginosa* PAO1 (R10, ligation)',
    'paeru_r9': '*P. aeruginosa* PAO1 (R9.4.1)', 'ecoli_cft073': '*E. coli* CFT073',
    'kpneu_mgh': '*K. pneumoniae* MGH 78578', 'mrsa252': '*S. aureus* MRSA252',
    'vpara': '*V. parahaemolyticus* ATCC 17802 (reads N50 1.2 kb)', 'barboris': '*B. arboris* LMG 14939',
    'barb_r132': '*B. arboris* R-132', 'barb_24066': '*B. arboris* LMG 24066',
    'bsp33096': '*Burkholderia* sp. LMG 33096', 'bsp_r72673': '*Burkholderia* sp. LMG 33096 (second run)',
    'bsp34232': '*Burkholderia* sp. LMG 34232', 'bmallei': '*B. mallei* ATCC 23344',
    'efaecalis': '*E. faecalis* V583', 'ecoli_25922': '*E. coli* ATCC 25922 (reads 8.3% divergent)',
    'senterica_b': '*S. enterica* ATCC 10708', 'abaum': '*A. baumannii* ATCC 19606 (reads 10.6% divergent)',
    'sturgid': '*S. turgidiscabies* (reads 7.3% divergent, 14x)',
    'seuropae': '*S. europaeiscabiei* (draft reference)', 'strep_stell': '*S. stelliscabiei* (draft reference)',
    'strep_puni': '*S. puniciscabiei* (draft reference)',
    'psph': '*P. savastanoi* pv. phaseolicola 1448A (structure only)', 'pst': '*P. syringae* pv. tomato DC3000',
    # fresh isolates
    'bsub_dsm10': '*B. subtilis* DSM 10', 'cdiff_nctc12727': '*C. difficile* NCTC 12727',
    'scap_dsm20326': '*S. capitis* DSM 20326', 'sxyl_cskr35': '*S. xylosus* CSKR35',
    'psia_lmg24279': '*P. siamensis* LMG 24279', 'cbei_dsm105335': '*C. beijerinckii* DSM 105335',
    'lpla_atcc14917': '*L. plantarum* ATCC 14917', 'xoo_whri5233': '*X. oryzae* WHRI 5233',
    'xsp_ncppb1131': '*Xanthomonas* sp. NCPPB 1131', 'cnec_atcc25207': '*C. necator* ATCC 25207',
    'sthe_dsm20745': '*S. thermophilus* DSM 20745', 'mrub_dsm1279': '*M. ruber* DSM 1279',
    'sone_mr1': '*S. oneidensis* MR-1', 'bmul_atcc17616_hifi': '*B. multivorans* ATCC 17616',
    'rpal_cga009_hifi': '*R. palustris* CGA009', 'nmen_fam18_hifi': '*N. meningitidis* FAM18',
    'xoo_m10_hifi': '*X. oryzae* M10',
    # final validation set
    'seon_abd4340': '*Seonamhaeicola* sp. ABD4_340', 'lsal_11254': '*Ligilactobacillus salivarius* 11254D021BM',
    'llac_s11599': '*Lactococcus lactis* S11-599', 'ppsy_bc2': '*Psychrobacillus psychrotolerans* BC2',
    'ypes_80': '*Yersinia pestis* 80', 'cpor_mei2o223': '*Citrobacter portucalensis* MEI2O2-23 (ligation)',
    'psep_tb5p01a': '*Pantoea septica* TB5-P01A', 'parth_whri8279': '*Pseudarthrobacter* sp. WHRI 8279',
    'csph_dsm158': '*Cereibacter sphaeroides* DSM 158', 'ader_km1812': '*Arsenicicoccus dermatophilus* KM18/12',
    'anae_jcm17239_hifi': '*Anaerosalibacter bizertensis* JCM 17239',
    'pmeg_imgn6_hifi': '*Priestia megaterium* IMGN6',
    'fcer_lth5742_hifi': '*Furfurilactobacillus cerevisiae* LTH 5742',
    'athe_ccm8960_hifi': '*Aneurinibacillus thermoaerophilus* CCM 8960',
    'psyc_ch211_hifi': '*Psychrobacter* sp. KFRI-CH2-11', 'kbal_lth7215_hifi': '*Kozakia baliensis* LTH 7215',
    'tcar_atcc25418_hifi': '*Trinickia caryophylli* ATCC 25418', 'rrub_dsm467_hifi': '*Rhodospirillum rubrum* DSM 467',
    'cthe_dsm15264_hifi': '*Caldimonas thermodepolymerans* DSM 15264', 'msol_sl370_hifi': '*Micromonospora soli* SL3-70',
    # depth series
    'sflex_half': '*S. flexneri*', 'kpneu_half': '*K. pneumoniae*', 'ehorm_half': '*E. hormaechei*',
    'linno_half': '*L. innocua*', 'prett_half': '*P. rettgeri*',
    'kpneu_20x': '*K. pneumoniae* 20x', 'sflex_20x': '*S. flexneri* 20x', 'ehorm_15x': '*E. hormaechei* 15x',
    'linno_12x': '*L. innocua* 12x', 'prett_12x': '*P. rettgeri* 12x',
    'kpneu_x4': '*K. pneumoniae* ~200x', 'sflex_x4': '*S. flexneri* ~200x',
}
GENOMES = {'sflex': 'Shigella_flexneri', 'kpneu': 'Klebsiella_pneumoniae', 'ehorm': 'Enterobacter_hormaechei',
           'linno': 'Listeria_innocua', 'prett': 'Providencia_rettgeri'}
TOOL_NAMES = {'picoasm': 'picoASM', 'flye': 'Flye', 'flye_medaka': 'Flye + Medaka', 'hifiasm': 'hifiasm',
              'autocycler': 'Autocycler', 'autocycler_medaka': 'Autocycler + Medaka',
              'autocycler_manual': 'Autocycler (manual curation)', 'hybracter': 'Hybracter', 'canu': 'Canu',
              'raven': 'Raven', 'miniasm': 'miniasm', 'dragonflye': 'Dragonflye', 'lja': 'LJA',
              'maeci': 'MAECI', 'necat': 'NECAT', 'metamdbg': 'metaMDBG', 'myloasm': 'Myloasm',
              'redbean': 'Redbean', 'nextdenovo': 'NextDenovo'}


# Datasets left out of every table: no reads could be downloaded for this one.
SKIP = {'xoo_whri5233'}


def errors(r):
    """The row's errors, 0 for a dataset scored for structure only."""
    return 0 if r['errors'] == 'NA' else int(r['errors'])


def read(path):
    with open(path) as f:
        return list(csv.DictReader(f, delimiter='\t'))


def cell(r):
    """`complete/replicons, errors`, then missing kb and misassemblies where there are any."""
    errors = 'n/a' if r['errors'] == 'NA' else f"{int(r['errors']):,}"
    s = f"{r['complete']}/{r['replicons']}, {errors}"
    missing, mis = int(r['missing_bases']), int(r['misassemblies'])
    if missing >= 100 or mis:
        s += f", {missing / 1000:.1f} kb" if missing >= 100 else ", 0 kb"
        if mis:
            s += f" / {mis}"
    return s


def table(results, datasets, tools, reads=None):
    print('| Isolate | ' + ' | '.join(TOOL_NAMES[t] for t in tools) + ' |')
    print('|---|' + '---|' * len(tools))
    for d in datasets:
        rows = {r['tool']: r for r in results if r['dataset'] == d and (reads is None or r['reads'] == reads)}
        if 'picoasm' not in rows:
            continue
        print(f'| {NAMES.get(d, d)} | ' + ' | '.join(cell(rows[t]) if t in rows else '' for t in tools) + ' |')
    print()


def in_set(results, *sets):
    """The datasets of `sets` that have results, in the order of their lists."""
    have = {r['dataset'] for r in results}
    seen = []
    for s in sets:
        with open(ROOT / 'scripts' / f'bact_datasets_{s}.tsv') as f:
            for line in f:
                name = line.split('\t')[0]
                if not line.startswith('#') and name in have and name not in seen and name not in SKIP:
                    seen.append(name)
    return seen


def subsets(results):
    paper = read(BENCH / 'autocycler_paper.tsv')
    by_tool = collections.defaultdict(list)
    for r in paper:
        by_tool[r['tool']].append(r)
    for r in results:
        if r['set'] == 'subsets' and r['tool'] == 'picoasm':
            g, n = r['dataset'].rsplit('_', 1)
            by_tool['picoasm'].append(dict(r, genome=GENOMES[g], subset=n))
    print('| Tool | All replicons complete | Perfect | Median errors | Total errors | Missing bases '
          '| Extra contigs | Misassemblies |')
    print('|---|---|---|---|---|---|---|---|')
    order = sorted(by_tool, key=lambda t: (statistics.median(int(r['errors']) for r in by_tool[t]),
                                           sum(int(r['errors']) for r in by_tool[t])))
    for t in order:
        rs = by_tool[t]
        complete = sum(r['complete'] == r['replicons'] for r in rs)
        perfect = sum(r['perfect'] == 'True' for r in rs)
        med = statistics.median(int(r['errors']) for r in rs)
        print(f'| {TOOL_NAMES[t]} | {complete}/{len(rs)} | {perfect} | {med:g} | '
              f"{sum(errors(r) for r in rs):,} | {sum(int(r['missing_bases']) for r in rs):,} | "
              f"{sum(int(r['extra_contigs']) for r in rs)} | {sum(int(r['misassemblies']) for r in rs)} |")
    print()
    print('Median errors per genome over its six subsets:')
    print()
    tools = ['picoasm', 'autocycler_medaka', 'autocycler', 'flye_medaka', 'flye']
    print('| Genome | ' + ' | '.join(TOOL_NAMES[t] for t in tools) + ' |')
    print('|---|' + '---|' * len(tools))
    for g in GENOMES.values():
        meds = [statistics.median(int(r['errors']) for r in by_tool[t] if r['genome'] == g) for t in tools]
        print(f"| *{g.replace('_', ' ')}* | " + ' | '.join(f'{m:g}' for m in meds) + ' |')
    print()
    pico = by_tool['picoasm']
    print(f"picoASM on the subsets: median {statistics.median(float(r['wall_s']) for r in pico):.0f} s wall, "
          f"{statistics.median(float(r['cpu_s']) for r in pico) / 60:.1f} CPU-min, median peak memory "
          f"{statistics.median(float(r['rss_gb']) for r in pico):.2f} GB (max "
          f"{max(float(r['rss_gb']) for r in pico):.2f} GB).")
    print()


def summary(results, sets, reads, tools):
    """Per tool over the datasets of `sets` with `reads` that every tool in `tools` was run on:
    datasets with every replicon complete, complete replicons, errors, missing bases."""
    ds = [d for d in in_set(results, *sets)
          if all(any(r['dataset'] == d and r['reads'] == reads and r['tool'] == t for r in results) for t in tools)]
    print(f"{', '.join(sets)} ({reads}, {len(ds)} datasets run with every tool):")
    for t in tools:
        rows = [r for r in results if r['dataset'] in ds and r['reads'] == reads and r['tool'] == t]
        complete = sum(r['complete'] == r['replicons'] for r in rows)
        reps = sum(int(r['complete']) for r in rows), sum(int(r['replicons']) for r in rows)
        wall = statistics.median(float(r['wall_s']) for r in rows)
        print(f"  {TOOL_NAMES[t]:15} all complete {complete:3}  replicons {reps[0]}/{reps[1]}  "
              f"errors {sum(errors(r) for r in rows):>9,}  missing {sum(int(r['missing_bases']) for r in rows):>10,}  "
              f"misassemblies {sum(int(r['misassemblies']) for r in rows):3}  median wall {wall:.0f} s")
    print()


def main():
    results = read(BENCH / 'results.tsv')
    print('## Autocycler paper, 30 read subsets\n')
    subsets(results)
    print('## Hall et al. 2024 isolates\n')
    table(results, in_set(results, 'unseen', 'heldout') + ['mtb', 'kpc2'], ['picoasm', 'flye', 'flye_medaka'], 'ont')
    print('## Wider validation\n')
    draft = ['seuropae', 'strep_stell', 'strep_puni']
    table(results, [d for d in in_set(results, 'wide', 'wide2', 'extra') if d not in draft], ['picoasm', 'flye_medaka'], 'ont')
    print('Draft references (missing bases and misassemblies measure the reference as much as the assembly):\n')
    table(results, draft, ['picoasm', 'flye_medaka'], 'ont')
    print('## Fresh isolates\n')
    table(results, in_set(results, 'fresh'), ['picoasm', 'flye_medaka'], 'ont')
    table(results, in_set(results, 'fresh'), ['picoasm', 'flye', 'hifiasm'], 'hifi')
    print('## Final validation set\n')
    table(results, in_set(results, 'final'), ['picoasm', 'flye_medaka'], 'ont')
    table(results, in_set(results, 'final'), ['picoasm', 'flye', 'hifiasm'], 'hifi')
    print('## Depth series\n')
    table(results, in_set(results, 'lowdepth', 'verylow', 'highdepth'), ['picoasm', 'flye', 'flye_medaka'], 'ont')
    print('## Totals\n')
    summary(results, ['final'], 'ont', ['picoasm', 'flye_medaka'])
    summary(results, ['final'], 'hifi', ['picoasm', 'flye', 'hifiasm'])
    summary(results, ['more'], 'ont', ['picoasm', 'flye_medaka'])
    summary(results, ['more'], 'hifi', ['picoasm', 'flye', 'hifiasm'])
    summary(results, ['fresh'], 'ont', ['picoasm', 'flye_medaka'])
    summary(results, ['unseen', 'heldout'], 'ont', ['picoasm', 'flye', 'flye_medaka'])
    summary(results, ['wide', 'wide2', 'extra'], 'ont', ['picoasm', 'flye_medaka'])


if __name__ == '__main__':
    main()
