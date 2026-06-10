# LAI HMM v0.4 CLI

Command-line interface and Python wrapper for local ancestry inference (LAI) from VCF variant genotypes, `hap_genotype` marker haplotype allele IDs, or both. The model uses reference population, species, clade, or group allele frequencies to infer ancestry along the genome in diploid, unphased samples.

The hidden states are all unordered diploid pairs of the requested clades. For example, four clades produce 10 possible diploid states and five clades produce 15 possible diploid states. Transitions are based on physical distance between adjacent markers.

## What this tool can run

`LAI_HMM_v0.4.py` supports four common workflows:

1. Build reference allele-frequency files from non-admixed reference samples.
2. Run LAI from VCF variants only.
3. Run LAI from `hap_genotype` haplotype allele IDs only.
4. Run LAI from both VCF variants and `hap_genotype` data, combining the two evidence sources marker-by-marker.

When both input types are provided, the model calculates variant-based and haplotype-ID-based emissions separately, then mixes them using marker-specific informativeness. 
This is useful when a haplotype allele ID is missing from the reference panel, but the variants within the haplotype sequence still contain useful information. 
In this case, haplotype evidence is down-weighted to zero when the observed haplotype allele ID is not represented in the reference lookup.

## Installation

Use Python 3.9 or newer. A conda environment is recommended, especially for VCF support through `cyvcf2`.

```bash
conda create -n lai-hmm -c conda-forge -c bioconda python=3.12 pandas numpy matplotlib cyvcf2
conda activate lai-hmm
```

Check that the command-line interface is available:

```bash
python LAI_HMM_v4.py --help
```



## Example data

The `example_data` directory includes small files for testing:

- `example.vcf.gz` and `example.vcf.gz.csi`
- `hap_genotype_example`
- `marker_positions.csv`
- `example_reference_membership.tsv`
- `reference_variant_profiles.tsv`
- `reference_hap_allele_frequencies.tsv`
- `reference_hap_allele_frequency_lookup.pkl`
- `reference_hap_allele_informativeness.tsv`
- `example.fai`



## Required and optional input files

### VCF genotype input

Important VCF notes:

- VCF should be indexed
- VCF input uses `cyvcf2`.
- VCF-only workflows may omit `INFO/MARKER`.
- When combining VCF with `hap_genotype`, VCF records must be annotated with `INFO/MARKER=<marker_name>` so variants can be grouped into the same marker names used by the haplotype matrix.
- In combined workflows, `INFO/MARKER` values should match the `hap_genotype` row names and `--marker-positions` marker names exactly after string conversion.
- If no marker position file is supplied for VCF-only runs, marker positions can be inferred from the VCF records, but a marker position file is still recommended for reproducible ordering and combined runs.

### `hap_genotype` format

A `hap_genotype` file is a marker-by-sample matrix of haplotype allele IDs.

- One row per marker.
- The first column is the marker name, usually `Locus`
- The second column, `Haplotypes`, is optional metadata that are ignored by the parser.
- Remaining columns are sample names.
- Each sample cell contains two allele IDs, such as `1/2`, `1|2`, or `1/2:read_counts`.
- Missing values may be `./.:0`, `./.`, `.`, `NA`, or blank.
- Gzipped text files are accepted when the filename ends in `.gz`.

Example:

```text
Locus	Haplotypes	SampleA	SampleB
Marker1	1(0.5);2(0.5);	1/2:10,8	2/2:12
Marker2	3(0.7);4(0.3);	./.:0	3/4:9,2
```

### Reference membership format

This is a metadata table mapping non-admixed reference samples IDs to clades, populations, species, or other groups.

Example (header optional):

```text
sample	clade
SampleA	EA
SampleB	Mus
SampleC	NA1
SampleD	Vv
```

Accepted sample column names include `sample`, `IID`, `id`, and `sample_id`. Accepted group column names include `clade`, `group`, `population`, `pop`, `species`, and `index`. 
If no recognized header is available, the first column is treated as the sample name and the second as the clade/group name.

The values in this file must match the sample names in the VCF and/or `hap_genotype` file. The clade/group names should also match the names supplied to `--clades`.

### Haplotype frequency lookup format

`--hap-frequencies` accepts either:

- the legacy pickle lookup, such as `reference_hap_allele_frequency_lookup.pkl`
- a delimited text table, such as `reference_hap_allele_frequencies.tsv` or a comma-delimited `.csv`

For text input, the file should contain one row per marker/allele pair, with:

- a marker column such as `marker`, `locus`, or `id`
- an allele column such as `allele_id` or `allele`
- one numeric column per clade listed in `--clades`

Example:

```text
marker	allele_id	EA	Mus	NA1	NA2	Vv
Marker1	1	0.95	0.05	0.00	0.00	0.00
Marker1	2	0.05	0.85	0.05	0.05	0.00
Marker2	3	0.10	0.10	0.70	0.10	0.00
```

The long-form `reference_hap_allele_frequencies.tsv` file written during `build-reference` can be fed back into `--hap-frequencies` directly.

### Marker position format

Use `--marker-positions` for a table describing marker order and physical positions.

Accepted column names are:
- marker: `marker`, `markers`, `locus`, `loci`, `id`
- chromosome: `chrom`, `chromosome`, `chr`
- position: `pos`, `position`, `bp`, `start`

Headerless three-column files are also accepted.
Example:

```text
marker,chrom,pos
Marker1,chr01,105000
Marker2,chr01,350000
Marker3,chr02,120000
```

The loader normalizes chromosome labels such as `chr01`, `chr1`, `Chromosome12`, and `CHR07` to numeric chromosome order for sorting. 
Rows on chromosome labels equivalent to `chr00`, `0`, or `00` are dropped.

### Chromosome lengths for plotting

Chromosome painting plots require `--chrom-lengths`. This can be either:

- A FASTA index file ending in `.fai` or `.fai.gz`, where column 1 is chromosome and column 2 is length.
- A delimited table with chromosome name and chromosome length columns.

Example:
```text
chrom	length
chr01	23000000
chr02	19000000
chr03	21000000
```

Plotting is enabled by default. If you do not provide chromosome lengths, include `--no-plots`.

## Quick start examples

### 1. Build reference files from VCF and hap_genotype data

Use this when you have reference samples and want to create clade-specific allele frequency files.

```bash
python LAI_HMM_v4.py \
  --step build-reference \
  --vcf example.vcf.gz \
  --hap-genotype hap_genotype_example \
  --reference-membership example_reference_membership.tsv \
  --clades EA,Mus,NA1,NA2,Vv \
  --reference-outdir ./
```

You can build only haplotype reference files by omitting `--vcf`, or only VCF reference files by omitting `--hap-genotype`.

### 2. Run the full pipeline and auto-build missing reference files

If precomputed reference files are missing and `--reference-membership` is supplied, the default `--step all` workflow builds reference files under `OUTDIR/reference` and then runs the HMM.

```bash
python LAI_HMM_v4.py \
  --vcf example.vcf.gz \
  --hap-genotype hap_genotype_example \
  --marker-positions marker_positions.csv \
  --reference-membership example_reference_membership.tsv \
  --clades EA,Mus,NA1,NA2,Vv \
  --all-samples \
  --chrom-lengths chrom_lengths.tsv \
  --outdir results
```

Add `--chrom-lengths chrom_lengths.tsv` and omit `--no-plots` if you want chromosome painting figures.

### 3. Run VCF-only LAI with precomputed reference profiles

```bash
python LAI_HMM_v4.py \
  --vcf example.vcf.gz \
  --variant-profiles reference_variant_profiles.tsv \
  --marker-positions marker_positions.csv \
  --outdir results_vcf_only
```

### 4. Run hap_genotype-only LAI with precomputed haplotype frequencies

`--hap-frequencies` accepts either the pickle lookup or the tabular `reference_hap_allele_frequencies.tsv` output from reference building.

```bash
python LAI_HMM_v4.py \
  --hap-genotype hap_genotype_example \
  --hap-frequencies reference_hap_allele_frequencies.tsv \
  --marker-positions marker_positions.csv \
  --outdir results_hap_only
```

### 5. Run combined VCF + hap_genotype LAI with four worker threads

Combined VCF + `hap_genotype` runs require haplotype informativeness scores in addition to the haplotype frequency lookup, because the model uses marker-specific informativeness to weight haplotype-ID evidence against variant evidence.

```bash
python LAI_HMM_v4.py \
  --vcf example.vcf.gz \
  --hap-genotype hap_genotype_example \
  --variant-profiles reference_variant_profiles.tsv \
  --hap-frequencies reference_hap_allele_frequency_lookup.pkl \
  --hap-informativeness reference_hap_allele_informativeness.tsv \
  --marker-positions marker_positions.csv \
  --threads 4 \
  --outdir results_combined
```

### 6. Run one sample and save chromosome painting plots

```bash
python LAI_HMM_v4.py \
  --vcf cohort.vcf.gz \
  --hap-genotype hap_genotype.tsv.gz \
  --variant-profiles reference/reference_variant_profiles.tsv \
  --hap-frequencies reference/reference_hap_allele_frequency_lookup.pkl \
  --hap-informativeness reference/reference_hap_allele_informativeness.tsv \
  --marker-positions marker_positions.csv \
  --chrom-lengths genome.fai \
  --sample SAMPLE1 \
  --plots \
  --outdir results_with_plots
```

## Workflow steps

The `--step` option controls which part of the workflow runs.

| Step | Use case | Required inputs | Main outputs |
| --- | --- | --- | --- |
| `all` | Default full pipeline. Builds missing references when `--reference-membership` is supplied, then runs the HMM. | `--vcf` and/or `--hap-genotype`; reference files or `--reference-membership`; usually `--marker-positions` | Per-sample posteriors, clade summary, log, optional plots, optional reference files |
| `build-reference` | Build clade-specific frequency files from reference samples. | `--reference-membership`; `--vcf` and/or `--hap-genotype` | Reference frequency, lookup, informativeness, and optional PCA files |
| `run-hmm` | Run only the HMM using precomputed reference files. | Genotype inputs and matching reference files | Per-sample posteriors, clade summary, log, optional plots |
| `variant-informativeness` | Add informativeness metrics to an existing variant profile table. | `--variant-profiles` | `variant_profiles_with_informativeness.tsv`, `variant_locus_informativeness.tsv` |
| `hap-informativeness` | Build haplotype allele-ID informativeness from an existing haplotype frequency lookup pickle or delimited table. | `--hap-frequencies` | `hap_allele_informativeness.tsv` |

Examples for individual steps:

```bash
# Build reference files only
python LAI_HMM_v4.py \
  --step build-reference \
  --hap-genotype hap_genotype_example \
  --reference-membership example_reference_membership.tsv \
  --clades EA,Mus,NA1,NA2,Vv \
  --reference-outdir ./

# Variant informativeness only
python LAI_HMM_v4.py \
  --step variant-informativeness \
  --variant-profiles reference/reference_variant_profiles.tsv \
  --clades EA,Mus,NA1,NA2,Vv \
  --outdir ./

# Haplotype allele-ID informativeness only
python LAI_HMM_v4.py \
  --step hap-informativeness \
  --hap-frequencies reference/reference_hap_allele_frequencies.tsv \
  --clades EA,Mus,NA1,NA2,Vv \
  --outdir ./

# HMM only, using precomputed reference files
python LAI_HMM_v4.py \
  --step run-hmm \
  --vcf example.vcf.gz \
  --hap-genotype hap_genotype_example \
  --marker-positions marker_positions.csv \
  --variant-profiles reference_variant_profiles.tsv \
  --hap-frequencies reference_hap_allele_frequency_lookup.pkl \
  --hap-informativeness reference_hap_allele_informativeness.tsv \
  --clades EA,Mus,NA1,NA2,Vv \
  --sample SampleA \
  --no-plots \
  --outdir results
```

## Command-line options

### Input and reference options

| Option | Description |
| --- | --- |
| `--vcf PATH` | Sample or cohort VCF/BCF. Indexed bgzipped VCF is recommended. |
| `--hap-genotype PATH` | `hap_genotype` matrix with marker rows and sample columns. Alias: `--hap_genotype`. |
| `--marker-positions PATH` | Marker position table with marker, chromosome, and position columns. Alias: `--marker_positions`. |
| `--variant-profiles PATH` | Clade-specific variant allele frequency/profile TSV. Alias: `--profiles`. |
| `--hap-frequencies PATH` | Clade-specific haplotype allele frequencies as either a pickle lookup or a delimited text table with marker, allele_id, and clade columns. Alias: `--hap_freq_lookup`. |
| `--reference-membership PATH` | Reference sample-to-clade table. Alias: `--reference_membership`. |
| `--reference-outdir DIR` | Directory for generated reference files. Defaults to `OUTDIR/reference` when references are built automatically. |
| `--step STEP` | One of `all`, `build-reference`, `run-hmm`, `variant-informativeness`, or `hap-informativeness`. Default: `all`. |
| `--strict-alleles PATH` | Optional pickle of strict diagnostic haplotype alleles. Alias: `--strict_alleles`. |
| `--hap-informativeness PATH` | Optional haplotype allele informativeness TSV. Required when combining VCF and `hap_genotype` evidence unless it is being built automatically. Alias: `--hap_inf`. |
| `--chrom-lengths PATH` | Optional chromosome length table or FASTA `.fai` used for plotting. Required when plots are enabled. Alias: `--chrom_lengths`. |
| `--pca` | During reference-building, run PCA on reference samples and save PC coordinates, metrics, and plots. |
| `--mus-hap-alleles PATH` | Optional Muscadine chromosome 20 diagnostic haplotype allele CSV/TSV. |
| `--nonmus-hap-alleles PATH` | Optional non-Muscadine chromosome 7 diagnostic haplotype allele CSV/TSV. |

### Sample and output options

| Option | Description |
| --- | --- |
| `--sample SAMPLE` | Sample name to run. Can be repeated or comma-separated, for example `--sample A --sample B,C`. |
| `--samples-file PATH` | Text file with one sample name per line. Blank lines and lines starting with `#` are ignored. |
| `--all-samples` | Run all selected samples. This is the default when neither `--sample` nor `--samples-file` is provided. |
| `--sample-source auto` | Default. When both VCF and `hap_genotype` are supplied and no explicit sample list is given, run only samples respreesented in both genotyping sets. |
| `--sample-source intersection` | Same practical behavior as `auto` for all-sample combined runs: use samples present in both inputs. |
| `--sample-source union` | Run samples present in either input. Samples found in only one input use the available evidence source. |
| `--outdir DIR` | Output directory. Default: `results`. |
| `--threads N` | Number of worker threads for per-sample parallel runs. Default: `1`. |
| `--plots` | Save chromosome painting figures. Enabled by default. Requires `--chrom-lengths`. |
| `--no-plots` | Skip chromosome painting figures. |
| `--show-plots` | Display plots interactively while saving. |
| `--no-save` | Run without writing CSV, plot, or log outputs. |

### Model and tuning options

| Option | Default | Description |
| --- | ---: | --- |
| `--clades EA,Mus,NA1,NA2,Vv` | `EA,Mus,NA1,NA2,Vv` | Comma-separated clade/group names in the order used by the reference profiles. |
| `--contigs chr01,chr02` | all VCF contigs | Limit VCF loading to selected contigs. |
| `--allow-missing-marker` | off | Legacy compatibility flag for VCF-only workflows. Combined VCF + `hap_genotype` runs still require `INFO/MARKER`. |
| `--lam-per-Mb FLOAT` | `0.05` | HMM switch rate per megabase. Higher values allow more frequent ancestry switches. Alias: `--lam_per_Mb`. |
| `--strict-boost FLOAT` | `5.0` | Multiplicative boost for strict diagnostic haplotype alleles. Must be at least `1.0` when strict allele data are supplied. Alias: `--strict_boost`. |
| `--cap-total-boost-per-marker FLOAT` | none | Optional cap on the total strict-allele boost per marker. Alias: `--cap_total_boost_per_marker`. |
| `--trans-temp FLOAT` | `1.0` | Transition temperature. Values above 1 flatten transition preferences; values below 1 sharpen them. Alias: `--trans_temp`. |
| `--alpha FLOAT` | `0.0` | Context smoothing blend. `0` disables smoothing; larger values blend marker emissions with neighboring-marker context. |
| `--windowsize INT` | `0` | Context smoothing half-window size. `0` disables smoothing. |
| `--e-geno FLOAT` | `0.01` | Genotyping or miscoding error rate used in emission probabilities. Alias: `--e_geno`. |
| `--e-homo FLOAT` | `0.1` | Homozygote dropout softening. A fraction of observed homozygotes is treated as potentially undercalled heterozygotes before any informativeness-based homozygote flattening. `0` disables this dropout softening. Alias: `--e_homo`. |
| `--b0 FLOAT` | `0.0` | Baseline tilt when mixing haplotype-ID and variant evidence. Positive values favor haplotype-ID evidence; negative values favor variant evidence. |
| `--certainty-softener FLOAT` | `1.0` | Posterior softening temperature. `1` disables softening; larger values flatten posterior certainty. Alias: `--certainty_softener`. |
| `--hom-soften-delta FLOAT` | `0.1` | Separation threshold for the `soften_mix_for_homozygote` step. Larger values require a bigger top-versus-second clade frequency gap before a homozygote is trusted strongly. Alias: `--hom_soften_delta`. |
| `--hom-soften-width FLOAT` | `0.6` | Width/slope around `--hom-soften-delta`. Larger values make the switch from softened to trusted homozygotes more gradual. Alias: `--hom_soften_width`. |
| `--hom-min-mix FLOAT` | `0.6` | Minimum weight retained on the dropout-adjusted homozygote emission during `soften_mix_for_homozygote`. Lower values flatten weak homozygotes more aggressively; `1` disables this separation-based softening. Alias: `--hom_min_mix`. |
| `--hom-neutral FLOAT` | `0.60` | Neutral homozygote probability target used when homozygous evidence is softened. Alias: `--hom_neutral`. |
| `--tau FLOAT` | `1.0` | Posterior-guided Viterbi blend. `0` uses regular Viterbi emissions; larger values blend emissions with posterior probabilities. |

### Logging options

| Option | Description |
| --- | --- |
| `--quiet` | Suppress progress messages. |
| `--debug` | Print extra diagnostic information. |

## Reference outputs

Reference-building writes files with the prefix `reference` by default.

For `hap_genotype` reference input:

- `reference_hap_allele_frequencies.tsv`: long table of clade-specific haplotype allele-ID frequencies. This file can be used directly with `--hap-frequencies`.
- `reference_hap_allele_frequency_lookup.pkl`: pickle lookup used by HMM runs. This remains supported as a compact legacy format.
- `reference_hap_allele_informativeness.tsv`: haplotype allele-ID informativeness scores.
- If `--pca` is used: `reference_hap_pca.tsv`, `reference_hap_pca.png`, and `reference_hap_differentiation_metrics.json`.

For VCF reference input:

- `reference_variant_profiles.tsv`: clade-specific VCF allele frequency/profile table with informativeness columns.
- `reference_variant_locus_informativeness.tsv`: marker/locus-level variant informativeness.
- If `--pca` is used: `reference_vcf_pca.tsv`, `reference_vcf_pca.png`, and `reference_vcf_differentiation_metrics.json`.

The differentiation metrics includes explained variance and clade differentiation metrics such as mean pairwise centroid distance, mean within-clade distance, and the centroid-to-within-clade distance ratio.

## HMM outputs

The HMM wrapper writes the following outputs unless `--no-save` is used:

- `posteriors/<sample>_<YYYY-MM-DD>.csv`: per-marker HMM calls and posterior probabilities.
- `clade_percentage_summary_<YYYY-MM-DD>.csv`: one row per sample with ancestry percentages, missingness, warnings, or per-sample errors.
- `run_<timestamp>.log`: run paths, clades, and model parameters.
- `figures/<sample>_<date>.png`: chromosome painting figure, written only when plotting is enabled and chromosome lengths are supplied.

Per-marker posterior files includes columns:

- `Marker`
- `UnorderedCall` (most likely pair of clades/groups)
- `DiplPostProb` (posterior probability of the most likely clade pair)
- `DiplCertainty` (certainty (scaled posterior probability) of the most likely clade pair)
- `Hap1Clade` (clade/group selected for haplotype 1)
- `Hap2Clade` (clade/group selected for haplotype 2)
- `Chrom`
- `Pos`
- `MissingMarkers` (reports whether data was missing at that markers)
- one posterior-probability column for each unordered diploid state, such as `{EA,EA}`, `{EA,Mus}`, and `{Mus,Vv}`

The clade percentage summary calculates genome-length-weighted ancestry percentages using interval length and posterior support.

## Tuning guidance

The default settings keep optional softening and smoothing mostly off:

```text
--alpha 0
--windowsize 0
--e-homo 0
--hom-min-mix 1
--certainty-softener 1
--tau 0
```

Suggested starting points:

### Default for high-confidence genotypes

Use defaults.

### Low-depth data or likely heterozygote undercalling (recommended for DArTag data)

Increase homozygote/dropout softening and optionally soften final posteriors.

```bash
python LAI_HMM_v4.py ... \
  --e-homo 0.1 \
  --hom-min-mix 0.6 \
  --certainty-softener 2 
```

### Noisy isolated marker calls

Enable local context smoothing across groups of markers (`--windowsize`). Use small values first. 
The `--alpha` parameter controls the degree of smoothing across the window.

```bash
python LAI_HMM_v4.py ... \
  --alpha 0.3 \
  --windowsize 5
```

### More frequent or less frequent ancestry switches

Adjust `--lam-per-Mb`:

- Increase `--lam-per-Mb` when true ancestry switches are expected to be frequent.
- Decrease `--lam-per-Mb` for smoother ancestry tracts.

```bash
python LAI_HMM_v4.py ... --lam-per-Mb 0.02
```

### Rebalance haplotype-ID versus variant evidence

When both VCF and `hap_genotype` inputs are supplied, `--b0` shifts the marker-specific mixture:

- Positive `--b0` gives more baseline weight to haplotype-ID evidence.
- Negative `--b0` gives more baseline weight to variant evidence.

```bash
python LAI_HMM_v4.py ... --b0 0.5 --no-plots
```

### Tune homozygote tempering

The HMM has two separate mechanisms for preventing homozygous observations from becoming overconfident too quickly:

- `--e-homo` handles possible technical undercalling. When the observed genotype is homozygous, the model first mixes in some probability that the true genotype was heterozygous but one allele dropped out or was undercalled. Increase `--e-homo` when low coverage, allele dropout, or assay bias makes false homozygotes plausible.
- The internal `soften_mix_for_homozygote()` step handles weak ancestry informativeness. This is controlled by `--hom-soften-delta`, `--hom-soften-width`, `--hom-min-mix`, and `--hom-neutral`. It does not assume the genotype call is wrong. Instead, it flattens homozygous evidence when the observed allele is common across multiple clades and therefore not very ancestry-specific.

Use `--e-homo` when the main concern is technical miscalling of heterozygotes as homozygotes.
Use the `--hom-*` options when genotype calls look technically fine, but shared or weakly diagnostic alleles are still producing ancestry posteriors that feel too confident.
Use both when both failure modes are present, because the model applies them in sequence: dropout mixing first, then informativeness-based softening.

Practical effects of the `soften_mix_for_homozygote()` controls:

- Increase `--hom-soften-delta` to demand stronger clade specificity before a homozygote is trusted strongly.
- Increase `--hom-soften-width` to make the transition from softened to trusted evidence more gradual.
- Decrease `--hom-min-mix` to flatten weak homozygotes more aggressively; setting it to `1` disables this separation-based flattening.
- Adjust `--hom-neutral` to change the neutral homozygote target that weak loci are mixed toward.

## Python API examples

The same code can be used from Python.

### Build reference files

```python
from LAI_HMM_v4 import build_reference_files

outputs = build_reference_files(
    outdir="reference",
    clades=["EA", "Mus", "NA1", "NA2", "Vv"],
    membership_path="example_reference_membership.tsv",
    hap_genotype_path="hap_genotype_example",
    vcf_path="example.vcf.gz",
    pca=True,
    verbose=True,
)

print(outputs)
```

### Run the HMM

```python
from LAI_HMM_v4 import run_lai_hmm

result = run_lai_hmm(
    vcf_path="example.vcf.gz",
    hap_genotype_path="hap_genotype_example",
    marker_positions_path="marker_positions.csv",
    variant_profiles_path="reference/reference_variant_profiles.tsv",
    hap_freq_lookup_path="reference/reference_hap_allele_frequencies.tsv",
    hap_informativeness_path="reference/reference_hap_allele_informativeness.tsv",
    clades=["EA", "Mus", "NA1", "NA2", "Vv"],
    sample="SampleA",
    outdir="results",
    make_plots=False,
)

print(result["summary_df"])
```

## Troubleshooting

### Missing reference files

For VCF input, provide `--variant-profiles` or supply `--reference-membership` so the CLI can build variant profiles.

For `hap_genotype` input, provide `--hap-frequencies` or supply `--reference-membership` so the CLI can build haplotype frequencies.

For combined VCF + `hap_genotype` input, also provide `--hap-informativeness` unless it is being built automatically.

### No VCF loci collected

Likely causes:

- VCF contig names do not match `--contigs`.
- In combined VCF + `hap_genotype` runs, VCF records are missing `INFO/MARKER`.
- The VCF index is missing or incompatible.

Try removing `--contigs`, checking VCF header contig names, and, for combined runs, confirming that `INFO/MARKER` is present and matches the haplotype marker names.

### Many missing markers

Check marker names across all files:

- `marker_positions.csv`
- VCF `INFO/MARKER` for combined runs
- `hap_genotype` row names
- reference frequency files

Marker names must match exactly after string conversion.

### Unexpected clades or state columns

Ensure `--clades` exactly matches the group names and order in the reference frequency files. The HMM creates one state for every unordered diploid pair of clades.

### Requested sample not found

Check that the sample name appears in at least one genotype input. When both VCF and `hap_genotype` inputs are present, explicitly requested samples may appear in either file, but all-sample combined runs use the intersection by default unless `--sample-source union` is specified.

### Combined run fails because haplotype informativeness is missing

Provide the informativeness file generated during reference-building:

```bash
--hap-informativeness reference/reference_hap_allele_informativeness.tsv
```

or regenerate it from the frequency lookup:

```bash
python LAI_HMM_v4.py \
  --step hap-informativeness \
  --hap-frequencies reference/reference_hap_allele_frequency_lookup.pkl \
  --clades EA,Mus,NA1,NA2,Vv \
  --outdir reference
```
