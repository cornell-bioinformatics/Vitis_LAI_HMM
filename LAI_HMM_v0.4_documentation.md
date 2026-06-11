# LAI HMM v0.4 CLI

## Introduction

This is a command-line interface and Python wrapper for local ancestry inference (LAI) from amplicon-based, multi-allelic marker data. The model uses non-admixed reference population allele frequencies to infer ancestry along the genome in diploid, unphased samples, accounting for genotyping error and heterozygous undercalling. Ancestry can be inferred for any populations/species/clades or other groupings, as long as there are a sufficient number of non-admixed samples for each clade/group with which to estimate the clade-specific allele frequencies. The chromosome painting pipeline outputs visualizations of the inferred ancestry tracts across the genome.

The hidden states in the hidden Markov model are all unordered diploid pairs of the specified clades/groups (e.g., four clades produce 10 possible diploid states, five clades produce 15 possible diploid states). The fewer the hidden states, the easier it is to infer. Transitions are based on physical distance between adjacent markers, with a prior switch rate that can be tuned according to the assumptions about the recombination frequency in your population.

The model also has tunable parameters to account for genotyping errors and heterozygous undercalling (due to allele drop-out at low sequencing depth, deletions, PCR non-amplification, etc.) which is common (~10%) in rhAmpSeq data and (possibly to a lesser extent) in DArTag data. Ideally, the pipeline should be tested on samples with known ancestry (e.g., F1 hybrids or previously characterized samples) so that these parameters can be tuned.

The pipeline can be parallelized to run on multiple threads for large sample sets.

## What this tool can run

`LAI_HMM_v0.4.py` supports four workflows:

1. Build reference allele-frequency files from non-admixed reference samples.
2. Run LAI from VCF variants only.
3. Run LAI from `hap_genotype`-format microhaplotype allele IDs only.
4. Run LAI from both VCF variants and `hap_genotype` data, combining the two evidence sources.


Combining both input types is recommended because the multiallelic microhaplotypes tend to be more informative of ancestry, but if a microhaplotype allele ID is missing from the reference panel, the individual variant sites within the microhaplotype sequence still contain useful information. When both input types are provided, the model calculates variant-based and haplotype-ID-based emissions separately, then mixes them using marker-specific informativeness.


## Installation

Use Python 3.9 or newer. A conda environment is recommended, especially for VCF support through `cyvcf2`.

```bash
conda create -n lai-hmm -c conda-forge -c bioconda python=3.12 pandas numpy matplotlib cyvcf2
conda activate lai-hmm
```

Check that the command-line interface is available:

```bash
python LAI_HMM_v0.4.py --help
```



## Example data

Example inputs are in the `example_data/` directory for testing the pipeline:

**Reference-panel inputs:**

- `example_data/example_refset.vcf.gz` and `example_data/example_refset.vcf.gz.csi`
- `example_data/hap_genotype_refset_example`
- `example_data/example_reference_membership.tsv`

**Sample inputs to infer ancestry for:**

- `example_data/example_samples.vcf.gz` and `example_data/example_samples.vcf.gz.csi`
- `example_data/hap_genotype_samples_example`

**Genome coordinate inputs:**

- `example_data/marker_positions.csv`
- `example_data/chrom_lengths.fai`

**Precomputed clade-specific allele frequency and informativeness files:**

- `example_data/reference_variant_profiles.tsv`
- `example_data/reference_hap_allele_informativeness.tsv`



## Required and optional input files

#### VCF genotype input

Important VCF notes:

- The VCF must be indexed and may be gzipped. 
- Ideally, VCFs should be annotated with an `INFO/MARKER=<marker_name>` field that specifies the amplicon marker name/ID, but VCF-only workflows may omit the `INFO/MARKER` field.
- When `INFO/MARKER` is present on any of the selected VCF records, all variants sharing the same marker name/ID are combined into one marker-level ancestry estimate. This helps integrate information across the (noisier) individual variant sites within an amplicon. In that mode, records lacking `INFO/MARKER` are omitted.
- If none of the selected VCF records use `INFO/MARKER` in a VCF-only run, ancestry is inferred at each individual variant position using a fallback marker label of `CHROM:POS`.
- VCF records must be annotated with `INFO/MARKER` when combining the VCF input and the `hap_genotype` microhaplotypes input. The `INFO/MARKER` values should match the `hap_genotype` row names and `--marker-positions` marker names exactly.

#### `hap_genotype` format

A `hap_genotype` file is a marker-by-sample matrix of microhaplotype allele IDs.

- One row per marker.
- The first column is the marker name, usually `Locus`
- The second column, `Haplotypes`, is optional metadata that is ignored by the parser.
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

#### Reference membership format

This text file maps non-admixed reference samples IDs to their corresponding clades/species/groups. It is used to calculate the clade-specific allele frequencies, and may be omitted if those frequencies have been calculated already.

Example (header optional):

```text
sample	clade
SampleA	EA
SampleB	Mus
SampleC	NA1
SampleD	Vv
```

Accepted sample column names include `sample`, `IID`, `id`, and `sample_id`. Accepted group column names include `clade`, `group`, `population`, `pop`, `species`, and `index`. If no recognized header is available, the first column is treated as the sample name and the second as the clade/group name.

The values in this file must match the sample names in the VCF and/or `hap_genotype` file. The clade/group names should also match the names supplied to `--clades`.

#### Haplotype frequency lookup format

This file specifies the clade-specific allele frequencies. It may be provided by the user or built using the `build-reference` step.

`--hap-frequencies` accepts either:

- a delimited text table, such as `reference_hap_allele_frequencies.tsv` or a comma-delimited `.csv`
- a pickled lookup, such as `reference_hap_allele_frequency_lookup.pkl`

For text input, the file should contain one row per marker/allele pair, with:

- a marker column such as `marker`, `locus`, or `id`
- an allele column such as `allele_id` or `allele`
- one numeric column per clade listed in `--clades`

Example:

```text
marker	allele_id	EA	Mus	NA1	NA2	Vv
Marker1	1	0.95	0.05	0.00	0.00	0.00
Marker1	2	0.05	0.85	0.05	0.05	0.00
Marker2	1	0.10	0.10	0.70	0.10	0.00
```

The long-form `reference_hap_allele_frequencies.tsv` file written during `build-reference` can be fed back into `--hap-frequencies` directly.


#### Haplotype allele informativeness format

Combined VCF + `hap_genotype` runs require haplotype informativeness scores such as `reference_hap_allele_informativeness.tsv` because the model uses marker-specific informativeness to weight haplotype-ID evidence against variant evidence.

Informativeness scores are calculated during the reference building step or can be calculated separately with `--step hap-informativeness` (see below). Both specific and normalized (0-1) scores are given.

Allele "informativeness of assignment" is calculated as in Rosenberg et al. 2003. 
If an allele is fixed in one clade and absent from all others, it is maximally informative of ancestry at that locus (normalized informativeness score of 1). 
If an allele is at the same frequency across all clades, then it is non-informative (normalized informativeness score of 0).


Example:

```text
marker	allele_id	specific	normalized_Inf
Marker1	1	1.45	0.905
Marker2	10	1.41	0.877
```


#### Marker position format

The `--marker-positions` input requires a text file describing marker order and physical positions.

Accepted column names are:

- marker: `marker`, `markers`, `locus`, `loci`, `id`
- chromosome: `chrom`, `chromosome`, `chr`
- position: `pos`, `position`, `bp`, `start`. 

Headerless three-column files are also accepted.

Common chromosome name prefixes such as `chr01`, `chr1`, `Chromosome12`, and `CHR07` are automatically normalized to numeric chromosome number for sorting. 

Example:

```text
marker,chrom,pos
Marker1,chr01,105000
Marker2,chr01,350000
Marker3,chr02,120000
```



#### Chromosome lengths

Generating the chromosome painting plots requires `--chrom-lengths`. This can be either:

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

### 1. Calculate clade-specific allele frequencies from reference samples

`--step build-reference` is a helper to create clade-specific allele frequency files from VCF and/or hap_genotype files.
You can build only haplotype reference files by omitting `--vcf`, or only VCF reference files by omitting `--hap-genotype`.
The optional `--pca` option runs principal component analysis of the reference samples and computes population differentiation metrics to help evaluate how distinct the clades/groups are.
Use `--pca` by itself to run PCA on whichever reference input types were provided, or specify `--pca hap`, `--pca vcf`, or `--pca both` to control the source explicitly.
The optional `--downsample` argument randomly keeps only a specified proportion of PCA feature columns before PCA, which can reduce runtime on large matrices. The default is `1` (use the full matrix). The same `--pca` and `--downsample` options also apply when reference files are built automatically during the default `--step all` workflow.

```bash
python LAI_HMM_v0.4.py \
  --step build-reference \
  --vcf example_data/example_refset.vcf.gz \
  --hap-genotype example_data/hap_genotype_refset_example.gz \
  --reference-membership example_data/example_reference_membership.tsv \
  --clades EA,Mus,NA,Vv \
  --reference-outdir reference \
  --pca both \
  --downsample 0.5
```


### 2. Run VCF-only LAI with precomputed allele frequencies

VCF-based runs require `--variant-profiles` which gives variant allele frequencies and informativeness scores.

```bash
python LAI_HMM_v0.4.py \
  --vcf example_data/example_samples.vcf.gz \
  --variant-profiles example_data/reference_variant_profiles.tsv \
  --marker-positions example_data/marker_positions.csv \
  --chrom-lengths example_data/chrom_lengths.fai \
  --clades EA,Mus,NA1,NA2,Vv \
  --outdir results_vcf_only
```

### 3. Run hap_genotype-only LAI with precomputed haplotype frequencies

Hap_genotype-based runs require `--hap-frequencies` and `hap-informativeness` which give microhaplotype allele frequencies and informativeness scores, respectively.

```bash
python LAI_HMM_v0.4.py \
  --hap-genotype example_data/hap_genotype_samples_example \
  --hap-frequencies example_data/reference_hap_allele_frequencies.tsv \
  --hap-informativeness example_data/reference_hap_allele_informativeness.tsv \
  --marker-positions example_data/marker_positions.csv \
  --chrom-lengths example_data/chrom_lengths.fai \
  --clades EA,Mus,NA1,NA2,Vv \
  --outdir results_hap_only
```

### 4. Run combined VCF + hap_genotype LAI for individual sample(s)

The `--sample` option can be used to specify an individual sample ID or a comma-separated list of sample IDs to run the pipeline on.
The `--samples-file` option can be used to input a list of samples.
The default, `--all-samples`, will run the pipeline on all of the samples present in both the VCF and hap_genotype file.


```bash
python LAI_HMM_v0.4.py \
  --vcf example_data/example_samples.vcf.gz \
  --hap-genotype example_data/hap_genotype_samples_example \
  --variant-profiles example_data/reference_variant_profiles.tsv \
  --hap-frequencies example_data/reference_hap_allele_frequencies.tsv \
  --hap-informativeness example_data/reference_hap_allele_informativeness.tsv \
  --marker-positions example_data/marker_positions.csv \
  --chrom-lengths example_data/chrom_lengths.fai \
  --clades EA,Mus,NA1,NA2,Vv \
  --sample SAMPLE1 \
  --outdir results
```

### 5. Run combined VCF + hap_genotype LAI for all samples with multiple threads

```bash
python LAI_HMM_v0.4.py \
  --vcf example_data/example_samples.vcf.gz \
  --hap-genotype example_data/hap_genotype_samples_example \
  --variant-profiles example_data/reference_variant_profiles.tsv \
  --hap-frequencies example_data/reference_hap_allele_frequencies.tsv \
  --hap-informativeness example_data/reference_hap_allele_informativeness.tsv \
  --marker-positions example_data/marker_positions.csv \
  --chrom-lengths example_data/chrom_lengths.fai \
  --clades EA,Mus,NA1,NA2,Vv \
  --threads 4 \
  --outdir results_combined
```


### 6. Run the full end-to-end pipeline including building reference files

If precomputed clade-specific allele frequency files are missing and `--reference-membership` is supplied, the default `--step all` workflow first builds reference files under `OUTDIR/reference` and then runs the HMM. 
In this case, the reference samples must exist in the same VCF/hap_genotype file as the test samples. If you also provide `--pca` and optionally `--downsample`, those settings are applied to the automatic reference-building step before the HMM run.

```bash
python LAI_HMM_v0.4.py \
  --vcf example_data/example_refset.vcf.gz \
  --hap-genotype example_data/hap_genotype_refset_example.gz \
  --reference-membership example_data/example_reference_membership.tsv \
  --all-samples \
  --marker-positions example_data/marker_positions.csv \
  --chrom-lengths example_data/chrom_lengths.fai \
  --clades EA,Mus,NA1,NA2,Vv \
  --outdir results
```


## Run individual steps of the workflow

The `--step` option controls which part of the workflow runs.

| Step | Description | Required inputs | Main outputs |
| --- | --- | --- | --- |
| `all` | Default full pipeline. Builds missing references when `--reference-membership` is supplied, then runs the HMM. | `--vcf` and/or `--hap-genotype`; reference files or `--reference-membership`; `--marker-positions` | Per-sample posteriors, clade summary, log, optional plots, optional reference files |
| `build-reference` | Build clade-specific allele frequency files from reference samples. | `--reference-membership`; `--vcf` and/or `--hap-genotype` | Reference frequency, lookup, informativeness, and optional PCA files |
| `run-hmm` | Run only the HMM using precomputed reference files. | Genotype inputs and matching reference files | Per-sample posteriors, clade summary, log, optional plots |
| `variant-informativeness` | Add informativeness metrics to an existing variant profile table. | `--variant-profiles` | `variant_profiles_with_informativeness.tsv`, `variant_locus_informativeness.tsv` |
| `hap-informativeness` | Build haplotype allele-ID informativeness from an existing haplotype frequency lookup pickle or delimited table. | `--hap-frequencies` | `hap_allele_informativeness.tsv` |

Examples for individual steps:

```bash
# Build reference files only
python LAI_HMM_v0.4.py \
  --step build-reference \
  --hap-genotype example_data/hap_genotype_refset_example.gz \
  --reference-membership example_data/example_reference_membership.tsv \
  --clades EA,Mus,NA1,NA2,Vv \
  --reference-outdir reference \
  --pca hap

# Calculate variant informativeness only from existing variant allele frequency information
python LAI_HMM_v0.4.py \
  --step variant-informativeness \
  --variant-profiles example_data/reference_variant_profiles.tsv \
  --clades EA,Mus,NA1,NA2,Vv \
  --outdir example_output

# Calculate microhaplotype allele informativeness only
python LAI_HMM_v0.4.py \
  --step hap-informativeness \
  --hap-frequencies example_output/reference/reference_hap_allele_frequencies.tsv \
  --clades EA,Mus,NA1,NA2,Vv \
  --outdir example_output

# LAI only (without generating chromosome painting figures), using precomputed reference files
python LAI_HMM_v0.4.py \
  --step run-hmm \
  --vcf example_data/example_samples.vcf.gz \
  --hap-genotype example_data/hap_genotype_samples_example \
  --marker-positions example_data/marker_positions.csv \
  --variant-profiles example_data/reference_variant_profiles.tsv \
  --hap-frequencies example_data/reference_hap_allele_frequencies.tsv \
  --hap-informativeness example_data/reference_hap_allele_informativeness.tsv \
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
| `--hap-genotype PATH` | `hap_genotype` matrix with marker rows and sample columns. |
| `--marker-positions PATH` | Marker position table with marker, chromosome, and position columns. |
| `--variant-profiles PATH` | Clade-specific variant allele frequency/profile TSV. |
| `--hap-frequencies PATH` | Clade-specific haplotype allele frequencies as either a delimited text table with marker, allele_id, and clade columns, or a pickle lookup. |
| `--reference-membership PATH` | Reference sample-to-clade table. |
| `--reference-outdir DIR` | Directory for generated reference files. Defaults to `OUTDIR/reference` when references are built automatically. |
| `--step STEP` | One of `all`, `build-reference`, `run-hmm`, `variant-informativeness`, or `hap-informativeness`. Default: `all`. |
| `--strict-alleles PATH` | Optional pickle of strict diagnostic haplotype alleles. |
| `--hap-informativeness PATH` | Optional haplotype allele informativeness TSV. Required when combining VCF and `hap_genotype` evidence unless it is being built automatically. |
| `--chrom-lengths PATH` | Chromosome length table or FASTA `.fai` used for plotting. Required when plots are enabled. |
| `--pca [auto\|hap\|vcf\|both]` | During reference-building, run PCA on reference samples and save PC coordinates, metrics, and plots. `--pca` by itself defaults to whichever of `--hap-genotype` and/or `--vcf` were provided. `auto` means “use whichever reference inputs were supplied.” `both` requires both `--hap-genotype` and `--vcf`. |
| `--downsample FLOAT` | During reference-building PCA, randomly retain this proportion of feature columns before PCA. Default: `1`. Example: `--downsample 0.5` keeps about half of the features. Values must be greater than `0` and less than or equal to `1`. |
| `--mus-hap-alleles PATH` | Optional Muscadine chromosome 20 diagnostic haplotype allele CSV/TSV. |
| `--nonmus-hap-alleles PATH` | Optional non-Muscadine chromosome 7 diagnostic haplotype allele CSV/TSV. |

### Sample and output options

| Option | Description |
| --- | --- |
| `--sample SAMPLE` | Sample name to run. Can be a comma-separated list, for example `--sample A` or  `--sample B,C`. |
| `--samples-file PATH` | Text file with one sample name per line. Blank lines and lines starting with `#` are ignored. |
| `--all-samples` | Run all selected samples. This is the default when neither `--sample` nor `--samples-file` is provided. |
| `--sample-source auto` | Default. When both VCF and `hap_genotype` are supplied and no explicit sample list is given, run only samples represented in both genotyping sets. |
| `--sample-source intersection` | Same practical behavior as `auto` for all-sample combined runs: use samples present in both inputs. |
| `--sample-source union` | Run samples present in either input. Samples found in only one input use the available evidence source. |
| `--outdir DIR` | Output directory. Default: `results`. |
| `--threads N` | Number of worker threads for per-sample parallel runs. Default: `1`. |
| `--plots` | Save chromosome painting figures. Enabled by default. Requires `--chrom-lengths`. |
| `--no-plots` | Skip chromosome painting figures. |
| `--no-save` | Run without writing CSV, plot, or log outputs. |

### Model and tuning options

| Option | Default | Description |
| --- | ---: | --- |
| `--clades STR` | `None` | Comma-separated clade/group names **in the order used by the reference profiles**. |
| `--contigs chr01,chr02` | all VCF contigs | Limit VCF loading to selected contigs. |
| `--lam-per-Mb FLOAT` | `0.05` | HMM switch rate per megabase. Higher values allow more frequent ancestry switches. |
| `--strict-boost FLOAT` | `5.0` | Multiplicative boost for strict diagnostic haplotype alleles. Must be at least `1.0` when strict allele data are supplied. |
| `--cap-total-boost-per-marker FLOAT` | none | Optional cap on the total strict-allele boost per marker. |
| `--trans-temp FLOAT` | `1.0` | Transition temperature. Values above 1 flatten transition preferences; values below 1 sharpen them. |
| `--alpha FLOAT` | `0.0` | Context smoothing blend. `0` disables smoothing; larger values blend marker emissions with neighboring-marker context. |
| `--windowsize INT` | `0` | Context smoothing half-window size. `0` disables smoothing. |
| `--e-geno FLOAT` | `0.01` | Genotyping error rate used in emission probabilities. |
| `--e-homo FLOAT` | `0.1` | Homozygote softening. A fraction of observed homozygotes is treated as potentially undercalled heterozygotes before any informativeness-based homozygote flattening. `0` disables this softening. |
| `--b0 FLOAT` | `0.0` | Baseline tilt when mixing haplotype-ID and variant evidence. Positive values favor haplotype-ID evidence; negative values favor variant evidence. |
| `--certainty-softener FLOAT` | `1.0` | Posterior softening temperature. `1` disables softening; larger values flatten posterior certainty (useful for adjusting opacity range for graphing purposes) |
| `--hom-soften-delta FLOAT` | `0.1` | Separation threshold for the `soften_mix_for_homozygote` step. Larger values require a bigger top-versus-second clade frequency gap before a homozygote is trusted strongly. |
| `--hom-soften-width FLOAT` | `0.6` | Width/slope around `--hom-soften-delta`. Larger values make the switch from softened to trusted homozygotes more gradual. |
| `--hom-min-mix FLOAT` | `0.6` | Minimum weight retained on the dropout-adjusted homozygote emission during `soften_mix_for_homozygote`. Lower values flatten weak homozygotes more aggressively; `1` disables this separation-based softening. |
| `--hom-neutral FLOAT` | `0.60` | Neutral probability to mix with when homozygous evidence is softened. |
| `--tau FLOAT` | `1.0` | Posterior-guided Viterbi blend. `0` uses regular Viterbi emissions; larger values blend emissions with posterior probabilities. |

### Logging options

| Option | Description |
| --- | --- |
| `--verbose` | Print progress messages. This is the default. |
| `--quiet` | Suppress progress messages. |
| `--debug` | Print extra diagnostic information. |

## Reference outputs

Reference-building writes files with the prefix `reference` by default.

For `hap_genotype` reference input:

- `reference_hap_allele_frequencies.tsv`: long table of clade-specific haplotype allele-ID frequencies. This file can be used directly with `--hap-frequencies`.
- `reference_hap_allele_informativeness.tsv`: haplotype allele-ID informativeness scores.
- If `--pca` includes haplotype PCA (`--pca`, `--pca auto`, `--pca hap`, or `--pca both` when `--hap-genotype` is supplied): `reference_hap_pca.tsv`, `reference_hap_pca.png`, and `reference_hap_differentiation_metrics.json`.

For VCF reference input:

- `reference_variant_profiles.tsv`: clade-specific VCF allele frequency/profile table with informativeness columns.
- `reference_variant_locus_informativeness.tsv`: marker/locus-level variant informativeness.
- If `--pca` includes VCF PCA (`--pca`, `--pca auto`, `--pca vcf`, or `--pca both` when `--vcf` is supplied): `reference_vcf_pca.tsv`, `reference_vcf_pca.png`, and `reference_vcf_differentiation_metrics.json`.

The differentiation metrics file includes the number of samples, the number of PCA features actually used after any `--downsample` filtering, explained variance, and clade differentiation metrics: mean pairwise centroid distance, mean within-clade distance, and the centroid-to-within-clade distance ratio.

## HMM outputs

The HMM wrapper writes the following outputs unless `--no-save` is used:

- `posteriors/<sample>_<YYYY-MM-DD>.csv`: per-marker HMM calls and posterior probabilities.
- `clade_percentage_summary_<YYYY-MM-DD>.csv`: one row per sample with ancestry percentages, missingness, warnings, or per-sample errors.
- `run_<timestamp>.log`: run paths, clades, and model parameters.
- `figures/<sample>_<date>.png`: chromosome painting figure.

Per-marker posterior files include columns:

- `Marker`
- `UnorderedCall` (most likely pair of clades/groups)
- `DiplPostProb` (posterior probability of the most likely clade pair)
- `DiplCertainty` (certainty (scaled posterior probability) of the most likely clade pair)
- `Hap1Clade` (clade/group selected for haplotype 1)
- `Hap2Clade` (clade/group selected for haplotype 2)
- `Chrom`
- `Pos`
- `MissingMarkers` (reports whether data was missing at that marker)
- one posterior-probability column for each unordered diploid state, such as `{EA,EA}`, `{EA,Mus}`, and `{Mus,Vv}`

The clade percentage summary calculates genome-length-weighted ancestry percentages using interval length and posterior support.

## Tuning guidance

The current defaults disable context smoothing and apply moderate homozygote tempering and posterior-guided Viterbi blending:

```text
--alpha 0
--windowsize 0
--e-homo 0.1
--hom-min-mix 0.6
--certainty-softener 1
--tau 1
```

Suggested starting points:

### Default for high-confidence genotypes

Use defaults.

### Low-depth data and more likely heterozygote undercalling 

Increase homozygote/dropout softening and optionally soften final posteriors.

```bash
python LAI_HMM_v0.4.py ... \
  --e-homo 0.1 \
  --hom-min-mix 0.5
  --hom_soften_delta 0.5
```

### Noisy isolated marker calls

Enable local context smoothing across groups of markers (`--windowsize`). Use small values first. 
The `--alpha` parameter controls the degree of smoothing across the window.

```bash
python LAI_HMM_v0.4.py ... \
  --alpha 0.3 \
  --windowsize 5
```

### More frequent or less frequent ancestry switches

Adjust `--lam-per-Mb`:

- Increase `--lam-per-Mb` when true ancestry switches are expected to be frequent.
- Decrease `--lam-per-Mb` for smoother ancestry tracts.

```bash
python LAI_HMM_v0.4.py ... --lam-per-Mb 0.02
```

### Rebalance haplotype-ID versus variant evidence

When both VCF and `hap_genotype` inputs are supplied, `--b0` shifts the marker-specific mixture:

- Positive `--b0` gives more baseline weight to haplotype-ID evidence.
- Negative `--b0` gives more baseline weight to variant evidence.

```bash
python LAI_HMM_v0.4.py ... --b0 0.5 --no-plots
```

### Tune homozygote tempering

The HMM has two separate mechanisms for preventing homozygous observations from becoming overconfident too quickly. Use `--e-homo` for the baseline rate of expected miscalling of heterozygotes as homozygotes due to technical or biological reasons.
Use the `--hom-*` options to tune further if the model is producing overly confident homozygous calls.

- `--e-homo` handles possible heterozygous undercalling. When the observed genotype is homozygous, the model first mixes in some probability that the true genotype was heterozygous but one allele dropped out or was undercalled. Increase `--e-homo` when low coverage, allele dropout, structural variation, or assay bias makes false homozygotes plausible.
- The internal `soften_mix_for_homozygote()` function allows further, more tune-able correction of homozygous overconfidence without losing information from strongly informative sites. 
It softens homozygous evidence by mixing the original homozygous probability with a neutral probability, using a logistic function of the gap between the most likely and second-most-likely clade assignments to decide how much of the original homozygous evidence to retain.
For example:

- mix = 1.0  -> keep original homozygous probability
- mix = 0.6  -> use 60% original probability + 40% neutral baseline

The mixture is controlled by `--hom-soften-width`, `--hom-soften-delta`, `--hom-min-mix`, and `--hom-neutral`:

- `--hom-min-mix` sets the lower bound on the homozygous evidence mix. Smaller values allow weak homozygous calls to be flattened more strongly toward `--hom-neutral`; setting it to `1` disables homozygote softening entirely.
- `--hom-soften-width` controls how gradually the mix weight increases as the top ancestry probability separates from the second-best probability (more informative). Larger values make the transition broader/more gradual (so clade informativeness has less of an effect on the weighting); smaller values make it sharper/more dependendent on the informativeness.
- `--hom-soften-delta` sets the separation threshold around which the transition is centered. When the difference between the most likely clade pair and the second most likely clade pair = `hom_soften_delta`, the logistic weight is 0.5. 
- `--hom-neutral` sets the neutral target that softened homozygous probabilities are mixed toward.

For example, with the default parameters, the homozygous clade calls probabilities are more strongly flattened towards neutral when the clade evidence is ambiguous:

- ambiguous clade evidence      -> keep ~78–80% original, mix ~20–22% neutral
- moderate separation between most likely and second most likely clade    -> keep ~85–88% original
- strong evidence for most likely clade      -> keep ~90–93% original


## Python API examples

### Build reference files

```python
import importlib.util
from pathlib import Path

module_path = Path("LAI_HMM_v0.4.py")
spec = importlib.util.spec_from_file_location("lai_hmm_v0_4", module_path)
lai = importlib.util.module_from_spec(spec)
assert spec.loader is not None
spec.loader.exec_module(lai)

outputs = lai.build_reference_files(
    outdir="reference",
    clades=["EA", "Mus", "NA1", "NA2", "Vv"],
    membership_path="example_data/example_reference_membership.tsv",
    hap_genotype_path="example_data/hap_genotype_refset_example",
    vcf_path="example_data/example_refset.vcf.gz",
    pca=True,
    verbose=True,
)

print(outputs)
```

### Run the HMM

```python
import importlib.util
from pathlib import Path

module_path = Path("LAI_HMM_v0.4.py")
spec = importlib.util.spec_from_file_location("lai_hmm_v0_4", module_path)
lai = importlib.util.module_from_spec(spec)
assert spec.loader is not None
spec.loader.exec_module(lai)

result = lai.run_lai_hmm(
    vcf_path="example_data/example_samples.vcf.gz",
    hap_genotype_path="example_data/hap_genotype_samples_example",
    marker_positions_path="example_data/marker_positions.csv",
    variant_profiles_path="example_data/reference_variant_profiles.tsv",
    hap_freq_lookup_path="example_data/reference_hap_allele_frequencies.tsv",
    hap_informativeness_path="example_data/reference_hap_allele_informativeness.tsv",
    clades=["EA", "Mus", "NA1", "NA2", "Vv"],
    sample="SampleA",
    outdir="results",
    make_plots=True,
)

print(result["summary_df"])
```

## Troubleshooting

##### Missing reference files

For VCF input, provide a pre-calculated `--variant-profiles` or supply `--reference-membership` to calculate variant profiles from the reference samples in the VCF.

For `hap_genotype` input, provide `--hap-frequencies` or supply `--reference-membership` so the CLI can build haplotype frequencies.

For combined VCF + `hap_genotype` input, also provide `--hap-informativeness` unless it is being built automatically.

##### No VCF loci collected

Likely causes:

- VCF contig names do not match `--contigs`.
- In combined VCF + `hap_genotype` runs, VCF records are missing `INFO/MARKER`.
- The VCF index is missing or incompatible.

Try removing `--contigs`, checking VCF header contig names, and, for combined runs, confirming that `INFO/MARKER` is present and matches the haplotype marker names.

##### Many missing markers

Check marker names across all files:

- `marker_positions.csv`
- VCF `INFO/MARKER` for combined runs or marker-annotated VCF-only runs, or `CHROM:POS` site labels for fully unlabeled VCF-only runs
- `hap_genotype` row names
- reference frequency files

Marker names must match exactly after string conversion.

##### Unexpected clades or state columns

Ensure `--clades` exactly matches the group names and order in the reference frequency files. The HMM creates one state for every unordered diploid pair of clades.

##### Requested sample not found

Check that the sample name appears in at least one genotype input. When both VCF and `hap_genotype` inputs are present, explicitly requested samples may appear in either file, but all-sample combined runs use the intersection by default unless `--sample-source union` is specified.

##### Combined run fails because haplotype informativeness is missing

Provide the informativeness file generated during reference-building:

```bash
--hap-informativeness reference/reference_hap_allele_informativeness.tsv
```

or regenerate it from the frequency lookup:

```bash
python LAI_HMM_v0.4.py \
  --step hap-informativeness \
  --hap-frequencies example_data/reference_hap_allele_frequencies.tsv \
  --clades EA,Mus,NA1,NA2,Vv \
  --outdir reference
```
