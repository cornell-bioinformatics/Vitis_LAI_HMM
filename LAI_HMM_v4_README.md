# LAI HMM v4 CLI

Command-line interface for local ancestry inference (LAI) from VCF variant genotypes and/or `hap_genotype` format marker haplotype allele IDs. Based on allele frequencies in reference populations. Infers ancestry patterns from populations, species, phylogenetic clades, or any other user-specified groups that are expected to differ in allele frequencies. Tunable genotyping error rates for sequencing error and heterozygous undercalling.


## Setup

Use Python 3.9 or newer. A conda environment is recommended because VCF support uses `cyvcf2`.

```bash
conda create -n lai-hmm -c conda-forge -c bioconda python=3.12 pandas numpy matplotlib cyvcf2
conda activate lai-hmm
```

Check the set up by calling the command line interface help menu:

```bash
python LAI_HMM_v4.py --help
```


## Example Files

The `example_data` directory includes small examples for testing:

- `example.vcf.gz` and `example.vcf.gz.csi`
- `hap_genotype_example`
- `marker_positions.csv`
- `example_reference_membership.tsv`


## hap_genotype Format

A `hap_genotype` file is a marker-by-sample matrix of haplotype allele IDs.

- One row per marker.
- First column is the marker name, usually `Locus`, `Marker`, or similar.
- Optional `Haplotypes` column lists observed allele IDs/frequencies and is ignored by the parser.
- Remaining columns are sample names.
- Each sample cell contains two allele IDs, such as `1/2`, `1|2`, or `1/2:read_counts`.
- Missing values can be `./.:0`, `./.`, `.`, `NA`, or blank.

Example:

```text
Locus	Haplotypes	SampleA	SampleB
Marker1	1(0.5);2(0.5);	1/2:10,8	2/2:12
Marker2	3(0.7);4(0.3);	./.:0	3/4:9,2
```

## Reference Membership Format
The HMM requires non-admixed reference samples to estimate the allele frequencies in each population/clade.
Reference-building needs a table mapping reference samples to clades/groups:

```text
SampleA	EA
SampleB	Mus
SampleC	NA
SampleD	Vv
```

Accepted sample column names include `sample`, `IID`, `id`, and `sample_id`. Accepted group column names include `clade`, `group`, `population`, `pop`, `species`, and `index`. If no column names are given, it assumes the first column is the sample name and the second is the group name.

## Build Reference Files
- This step calculates the clade-specific allele frequencies and the allele informativeness for the reference samples in a VCF and/or hap_genotype file.
- Specific and normalized allele "informativeness for assignment" scores are calculated based on Rosenberg et al. 2003 Am J Hum Genet 73(6).
- Listing the clade names is optional; if not specified, it will use the clade names found in the reference-membership list.
- Optional --pca argument runs a principal component analysis of the reference samples and outputs a plot and clade differentiation metrics.
- This step can take some time to run depending on the size of the vcf and/or hap_genotype files.   
```bash
python LAI_HMM_v4.py \
  --step build-reference \
  --hap-genotype hap_genotype_example \
  --vcf example.vcf.gz \
  --reference-membership example_reference_membership.tsv \
  --clades EA,Mus,NA1,NA2,Vv \
  --verbose 
  --pca
  --outdir example_reference
```

Outputs include:

- haplotype allele frequency table
- haplotype allele frequency lookup pickle
- haplotype allele informativeness table
- VCF allele frequency/profile table with `specific` and `normalized_Inf`
- marker/locus informativeness table
- PCA coordinates and PCA plot for reference samples
- Population differentiation metrics, including a centroid-to-within-clade distance ratio

You can build only hap_genotype reference files by omitting `--vcf`, or only VCF reference files by omitting `--hap-genotype`.

## Run The Full Pipeline

If precomputed reference files are missing and `--reference-membership` is supplied, the wrapper builds them automatically under `OUTDIR/reference`.

```bash
python LAI_HMM_v4.py \
  --vcf example.vcf.gz \
  --hap-genotype hap_genotype_example \
  --marker-positions marker_positions.csv \
  --reference-membership example_reference_membership.tsv \
  --clades EA,Mus,NA,Vv \
  --all-samples \
  --outdir results
```


## Run Individual Steps

Build reference files:

```bash
python LAI_HMM_v4.py --step build-reference \
  --hap-genotype hap_genotype_example \
  --reference-membership example_reference_membership.tsv \
  --clades EA,Mus,NA,Vv \
  --outdir reference
```

Calculate variant informativeness from an existing variant profile:

```bash
python LAI_HMM_v4.py --step variant-informativeness \
  --variant-profiles reference_variant_profiles.tsv \
  --clades EA,Mus,NA,Vv \
  --outdir reference_only
```

Calculate haplotype allele-ID informativeness from an existing lookup pickle:

```bash
python LAI_HMM_v4.py --step hap-informativeness \
  --hap-frequencies reference_hap_allele_frequency_lookup.pkl \
  --clades EA,Mus,NA,Vv \
  --outdir reference_only
```

Run only the HMM using precomputed reference files:

```bash
python LAI_HMM_v4.py --step run-hmm \
  --vcf example.vcf.gz \
  --hap-genotype hap_genotype_example \
  --marker-positions marker_positions.csv \
  --variant-profiles reference_variant_profiles.tsv \
  --hap-frequencies reference_hap_allele_frequency_lookup.pkl \
  --hap-informativeness reference_hap_allele_informativeness.tsv \
  --clades EA,Mus,NA,Vv \
  --sample SampleA \
  --outdir results
```

## Tuning Parameters

The wrapper defaults keep optional softening/smoothing off:

- `--alpha 0` and `--windowsize 0`: no context/window smoothing.
- `--e-homo 0`: no homozygous dropout/heterozygote undercalling softening.
- `--hom-min-mix 1`: no homozygote frequency-separation softening.
- `--certainty-softener 1`: no posterior certainty softening.
- `--tau 0`: no posterior-guided Viterbi blending.

Adjust these to your genotyping assay/error model, for example:

```bash
python LAI_HMM_v4.py ... --alpha 0.4 --windowsize 6 --e-homo 0.1 --certainty-softener 3
```

## Main Outputs

The HMM wrapper writes:

- `posteriors/<sample>_<date>.csv`: per-marker calls and posterior probabilities.
- `clade_percentage_summary_<date>.csv`: one row per sample with clade percentages, missingness, and warnings.
- `run_<timestamp>.log`: paths and model parameters.
- `figures/<sample>_<date>.png`: only when `--plots` and `--chrom-lengths` are supplied.

## Troubleshooting

- `ImportError: cyvcf2 is required`: install `cyvcf2` or run hap_genotype-only steps.
- Missing reference files: provide `--variant-profiles` / `--hap-frequencies`, or provide `--reference-membership` so the CLI can build them.
- Many missing markers: check that marker names match across `marker_positions.csv`, VCF `INFO/MARKER`, hap_genotype rows, and reference files.
- Unexpected clades: ensure `--clades` exactly matches the group names and order in the reference frequency files.
