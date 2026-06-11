#!/usr/bin/env python3

"""
LAI_HMM_v0.4.py

Pipeline to infer population/species/clade of origin across segments of the genome from marker data (VCF and/or hap_genotype-format matrix of multiallelic genotypes) 

Likelihoods derived from reference population-specific allele frequencies.

If both input types are given: Per-marker emission probabilities are the product (sum in log-space) of per-variant (from VCF), mixed with haplotype allele ID emission probabilities (), according to their respective informativeness of assignment.

Tunable parameters to account for allele drop out/heterozygous undercalling, genotyping error, and expected recombination frequency.

Assumes diploidy and unphased data.

Inputs
------
1) VCF (INFO field MARKER=<marker_name> is required only when combining with the multiallelic genotype matrix)
2) Reference population-specific variant frequencies
AND/OR
3) Sample multiallelic genotype matrix (aka "hap_genotype" file; sample names in columns and marker names in rows, each cell contains two haplotype allele IDs)
4) Reference population-specific haplotype allele ID frequencies


Notes
-----
- State space is all the possible unordered diploid pairs of k populations/species/clades. E.g. k=4 -> 10 states.
- Transition model and decoding are based on physical distance between markers on the same chromosome.

"""

from dataclasses import dataclass, asdict
from email import parser
from typing import Set, Dict, Tuple, List, Iterable, Optional, Any, Sequence
import pandas as pd
import numpy as np
import math
import sys
import gzip
import re
import json
import time
import pickle
try:
    import pysam
except ImportError:
    pysam = None
from typing import List, Optional, Tuple
import os, subprocess, tempfile, shutil
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed
from itertools import combinations_with_replacement
try:
    from cyvcf2 import VCF
except ImportError:
    VCF = None
#from GenotypeDF_CZ_CGPT import GenotypeDF
from collections import defaultdict

_default_tmp = Path(os.environ.get("TMPDIR", f"/workdir/{os.environ.get('USER', 'ssv42')}/tmp"))
try:
    _default_tmp.mkdir(parents=True, exist_ok=True)
    os.environ.setdefault("TMPDIR", str(_default_tmp))
    os.environ.setdefault("MPLCONFIGDIR", str(_default_tmp / "matplotlib"))
except Exception:
    pass

DEFAULT_CLADES = ("EA", "Mus", "NA1", "NA2", "Vv")
LEGACY_CLADE_COLORS = {
    "EA": "#2c4dff",
    "Mus": "#8b1c62",
    "NA1": "#4daf4a",
    "NA2": "#e1ad01",
    "Vv": "#00cdcd",
}
HMM_DEFAULTS = {
    "lam_per_Mb": 0.05,
    "strict_boost": 5.0,
    "cap_total_boost_per_marker": None,
    "trans_temp": 1.0,
    "alpha": 0.0,
    "windowsize": 0,
    "e_geno": 0.01,
    "e_homo": 0.1,
    "b0": 0.0,
    "certainty_softener": 1.0,
    "hom_soften_delta": 0.1,
    "hom_soften_width": 0.6,
    "hom_min_mix": 0.6,
    "hom_neutral": 0.5,
    "tau": 1.0,
}
CHR20_START_BP = 21_600_000

# ---------------------------
# Data loading functions
# ---------------------------

# Open gzipped files
def _open_text_maybe_gzip(path: str):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "rt")

def chrom_to_numeric(chrom):
    """
    Convert chromosome strings such as:
        chr01 -> 1
        chr1 -> 1
        Chromosome12 -> 12
        chromosome_5 -> 5
        CHR07 -> 7
    to numeric values.
    """
    if pd.isna(chrom):
        return pd.NA
    chrom = str(chrom).strip()
    # Remove common chromosome name prefixes
    chrom = re.sub(
        r"^(chromosome|chrom|chr)[_\-\s]*",
        "", chrom, flags=re.IGNORECASE)

    # Extract first integer
    match = re.search(r"\d+", chrom)
    if match:
        return int(match.group())
    return pd.NA

# VCF parsing with cyvcf2
@dataclass
class VCFVariant:
    chrom: str
    pos: int
    id: str
    ref: str
    alts: List[str]
    marker: str
    gt_alleles: Tuple[str, str]  # concrete allele strings (REF/ALT), unphased


@dataclass
class VCFChunk:
    samples: List[str]
    sample_to_idx: Dict[str, int]
    locus_df: pd.DataFrame               # columns: chrom,pos,id,marker,ref,alts (list[str])
    G1: np.ndarray                       # (n_loci, N) int8: -1=missing, 0=REF, 1..k=ALT index
    G2: np.ndarray                       # (n_loci, N) same

    def get_sample_variants(self, sample: str) -> List["VCFVariant"]:
        j = self.sample_to_idx[sample]
        g1 = self.G1[:, j]
        g2 = self.G2[:, j]
        mask = ~((g1 < 0) & (g2 < 0))
        if not np.any(mask):
            return []
        sub = self.locus_df.loc[mask]
        g1s, g2s = g1[mask], g2[mask]
        chroms = sub["chrom"].to_numpy()
        poss   = sub["pos"].to_numpy()
        ids    = sub["id"].to_numpy()
        markers= sub["marker"].to_numpy()
        refs   = sub["ref"].to_numpy()
        alts   = sub["alts"].to_list()
        out: List[VCFVariant] = []
        for i in range(len(sub)):
            ref = refs[i]; alt_list = alts[i]
            a1 = None if g1s[i] < 0 else (ref if g1s[i] == 0 else (alt_list[g1s[i]-1] if g1s[i]-1 < len(alt_list) else None))
            a2 = None if g2s[i] < 0 else (ref if g2s[i] == 0 else (alt_list[g2s[i]-1] if g2s[i]-1 < len(alt_list) else None))
            out.append(VCFVariant(
                chrom=str(chroms[i]),
                pos=int(poss[i]),
                id=str(ids[i]),
                ref=str(ref),
                alts=[str(a) for a in alt_list],
                marker=str(markers[i]),
                gt_alleles=(a1, a2),
            ))
        return out


def load_vcf_chunk_matrix(
    vcf_path: str,
    contig_whitelist: Optional[Set[str]] = None,   # optionally select certain contigs
    sample_names: Optional[Sequence[str]] = None,  # optionally select sample columns
    require_marker: bool = True,                   # enforce INFO/MARKER presence
    debug: bool = False) -> VCFChunk:
    """
    Load an already chunked, indexed VCF/BCF (e.g., ~200 samples) and build:
      - locus_df: per-locus metadata (shared by all samples)
      - G1, G2: int8 allele-index matrices (n_loci x n_samples)
    Multiallelic-safe via allele indices (0=REF, 1..k=ALT idx, -1=missing).
    """
    if VCF is None:
        raise ImportError("cyvcf2 is required for VCF input. Install it with `pip install cyvcf2`.")

    if sample_names is None:
        vcf = VCF(vcf_path, gts012=True)
    else:
        vcf = VCF(vcf_path, gts012=True, samples=list(sample_names))
    samples = list(vcf.samples)
    if not samples:
        vcf.close(); raise ValueError("VCF has no samples.")
    N = len(samples)

    # Pick contigs present in header
    header_contigs = set(vcf.seqnames)
    if contig_whitelist is not None:
        contigs = [c for c in header_contigs if c in contig_whitelist]
    else:
        contigs = header_contigs
    if not contigs:
        vcf.close()
        raise ValueError(f"VCF header contigs empty or do not match expectation")
    use_marker_annotations = require_marker or vcf_has_any_marker_annotations(vcf, contigs)

    # -------- Pass 1: enumerate loci and build row maps (robust to ordering) --------
    loci: List[Tuple[str,int,str,str,str,List[str]]] = []
    # per-contig: map (pos, id_or_marker) -> row index
    rowmap: Dict[str, Dict[Tuple[int, str], int]] = {c: {} for c in contigs}

    unnamed_marker_count = 0
    #errors = []
    for chrom in contigs:
        for rec in vcf(f"{chrom}"):
            marker = resolve_vcf_marker_name(
                rec,
                require_marker=require_marker,
                use_marker_annotations=use_marker_annotations,
            )

            if marker is None:
                unnamed_marker_count = unnamed_marker_count + 1
                #errors.append(f"{rec.CHROM}:{rec.POS}")
                continue
            rid = rec.ID if rec.ID not in (None, ".", "") else f"{rec.CHROM}:{rec.POS}"
            alt_list = [str(a).upper() for a in (rec.ALT or [])]
            row = len(loci)
            loci.append((
                str(rec.CHROM),
                int(rec.POS),
                str(rid),
                str(marker),
                str(rec.REF).upper(),
                alt_list
            ))
            # row key prefers ID if present, else marker
            key = (int(rec.POS), str(rid) if rec.ID not in (None, ".", "") else str(marker))
            rowmap[chrom][key] = row

    if debug == True:
        print(unnamed_marker_count, "variants without marker name annotation")
    
    if not loci:
        vcf.close()
        raise ValueError(
            "No loci collected (check contigs and VCF content; combined VCF + hap_genotype "
            "runs also require INFO/MARKER)."
        )

    locus_df = pd.DataFrame(loci, columns=["chrom","pos","id","marker","ref","alts"])
    n_loci = len(locus_df)

    # -------- Pass 2: allocate and fill matrices --------
    G1 = np.full((n_loci, N), -1, dtype=np.int8)
    G2 = np.full((n_loci, N), -1, dtype=np.int8)

    def enc(ix: int, n_alts: int) -> int:
        if ix is None or ix < 0: return -1
        if ix == 0: return 0
        return ix if ix <= n_alts else -1

    # Re-iterate records; map to row via (pos, id_or_marker)
    for chrom in contigs:
        rmap = rowmap[chrom]
        for rec in vcf(f"{chrom}"):
            marker = resolve_vcf_marker_name(
                rec,
                require_marker=require_marker,
                use_marker_annotations=use_marker_annotations,
            )
            if marker is None:
                continue
            rid_present = rec.ID not in (None, ".", "")
            key = (int(rec.POS), str(rec.ID) if rid_present else str(marker))
            row = rmap.get(key)
            if row is None:
                continue  # shouldn’t happen; safety for edge cases

            n_alts = len(rec.ALT or [])
            gts = rec.genotypes
            if not gts:
                continue
            for j in range(N):
                g = gts[j]
                if len(g) >= 2:
                    G1[row, j] = enc(g[0], n_alts)
                    G2[row, j] = enc(g[1], n_alts)

    vcf.close()

    sample_to_idx = {s: i for i, s in enumerate(samples)}
    return VCFChunk(samples=samples,
                    sample_to_idx=sample_to_idx,
                    locus_df=locus_df,
                    G1=G1, G2=G2)


# ---------------------------
# Vitis-specific functions
# ---------------------------

# Check for chromosome 20 presence/absence:

# Count how many muscadine chr 20 specific alleles are present
def check_for_chr_20(
    mus_hap_alleles: pd.DataFrame,
    NONMus_hap_alleles: Optional[pd.DataFrame],
    haplotypes: pd.Series,
):
    mus_matches = []
    
    for marker, group in mus_hap_alleles.groupby("Marker"):
        alleles = group["Allele"].astype(str).tolist()
        hap = haplotypes.get(marker)
    
        if not isinstance(hap, list):
            num_matching = 0
        else:
            num_matching = sum(a in hap for a in alleles)
    
        mus_matches.append({"Marker": marker, "num_matching": num_matching})
    
    mus_matches = pd.DataFrame(mus_matches)
    
    if sum(mus_matches["num_matching"]) > 6:
        chr20presence = True
    
        if mus_matches["num_matching"].value_counts().get(2, 0) > 2:
            homozygous = True
        else: 
            homozygous = False
            
    else: 
        chr20presence = False; homozygous = False

    nonmus_check_available = (
        NONMus_hap_alleles is not None
        and {"Marker", "Allele"}.issubset(NONMus_hap_alleles.columns)
    )
    if chr20presence and not homozygous and nonmus_check_available: #check for chr. 7/20 hemizygosity
        Euv_matches = []

        for marker, group in NONMus_hap_alleles.groupby("Marker"):
            alleles = group["Allele"].astype(str).tolist()
            hap = haplotypes.get(marker)
            if not isinstance(hap, list):
                num_matching = 0
            else:
                num_matching = sum(a in hap for a in alleles)
        
            Euv_matches.append({"Marker": marker, "num_matching": num_matching})
        Euv_matches = pd.DataFrame(Euv_matches)
        if sum(Euv_matches["num_matching"]) > 2:
            #print("possible chr.7/chr.20 hemizygosity")
            homozygous = "hemizygous"
    
    return(chr20presence, homozygous, mus_matches)


# --------------------------------
# Helpers for combining emissions
# --------------------------------


def _row_norm(logB, floor=-60.0):
    """Row max-center and floor; returns a new array."""
    X = np.array(logB, copy=True)
    X[~np.isfinite(X)] = floor
    m = np.max(X, axis=1, keepdims=True)
    X = np.clip(X - m, floor, 0.0)
    return X


def hap_known_from_reference(test_series, freq_lookup, marker_order):
    """
    test_series: pd.Series indexed by Marker -> (alleleA, alleleB) strings
    freq_lookup: dict[(Marker, Allele)] -> np.array([EA, Mus, NA1, NA2, Vv])
    Returns: boolean array (T,) aligned to marker_order
    """
    known = np.zeros(len(marker_order), dtype=bool)
    for t, m in enumerate(marker_order):
        alleles = test_series.get(m, None)
        if isinstance(alleles, (list, tuple)) and len(alleles) == 2:
            a, b = str(alleles[0]), str(alleles[1])
            # mark known if either allele has a reference profile
            known[t] = ((m, a) in freq_lookup) or ((m, b) in freq_lookup)
    return known



def _fuse_noisy_or(values: Iterable[float]) -> float:
    """Bounded [0,1] combine: 1 - Π(1 - v) for diploid informativeness score fuse"""
    x = 1.0
    for v in values:
        v = float(np.clip(v, 0.0, 1.0))
        x = x * (1.0 - v)
    return 1.0 - x
    

def _normalize_hap_value(v: Any) -> Tuple[str, str]:
    """
    Normalize a haplotype field into a pair of allele IDs or (None, None).
    Accepts: ('H1','H2'), ['H1','H2'], 'H1|H2', 'H1/H2',
    'H1/H2:read_counts', 'H1', None, np.nan, floats.
    Returns strings (or None) without extra whitespace.
    """
    if v is None:
        return (None, None)
    # Handle NaN-ish
    try:
        if isinstance(v, float) and np.isnan(v):
            return (None, None)
    except Exception:
        pass

    if isinstance(v, (tuple, list)):
        if len(v) >= 2:
            a1, a2 = v[0], v[1]
        elif len(v) == 1:
            a1, a2 = v[0], None
        else:
            return (None, None)
    elif isinstance(v, str):
        s = v.strip()
        if s in ("", ".", "./.", ".|.", "NA", "NaN", "nan", "None", "none"):
            return (None, None)
        # hap_genotype cells commonly store read/count details after the allele pair.
        if ":" in s:
            s = s.split(":", 1)[0].strip()
        if "|" in s:
            a1, a2 = s.split("|", 1)
        elif "/" in s:
            a1, a2 = s.split("/", 1)
        else:
            a1, a2 = s, None
    else:
        # last resort: try tuple-cast
        try:
            a1, a2 = tuple(v)
            if len((a1, a2)) != 2:
                return (None, None)
        except Exception:
            return (None, None)

    missing = {"", ".", "./.", ".|.", "NA", "NaN", "nan", "None", "none"}
    a1 = None if a1 is None or str(a1).strip() in missing else str(a1).strip()
    a2 = None if a2 is None or str(a2).strip() in missing else str(a2).strip()
    return (a1, a2)


def hap_strength_metrics_sample(
    haplotypes: pd.Series,
    hap_inf: pd.DataFrame,
    marker_order: List[str],
    hap_cols=("marker","allele_id","specific","normalized_Inf"),
    clades: Optional[Sequence[str]] = None,
) -> np.ndarray:
    """
    Returns per-marker array length T of sample-specific hap-ID informativeness in [0,1],
    using hap_inf['normalized_Inf'] for each of the individual's two haplotype allele IDs
    and fusing them with noisy-OR.
    """
    # Ensure required columns & normalized_Inf
    if not {"marker","allele_id"}.issubset(hap_inf.columns):
        raise ValueError("hap_inf must have columns: 'marker' and 'allele_id'")
    if "normalized_Inf" not in hap_inf.columns:
        if "specific" in hap_inf.columns:
            hap_inf = hap_inf.copy()
            K = len(clades) if clades is not None else _infer_k_from_freq_columns(hap_inf)
            hap_inf["normalized_Inf"] = hap_inf["specific"] / np.log(float(K))
        else:
            raise ValueError("hap_inf needs either 'normalized_Inf' or 'specific'.")

    # Standardize key types to str
    key_df = hap_inf.copy()
    key_df["marker"] = key_df["marker"].astype(str)
    key_df["allele_id"] = key_df["allele_id"].astype(str)

    hap_lookup = pd.Series(
        key_df["normalized_Inf"].to_numpy(),
        index=pd.MultiIndex.from_frame(key_df[["marker","allele_id"]])
    ).to_dict()

    # If haplotypes index aren't strings, cast for lookup consistency
    # Expect haplotypes to be a Series mapping marker -> (allele1, allele2 or string)
    out = np.zeros(len(marker_order), dtype=float)
    for t, m in enumerate(marker_order):
        mk = str(m)
        a1, a2 = _normalize_hap_value(haplotypes.get(mk, haplotypes.get(m, None)))
        v1 = hap_lookup.get((mk, a1), 0.0) if a1 else 0.0
        v2 = hap_lookup.get((mk, a2), 0.0) if a2 else 0.0
        out[t] = _fuse_noisy_or([v1, v2])
    return out


def _variant_lookup_from_prof(prof_df: pd.DataFrame) -> Dict[Tuple[str,int,str], float]:
    """
    Build {(CHROM, POS, ALLELE) -> normalized_Inf} dict from prof_df.
    Accepts either 'specific' (nats) or 'normalized_Inf' (0..1). If only 'specific' is present,
    we compute normalized_Inf = specific / ln(5).
    """
    need = {"CHROM","POS","ALLELE"}
    if not need.issubset(prof_df.columns):
        raise ValueError("prof_df must have columns CHROM, POS, ALLELE")

    if "normalized_Inf" in prof_df.columns:
        col = "normalized_Inf"
        work = prof_df.copy()
    elif "specific" in prof_df.columns:
        work = prof_df.copy()
        K = _infer_k_from_freq_columns(work)
        work["normalized_Inf"] = work["specific"] / np.log(float(K))
        col = "normalized_Inf"
    else:
        raise ValueError("prof_df needs 'normalized_Inf' or 'specific'.")

    # standardize key types
    work["CHROM"]  = work["CHROM"].astype(str)
    work["POS"]    = work["POS"].astype(int)
    work["ALLELE"] = work["ALLELE"].astype(str)

    # build the lookup
    lut = (
        work
        .set_index(["CHROM","POS","ALLELE"])[col]
        .astype(float)
        .to_dict()
    )
    return lut
    

def sum_In_per_marker_sample(
    variants: List["VCFVariant"],
    marker_order: List[str],
    prof_df: pd.DataFrame) -> np.ndarray:
    """
    Returns per-marker array length T of sample-specific variant informativeness
    in [0,1], by:
      (a) fusing the two alleles at each variant site (noisy-OR),
      (b) fusing across all sites within the marker (noisy-OR).
    """
    varLUT = _variant_lookup_from_prof(prof_df)

    # group sample variants by marker
    by_marker: Dict[str, List["VCFVariant"]] = {}
    for v in variants:
        key = str(v.marker).strip() if v.marker is not None else ""
        if key:
            by_marker.setdefault(key, []).append(v)

    out = np.zeros(len(marker_order), dtype=float)
    for t, m in enumerate(marker_order):
        site_scores = []
        for v in by_marker.get(m, []):
            chrom = str(v.chrom)
            pos   = int(v.pos)
            # gt_alleles is a tuple like ('A','A') or ('A','G')
            a1, a2 = v.gt_alleles if isinstance(v.gt_alleles, (tuple,list)) and len(v.gt_alleles)==2 else (None, None)
            s1 = varLUT.get((chrom, pos, str(a1)), 0.0) if a1 else 0.0
            s2 = varLUT.get((chrom, pos, str(a2)), 0.0) if a2 else 0.0
            site_scores.append(_fuse_noisy_or([s1, s2]))  # per-site diploid fuse
        out[t] = _fuse_noisy_or(site_scores) if site_scores else 0.0
    return out



def hap_strength_metrics(logB_fromhaps, floor=-60.0):
    """
    Input: logB_fromhaps (T, S)
    Returns per-marker arrays: ratio>=1, margin in [0,1], ent_norm in [0,1], hap_known bool
    """
    X = np.array(logB_fromhaps, copy=True)
    X[~np.isfinite(X)] = floor

    # Row max-center for numerical stability
    X = X - X.max(axis=1, keepdims=True)
    X = np.clip(X, floor, 0.0)

    # Probabilities
    P = np.exp(X)
    P /= P.sum(axis=1, keepdims=True)

    # Top1 / Top2
    top2 = np.partition(P, -2, axis=1)[:, -2]
    top1 = P.max(axis=1)
    ratio = top1 / np.maximum(top2, 1e-12)

    # Margin certainty (0..1)
    S = P.shape[1]
    margin = np.clip((top1 - 1.0/S) / (1.0 - 1.0/S), 0.0, 1.0)

    # Normalized entropy (0..1; 1 ~ uniform)
    ent = -(P * np.log(np.clip(P, 1e-300, 1))).sum(axis=1)
    ent_norm = ent / np.log(S)

    # Is the hap-ID emission informative? (non-uniform row)
    nearly_uniform = ent_norm > 0.98
    literal_uniform = np.all(np.isclose(P, 1.0/S, atol=1e-9), axis=1)
    hap_known = ~(nearly_uniform | literal_uniform)

    return ratio, margin, ent_norm, hap_known


def sum_In_per_marker(variants: List["VCFVariant"],
                      marker_order: List[str],
                      variant_inf: pd.Series) -> np.ndarray:
    """
    Sum per-site informativeness over variants belonging to each marker.
    variant_inf_index: Series indexed by (CHROM, POS) -> I_n (float).
    Returns array length T (len(marker_order)).
    """
    # group variants by marker once
    by_marker: Dict[str, List["VCFVariant"]] = {}
    for v in variants:
        key = str(v.marker).strip() if v.marker is not None else ""
        if key:
            by_marker.setdefault(key, []).append(v)

    T = len(marker_order)
    out = np.zeros(T, dtype=float)
    for t, m in enumerate(marker_order):
        total = 0.0
        for v in by_marker.get(m, []):
            key = (str(v.chrom), int(v.pos))
            try:
                total += float(variant_inf.loc[key])
            except KeyError:
                pass
        out[t] = total
    return out



# ---------------------------
# Variant frequency profiles
# ---------------------------

def _infer_k_from_freq_columns(df: pd.DataFrame) -> int:
    """Infer number of clade frequency columns from a profile-like table."""
    reserved = {
        "chrom", "pos", "id", "marker", "ref", "alts", "allele_index", "allele",
        "specific", "normalized_inf", "informativeness", "marker_informativeness",
    }
    clade_cols = [c for c in df.columns if str(c).strip().lower() not in reserved]
    numeric_cols = []
    for c in clade_cols:
        vals = pd.to_numeric(df[c], errors="coerce")
        if vals.notna().any():
            numeric_cols.append(c)
    if not numeric_cols:
        raise ValueError("Could not infer clade frequency columns.")
    return len(numeric_cols)


def _validate_clades(clades: Sequence[str]) -> List[str]:
    out = [str(c).strip() for c in clades if str(c).strip()]
    if len(out) < 2:
        raise ValueError("Provide at least two clade/group names.")
    if len(set(out)) != len(out):
        raise ValueError(f"Clade/group names must be unique: {out}")
    return out


def _use_legacy_clade_colors(clades: Sequence[str]) -> bool:
    clade_list = [str(c).strip() for c in clades if str(c).strip()]
    return len(clade_list) == len(DEFAULT_CLADES) and set(clade_list) == set(DEFAULT_CLADES)


def _colorblind_palette(n_colors: int) -> List[str]:
    if n_colors <= 0:
        return []

    try:
        import matplotlib.pyplot as plt

        prop_cycle = plt.style.library.get("tableau-colorblind10", {}).get("axes.prop_cycle")
        palette = list(prop_cycle.by_key().get("color", [])) if prop_cycle is not None else []
    except Exception:
        palette = []

    if not palette:
        palette = [
            "#006BA4", "#FF800E", "#ABABAB", "#595959", "#5F9ED1",
            "#C85200", "#898989", "#A2C8EC", "#FFBC79", "#CFCFCF",
        ]

    if n_colors <= len(palette):
        return palette[:n_colors]

    try:
        from matplotlib import cm, colors as mcolors

        extra = [mcolors.to_hex(cm.get_cmap("tab20", n_colors)(i)) for i in range(n_colors)]
    except Exception:
        extra = []

    merged: List[str] = []
    for color in palette + extra:
        if color not in merged:
            merged.append(color)
        if len(merged) >= n_colors:
            return merged[:n_colors]

    repeats = (n_colors + len(palette) - 1) // len(palette)
    return (palette * repeats)[:n_colors]


def get_clade_colors(clades: Sequence[str]) -> Dict[str, str]:
    clade_list = [str(c).strip() for c in clades if str(c).strip()]
    if not clade_list:
        return {}
    if _use_legacy_clade_colors(clade_list):
        return {clade: LEGACY_CLADE_COLORS[clade] for clade in clade_list}
    return dict(zip(clade_list, _colorblind_palette(len(clade_list))))


def allele_informativeness_by_locus(
    cladefreqs: pd.DataFrame,
    clades: Sequence[str],
    by_cols: Sequence[str] = ("CHROM", "POS"),
    eps: float = 1e-12,
) -> pd.DataFrame:
    """
    Rosenberg-style locus informativeness for assignment, in natural-log units.
    """
    clades = _validate_clades(clades)
    missing = [c for c in clades if c not in cladefreqs.columns]
    if missing:
        raise ValueError(f"Frequency table is missing clade columns: {missing}")
    df = cladefreqs.copy()
    for c in clades:
        df[c] = pd.to_numeric(df[c], errors="coerce").fillna(0.0)
    weights = np.ones(len(clades), dtype=float) / len(clades)

    rows = []
    for key, locus in df.groupby(list(by_cols), dropna=False):
        freqs = locus[list(clades)].to_numpy(dtype=float)
        freqs = np.clip(freqs, eps, 1.0)
        p_bar = freqs @ weights
        h_c = -np.sum(freqs * np.log(freqs), axis=0)
        h_bar = -np.sum(p_bar * np.log(np.clip(p_bar, eps, 1.0)))
        item = {col: val for col, val in zip(by_cols, key if isinstance(key, tuple) else (key,))}
        item["informativeness"] = float(h_bar - np.dot(weights, h_c))
        item["normalized_Inf"] = float(item["informativeness"] / np.log(float(len(clades))))
        rows.append(item)
    return pd.DataFrame(rows)


def variant_allele_informativeness(
    cladefreqs: pd.DataFrame,
    clades: Optional[Sequence[str]],
    by_cols: Sequence[str] = ("CHROM", "POS"),
    eps: float = 1e-12,
) -> pd.DataFrame:
    """
    Append per-variant-allele informativeness columns:
    specific = KL(p(clade | allele) || uniform), normalized_Inf = specific / ln(K).
    Clade names are either provided or inferred from numeric columns in the table. 
    Requires allele-specific frequencies and a column "ALLELE" with the allele string.
    """
    if clades:
        clades = _validate_clades(clades)
        missing = [c for c in clades if c not in cladefreqs.columns]
        if missing:
            raise ValueError(f"Frequency table is missing clade columns: {missing}")
    else:
        clades = [c for c in cladefreqs.columns if c not in ["CHROM", "POS", "ID", "MARKER", "REF", "ALTS", "ALLELE_INDEX", "ALLELE", "specific", "n"] and pd.to_numeric(cladefreqs[c], errors="coerce").notna().any()]
        print(f"Inferred clade columns: {clades}")
    df = cladefreqs.copy()
    for c in clades:
        df[c] = pd.to_numeric(df[c], errors="coerce").fillna(0.0)
    K = len(clades)
    w = 1.0 / K
    log_w = np.log(w)
    lnK = np.log(float(K))
    col_sums = df.groupby(list(by_cols))[list(clades)].transform("sum").replace(0.0, np.nan)
    P_ac = df[list(clades)].clip(lower=0.0).div(col_sums).fillna(0.0)
    p_a = P_ac.mean(axis=1).clip(lower=eps)
    post = P_ac.div(p_a, axis=0) * w
    post_sum = post.sum(axis=1).replace(0.0, np.nan)
    post = post.div(post_sum, axis=0).fillna(0.0).clip(lower=eps)
    post = post.div(post.sum(axis=1), axis=0)
    specific = (post * (np.log(post) - log_w)).sum(axis=1)
    df["specific"] = specific.astype(float).values
    df["normalized_Inf"] = (specific / lnK).astype(float).values
    return df


def hap_allele_informativeness(
    allele_freq_lookup: Dict[Tuple[str, str], np.ndarray],
    clades: Sequence[str],
    eps: float = 1e-12,
) -> pd.DataFrame:
    """
    Compute per-hap_genotype allele-ID informativeness from a frequency lookup.
    """
    clades = _validate_clades(clades)
    
    K = len(clades)
    w = np.full(K, 1.0 / K, dtype=float)
    by_marker: Dict[str, List[str]] = {}
    for (m, a), vec in allele_freq_lookup.items():
        vec = np.asarray(vec, dtype=float)
        if vec.shape[0] != K:
            raise ValueError(f"{(m, a)} has vector length {vec.shape[0]} but expected {K}.")
        by_marker.setdefault(str(m), []).append(str(a))

    rows = []
    for m, alleles in by_marker.items():
        alleles = sorted(set(alleles))
        P_ac = np.zeros((len(alleles), K), dtype=float)
        for i, a in enumerate(alleles):
            P_ac[i] = np.maximum(np.asarray(allele_freq_lookup[(m, a)], dtype=float), 0.0)
        col_sums = P_ac.sum(axis=0, keepdims=True)
        col_sums = np.where(col_sums <= 0, 1.0, col_sums)
        P_ac = P_ac / col_sums
        p_a = np.clip(P_ac @ w, eps, 1.0)
        post = (P_ac * w) / p_a[:, None]
        post = np.clip(post, eps, 1.0)
        post = post / post.sum(axis=1, keepdims=True)
        specific = np.sum(post * (np.log(post) - np.log(w)[None, :]), axis=1)
        for i, a in enumerate(alleles):
            rows.append({
                "marker": m,
                "allele_id": a,
                "specific": float(specific[i]),
                "normalized_Inf": float(specific[i] / np.log(float(K))),
            })
    return pd.DataFrame.from_records(rows)


def load_variant_profiles(tsv_path: str, clades: Optional[Sequence[str]] = None):
    """
    Returns:
      df: clade-specific allele frequencies DataFrame
      freq_map: dict[(chrom,pos,allele)] -> np.array([EA, Mus, NA1, NA2, Vv])
    """
    df = pd.read_csv(tsv_path, sep="\t", dtype=str)
    clades = _validate_clades(clades or DEFAULT_CLADES)
    # coerce types
    df["CHROM"] = df["CHROM"].astype(str)
    df["POS"] = pd.to_numeric(df["POS"], errors="coerce").astype("Int64")
    missing = [c for c in clades if c not in df.columns]
    if missing:
        raise ValueError(f"Variant profile file is missing clade columns: {missing}")
    for c in clades:
        df[c] = pd.to_numeric(df[c], errors="coerce")

    # normalize text fields
    df["ALLELE"] = df["ALLELE"].astype(str).str.upper()
    df["MARKER"] = df["MARKER"].astype(str)
    df = df.dropna(subset=["POS"])

    # build frequency map
    freq_map: Dict[Tuple[str,int,str], np.ndarray] = {}
    for _, row in df.iterrows():
        chrom = str(row["CHROM"])
        pos = int(row["POS"])
        allele = str(row["ALLELE"]).upper()
        vec = np.array([row[c] for c in clades], dtype=float)
        # clip for numeric stability
        eps = 1e-12
        vec = np.clip(vec, eps, 1.0 - eps)
        freq_map[(chrom, pos, allele)] = vec
    return df, freq_map

# ---------------------------
# Marker ordering / positions
# ---------------------------

def marker_df_from_chunk(chunk: VCFChunk) -> pd.DataFrame:
    df = chunk.locus_df
    # ensure numeric chrom for ordering; keep string too
    df["Chrom_numeric"] = df["chrom"].astype(str).map(chrom_to_numeric)
    df = df.dropna(subset=["Chrom_numeric"])
    # Collapse to one row per marker using median pos
    marker_df = (df.groupby("marker", as_index=False)
                   .agg({"chrom":"first",
                         "pos":lambda x: int(np.median(x)),
                         "Chrom_numeric":"first"})
                   .sort_values(["Chrom_numeric","pos"])
                   .reset_index(drop=True))
    #check uniqueness
    assert marker_df["marker"].is_unique
    return marker_df


def build_variant_pos_df(variants: List[VCFVariant]) -> pd.DataFrame:
    """
    Build a per-marker position table from parsed variants.
    Ensures columns: Marker, Chrom, Pos, Chrom_numeric.
    """
    import re
    rows = []

    for v in variants:
        marker = v.marker
        chrom = v.chrom
        pos = v.pos
        cn = chrom_to_numeric(chrom)
        rows.append({
            "Marker": marker,
            "Chrom": str(chrom),
            "Pos": int(pos),
            "Chrom_numeric": cn
        })

    # Always create the expected columns, even if rows is empty
    variant_df = pd.DataFrame(rows, columns=["Marker", "Chrom", "Pos", "Chrom_numeric"])

    # Drop rows with unknown chromosomes
    bad = variant_df["Chrom_numeric"].isna()
    if bad.any():
        variant_df = variant_df.loc[~bad]

    if variant_df.empty:
        raise ValueError("No usable markers found to build positions.")

    variant_df = variant_df.sort_values(["Chrom_numeric", "Pos"]).reset_index(drop=True)

    # Collapse to one row per marker using the median position
    marker_df = (variant_df.groupby("Marker", as_index=False)
    .agg({"Chrom": "first",
        "Pos": lambda x: round(x.median()),
        "Chrom_numeric": "first"})
    .sort_values(["Chrom_numeric", "Pos"]).reset_index(drop=True) )

    return variant_df, marker_df


def _guess_table_sep(path: str, sep: Optional[str] = None) -> str:
    if sep is not None:
        return sep
    lower = str(path).lower()
    if lower.endswith((".csv", ".csv.gz")):
        return ","
    return "\t"


def read_table(path: str, sep: Optional[str] = None, **kwargs) -> pd.DataFrame:
    """Read a delimited text table, defaulting to comma for .csv and tab otherwise."""
    return pd.read_csv(path, sep=_guess_table_sep(path, sep), **kwargs)


def load_pickle(path: str) -> Any:
    """Load a pickle file used for haplotype frequency or strict diagnostic dictionaries."""
    with open(path, "rb") as f:
        return pickle.load(f)


def normalize_vcf_marker_value(marker: Any) -> Optional[str]:
    """Normalize a VCF INFO/MARKER value to a clean string or None."""
    if isinstance(marker, (list, tuple)):
        marker = marker[0] if marker else None
    if marker is None:
        return None
    marker = str(marker).strip()
    return None if marker in {"", ".", "None"} else marker


def variant_site_label(chrom: Any, pos: Any) -> str:
    """Fallback site label for unlabeled VCF-only records."""
    return f"{str(chrom)}:{int(pos)}"


def vcf_has_any_marker_annotations(vcf: Any, contigs: Sequence[str]) -> bool:
    """Return True if any selected VCF record has a usable INFO/MARKER value."""
    for chrom in contigs:
        for rec in vcf(f"{chrom}"):
            if normalize_vcf_marker_value(rec.INFO.get("MARKER")) is not None:
                return True
    return False


def resolve_vcf_marker_name(
    rec: Any,
    require_marker: bool = True,
    use_marker_annotations: bool = False,
) -> Optional[str]:
    """
    Resolve the marker label for a VCF record.

    If marker annotations are in use for this VCF, unlabeled records are omitted.
    Otherwise, unlabeled VCF-only records fall back to CHROM:POS so each variant
    position is treated as its own ancestry site.
    """
    marker = normalize_vcf_marker_value(rec.INFO.get("MARKER"))
    if marker is not None:
        return marker
    if require_marker or use_marker_annotations:
        return None
    return variant_site_label(rec.CHROM, rec.POS)


def resolve_require_marker(
    require_marker: Optional[bool],
    vcf_path: Optional[str],
    hap_genotype_path: Optional[str],
) -> bool:
    """
    Resolve whether VCF INFO/MARKER annotations are required.

    Auto mode (require_marker=None) requires MARKER only when VCF records need to
    align to hap_genotype marker names. Explicit False is ignored in combined mode
    because the HMM cannot merge VCF and haplotype evidence without marker labels.
    """
    if not vcf_path:
        return False
    if hap_genotype_path:
        return True
    if require_marker is None:
        return False
    return bool(require_marker)


def load_hap_frequency_lookup(
    path: str,
    clades: Optional[Sequence[str]] = None,
    sep: Optional[str] = None,
) -> Dict[Tuple[str, str], np.ndarray]:
    """
    Load haplotype allele-frequency lookup data from either a pickle or a text table.

    Text tables must contain marker/locus and allele-id columns, plus one column per
    clade with numeric frequencies. The long-form TSV written by build-reference is
    accepted directly.
    """
    lower = str(path).lower()
    if lower.endswith((".pkl", ".pickle")):
        data = load_pickle(path)
        if not isinstance(data, dict):
            raise ValueError(f"Expected a dictionary in haplotype frequency pickle: {path}")
        return data

    df = read_table(path, sep=sep, keep_default_na=False)
    lower_to_col = {str(c).strip().lower(): c for c in df.columns}
    marker_aliases = ("marker", "markers", "locus", "loci", "id")
    allele_aliases = (
        "allele_id",
        "alleleid",
        "allele",
        "hap_allele",
        "haplotype_allele",
        "haplotype_allele_id",
    )

    marker_col = next((lower_to_col[a] for a in marker_aliases if a in lower_to_col), None)
    allele_col = next((lower_to_col[a] for a in allele_aliases if a in lower_to_col), None)
    if (marker_col is None or allele_col is None) and sep is None:
        alt_sep = "," if _guess_table_sep(path, sep) == "\t" else "\t"
        alt_df = pd.read_csv(path, sep=alt_sep, keep_default_na=False)
        alt_lower_to_col = {str(c).strip().lower(): c for c in alt_df.columns}
        alt_marker_col = next((alt_lower_to_col[a] for a in marker_aliases if a in alt_lower_to_col), None)
        alt_allele_col = next((alt_lower_to_col[a] for a in allele_aliases if a in alt_lower_to_col), None)
        if alt_marker_col is not None and alt_allele_col is not None:
            df = alt_df
            lower_to_col = alt_lower_to_col
            marker_col = alt_marker_col
            allele_col = alt_allele_col
    if marker_col is None or allele_col is None:
        raise ValueError(
            "Haplotype frequency table must include marker/locus and allele_id/allele columns."
        )

    df = df.rename(columns={marker_col: "marker", allele_col: "allele_id"}).copy()
    lower_to_col = {str(c).strip().lower(): c for c in df.columns}

    if clades:
        clades = _validate_clades(clades)
        missing = [c for c in clades if c not in df.columns and c.lower() not in lower_to_col]
        if missing:
            raise ValueError(f"Haplotype frequency table is missing clade columns: {missing}")
        rename = {}
        for c in clades:
            if c not in df.columns:
                rename[lower_to_col[c.lower()]] = c
        if rename:
            df = df.rename(columns=rename)
        clade_cols = list(clades)
    else:
        reserved = {"marker", "allele_id", "specific", "normalized_inf", "informativeness", "n"}
        clade_cols = [
            c for c in df.columns
            if str(c).strip().lower() not in reserved
            and pd.to_numeric(df[c], errors="coerce").notna().any()
        ]
        if not clade_cols:
            raise ValueError(
                "Could not infer clade columns from haplotype frequency table; provide named clade columns."
            )

    df["marker"] = df["marker"].astype(str).str.strip()
    df["allele_id"] = df["allele_id"].astype(str).str.strip()
    df = df[(df["marker"] != "") & (df["allele_id"] != "")].copy()

    dup_keys = df[df[["marker", "allele_id"]].duplicated()][["marker", "allele_id"]]
    if not dup_keys.empty:
        examples = dup_keys.head(5).astype(str).agg("/".join, axis=1).tolist()
        raise ValueError(
            "Haplotype frequency table contains duplicate marker/allele_id rows, "
            f"for example: {examples}"
        )

    for c in clade_cols:
        df[c] = pd.to_numeric(df[c], errors="coerce").fillna(0.0)

    lookup: Dict[Tuple[str, str], np.ndarray] = {}
    for _, row in df.iterrows():
        lookup[(str(row["marker"]), str(row["allele_id"]))] = np.asarray(
            [row[c] for c in clade_cols],
            dtype=float,
        )
    return lookup


def load_marker_positions(path: str, sep: Optional[str] = None, drop_chr00: bool = True) -> pd.DataFrame:
    """
    Load marker positions and normalize columns for the HMM.

    Accepted input columns are case-insensitive variants of:
      marker, chrom, pos

    Headerless three-column files are also accepted in marker/chrom/pos order.
    The returned DataFrame contains marker, chrom, pos, and Chrom_numeric.
    """
    table_sep = _guess_table_sep(path, sep)
    marker_aliases = {"marker", "markers", "locus", "loci", "id"}
    chrom_aliases = {"chrom", "chromosome", "chr"}
    pos_aliases = {"pos", "position", "bp", "start"}

    df = pd.read_csv(path, sep=table_sep)
    lower_to_col = {str(c).strip().lower(): c for c in df.columns}

    has_named_cols = (
        any(a in lower_to_col for a in marker_aliases)
        and any(a in lower_to_col for a in chrom_aliases)
        and any(a in lower_to_col for a in pos_aliases)
    )
    if not has_named_cols:
        df = pd.read_csv(path, sep=table_sep, header=None, names=["marker", "chrom", "pos"])
        lower_to_col = {str(c).strip().lower(): c for c in df.columns}

    rename = {}
    for alias in marker_aliases:
        if alias in lower_to_col:
            rename[lower_to_col[alias]] = "marker"
            break
    for alias in chrom_aliases:
        if alias in lower_to_col:
            rename[lower_to_col[alias]] = "chrom"
            break
    for alias in pos_aliases:
        if alias in lower_to_col:
            rename[lower_to_col[alias]] = "pos"
            break
    for alias in ("chrom_numeric", "chrom_num", "chromnumber", "chrom_number"):
        if alias in lower_to_col:
            rename[lower_to_col[alias]] = "Chrom_numeric"
            break

    df = df.rename(columns=rename)
    required = {"marker", "chrom", "pos"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"Marker position file is missing required columns: {sorted(missing)}")

    marker_df = df.copy()
    marker_df["marker"] = marker_df["marker"].astype(str)
    marker_df["chrom"] = marker_df["chrom"].astype(str)
    marker_df["pos"] = pd.to_numeric(marker_df["pos"], errors="coerce")
    if "Chrom_numeric" not in marker_df.columns:
        marker_df["Chrom_numeric"] = marker_df["chrom"].apply(chrom_to_numeric)
    marker_df["Chrom_numeric"] = pd.to_numeric(marker_df["Chrom_numeric"], errors="coerce")
    marker_df = marker_df.dropna(subset=["marker", "chrom", "pos", "Chrom_numeric"]).copy()
    marker_df["pos"] = marker_df["pos"].astype(int)
    marker_df["Chrom_numeric"] = marker_df["Chrom_numeric"].astype(int)

    marker_df = marker_df.sort_values(["Chrom_numeric", "pos", "marker"]).reset_index(drop=True)
    if marker_df["marker"].duplicated().any():
        marker_df = (
            marker_df.groupby("marker", as_index=False)
            .agg({"chrom": "first", "pos": lambda x: int(np.median(x)), "Chrom_numeric": "first"})
            .sort_values(["Chrom_numeric", "pos", "marker"])
            .reset_index(drop=True)
        )
    return marker_df[["marker", "chrom", "pos", "Chrom_numeric"]]


def normalize_marker_positions(marker_df: pd.DataFrame) -> pd.DataFrame:
    """Normalize an in-memory marker-position DataFrame to the columns used by the HMM."""
    df = marker_df.copy()
    lower_to_col = {str(c).strip().lower(): c for c in df.columns}
    rename = {}
    for aliases, target in (
        (("marker", "markers", "locus", "loci", "id"), "marker"),
        (("chrom", "chromosome", "chr"), "chrom"),
        (("pos", "position", "bp", "start"), "pos"),
        (("chrom_numeric", "chrom_num", "chromnumber", "chrom_number"), "Chrom_numeric"),
    ):
        for alias in aliases:
            if alias in lower_to_col:
                rename[lower_to_col[alias]] = target
                break
    df = df.rename(columns=rename)
    missing = {"marker", "chrom", "pos"} - set(df.columns)
    if missing:
        raise ValueError(f"marker_df is missing required columns: {sorted(missing)}")

    df["marker"] = df["marker"].astype(str)
    df["chrom"] = df["chrom"].astype(str)
    df["pos"] = pd.to_numeric(df["pos"], errors="coerce")
    if "Chrom_numeric" not in df.columns:
        df["Chrom_numeric"] = df["chrom"].apply(chrom_to_numeric)
    df["Chrom_numeric"] = pd.to_numeric(df["Chrom_numeric"], errors="coerce")
    df = df.dropna(subset=["marker", "chrom", "pos", "Chrom_numeric"]).copy()
    df["pos"] = df["pos"].astype(int)
    df["Chrom_numeric"] = df["Chrom_numeric"].astype(int)
    df = df[~df["chrom"].astype(str).str.lower().isin({"chr00", "0", "00"})]
    df = df.sort_values(["Chrom_numeric", "pos", "marker"]).reset_index(drop=True)
    if df["marker"].duplicated().any():
        df = (
            df.groupby("marker", as_index=False)
            .agg({"chrom": "first", "pos": lambda x: int(np.median(x)), "Chrom_numeric": "first"})
            .sort_values(["Chrom_numeric", "pos", "marker"])
            .reset_index(drop=True)
        )
    return df[["marker", "chrom", "pos", "Chrom_numeric"]]


def load_hap_genotypes(
    path: str,
    sep: Optional[str] = None,
    marker_col: Optional[str] = None,
    sample_names: Optional[Sequence[str]] = None,
    metadata_cols: Sequence[str] = ("Haplotypes",)
) -> pd.DataFrame:
    """
    Load a hap_genotype matrix as marker-indexed sample columns.

    The expected shape is one marker per row, sample names in columns, and genotype
    cells such as 1/2, 1|2, or 1/2:read_counts. Metadata columns like Haplotypes
    are ignored. Values are normalized to (allele1, allele2) tuples.
    """
    df = read_table(path, sep=sep, dtype=str)
    if df.empty:
        raise ValueError(f"hap_genotype file is empty: {path}")

    if marker_col is None:
        candidates = ["marker", "Marker", "locus", "Locus", "id", "ID"]
        marker_col = next((c for c in candidates if c in df.columns), df.columns[0])
    if marker_col not in df.columns:
        raise ValueError(f"Marker column {marker_col!r} was not found in {path}")

    metadata_lower = {c.lower() for c in metadata_cols}
    sample_cols = [
        c for c in df.columns
        if c != marker_col and str(c).strip().lower() not in metadata_lower
    ]
    if sample_names is not None:
        requested = set(sample_names)
        sample_cols = [c for c in sample_cols if c in requested]
        missing = sorted(requested - set(sample_cols))
        if missing:
            raise ValueError(f"Samples missing from hap_genotype file: {missing[:10]}")
    if not sample_cols:
        raise ValueError("No sample columns found in hap_genotype file.")

    gt = df[[marker_col] + sample_cols].copy()
    gt[marker_col] = gt[marker_col].astype(str)
    gt = gt.set_index(marker_col)
    gt.index.name = "marker"
    for col in sample_cols:
        gt[col] = gt[col].map(_normalize_hap_value)
    return gt


def load_reference_membership(
    path: str,
    sample_col: Optional[str] = None,
    clade_col: Optional[str] = None,
) -> pd.DataFrame:
    """
    Load a sample-to-clade/group table.

    Accepted sample column aliases: sample, Sample, IID, id.
    Accepted group column aliases: clade, group, population, pop, species, index.
    If no column names: assumes the first column is the sample name and the second is the clade/group.
    """
    sample_aliases = ("sample", "iid", "id", "sample_id")
    clade_aliases = ("clade", "group", "population", "pop", "species", "index")

    df = read_table(path, dtype=str, keep_default_na=False)
    lower_to_col = {str(c).strip().lower(): c for c in df.columns}
    has_named_sample_col = any(alias in lower_to_col for alias in sample_aliases)
    has_named_clade_col = any(alias in lower_to_col for alias in clade_aliases)

    if sample_col is not None or clade_col is not None or (has_named_sample_col and has_named_clade_col):
        if sample_col is None:
            sample_col = next((lower_to_col[a] for a in sample_aliases if a in lower_to_col), None)
        if clade_col is None:
            clade_col = next((lower_to_col[a] for a in clade_aliases if a in lower_to_col), None)
        if sample_col is None or clade_col is None:
            raise ValueError("Reference membership file needs sample/IID and clade/group/population columns.")
        out = df[[sample_col, clade_col]].rename(columns={sample_col: "sample", clade_col: "clade"}).copy()
    else:
        df = read_table(path, dtype=str, keep_default_na=False, header=None)
        if df.shape[1] < 2:
            raise ValueError("Reference membership file needs sample/IID and clade/group/population columns.")
        out = df.iloc[:, [0, 1]].copy()
        out.columns = ["sample", "clade"]

    out["sample"] = out["sample"].astype(str).str.strip()
    out["clade"] = out["clade"].astype(str).str.strip()
    out = out[(out["sample"] != "") & (out["clade"] != "")]
    if out["sample"].duplicated().any():
        #print(f"Note: removing {sum(out['sample'].duplicated())} duplicate reference samples")
        out = out.drop_duplicates()
    return out


def membership_clade_counts(
    membership: pd.DataFrame,
    clades: Optional[Sequence[str]] = None,
) -> pd.Series:
    """Return per-clade sample counts from the membership table."""
    counts = membership["clade"].value_counts()
    if clades is not None:
        counts = counts.reindex(list(clades), fill_value=0)
    else:
        counts = counts.sort_index()
    counts.name = "n_samples"
    return counts


def select_reference_samples(
    available_samples: Sequence[str],
    sample_to_clade: Dict[str, str],
    clades: Sequence[str],
    source_name: str,
) -> Tuple[List[str], pd.Series]:
    """
    Select samples present in both the input data and the membership table.

    Raises a clear error when the overlap is empty or does not cover every
    requested clade, because downstream reference frequencies and PCA would
    otherwise be silently misleading.
    """
    ref_samples = [s for s in available_samples if s in sample_to_clade]
    if not ref_samples:
        raise ValueError(
            f"No reference membership samples were found in the {source_name}. "
            f"Check that the membership file matches the {source_name} sample names."
        )

    matched_counts = pd.Series(
        [sample_to_clade[s] for s in ref_samples],
        dtype="object",
    ).value_counts().reindex(list(clades), fill_value=0)
    matched_counts.name = "matched_samples"

    missing_clades = [c for c in clades if int(matched_counts.get(c, 0)) == 0]
    if missing_clades:
        raise ValueError(
            f"Only {len(ref_samples)} of {len(sample_to_clade)} reference membership samples were found "
            f"in the {source_name}, and the matched set is missing these requested clades: "
            f"{missing_clades}\n"
            f"Reference membership samples present by clade:\n{matched_counts.to_string()}\n"
            f"This usually means the membership file does not correspond to the {source_name} sample names."
        )
    return ref_samples, matched_counts


def calculate_hap_reference_frequencies(
    hap_genotype_path: str,
    membership_path: str,
    clades: Sequence[str],
    sep: Optional[str] = None,
    marker_col: Optional[str] = None,
    verbose: bool = False,
) -> Tuple[pd.DataFrame, Dict[Tuple[str, str], np.ndarray], pd.DataFrame]:
    """
    Calculate clade-specific haplotype allele-ID frequencies from a hap_genotype matrix.
    """
    clades = _validate_clades(clades)
    membership = load_reference_membership(membership_path)
    hap_gt = load_hap_genotypes(hap_genotype_path, sep=sep, marker_col=marker_col)
    sample_clade = dict(zip(membership["sample"], membership["clade"]))
    missing_clades = sorted(set(sample_clade.values()) - set(clades))
    if missing_clades:
        raise ValueError(f"Membership contains clades not listed in --clades: {missing_clades}")
    ref_samples, matched_counts = select_reference_samples(
        hap_gt.columns,
        sample_clade,
        clades,
        "hap_genotype columns",
    )
    if verbose:
        print(
            f"Found {len(ref_samples)} of {len(membership)} reference membership samples "
            f"in the hap_genotype file for haplotype allele frequency estimation."
        )

    rows = []
    lookup: Dict[Tuple[str, str], np.ndarray] = {}
    for marker, row in hap_gt[ref_samples].iterrows():
        counts = {c: defaultdict(int) for c in clades}
        totals = {c: 0 for c in clades}
        for sample in ref_samples:
            c = sample_clade[sample]
            a1, a2 = row[sample]
            for a in (a1, a2):
                if a is None:
                    continue
                counts[c][str(a)] += 1
                totals[c] += 1
        alleles = sorted({a for c in clades for a in counts[c]})
        for allele in alleles:
            freqs = []
            item = {"marker": str(marker), "allele_id": str(allele)}
            for c in clades:
                f = counts[c][allele] / totals[c] if totals[c] else 0.0
                item[c] = float(f)
                freqs.append(float(f))
            rows.append(item)
            lookup[(str(marker), str(allele))] = np.asarray(freqs, dtype=float)
    return pd.DataFrame(rows), lookup, membership


def calculate_variant_reference_frequencies(
    vcf_path: str,
    membership_path: str,
    clades: Sequence[str],
    contigs: Optional[Set[str]] = None,
    require_marker: bool = True,
) -> Tuple[pd.DataFrame, Dict[Tuple[str, int, str], np.ndarray], pd.DataFrame]:
    """
    Calculate clade-specific REF/ALT allele frequencies from a reference VCF.
    """
    if VCF is None:
        raise ImportError("cyvcf2 is required for VCF reference-frequency calculation.")
    clades = _validate_clades(clades)
    membership = load_reference_membership(membership_path)
    header_samples = get_vcf_samples(vcf_path)
    sample_clade = dict(zip(membership["sample"], membership["clade"]))
    missing_clades = sorted(set(sample_clade.values()) - set(clades))
    if missing_clades:
        raise ValueError(f"Membership contains clades not listed in --clades: {missing_clades}")
    ref_samples, _ = select_reference_samples(
        header_samples,
        sample_clade,
        clades,
        "VCF header",
    )

    vcf = VCF(vcf_path, gts012=True, samples=ref_samples)
    sample_clades = [sample_clade[s] for s in vcf.samples]
    contig_list = [c for c in vcf.seqnames if contigs is None or c in contigs]
    use_marker_annotations = require_marker or vcf_has_any_marker_annotations(vcf, contig_list)
    rows = []
    freq_map: Dict[Tuple[str, int, str], np.ndarray] = {}

    try:
        for chrom in contig_list:
            for rec in vcf(f"{chrom}"):
                marker = resolve_vcf_marker_name(
                    rec,
                    require_marker=require_marker,
                    use_marker_annotations=use_marker_annotations,
                )
                if marker is None:
                    continue
                alleles = [str(rec.REF).upper()] + [str(a).upper() for a in (rec.ALT or [])]
                counts = {c: np.zeros(len(alleles), dtype=float) for c in clades}
                totals = {c: 0.0 for c in clades}
                for j, gt in enumerate(rec.genotypes):
                    c = sample_clades[j]
                    for ix in gt[:2]:
                        if ix is None or ix < 0 or ix >= len(alleles):
                            continue
                        counts[c][ix] += 1.0
                        totals[c] += 1.0
                for ix, allele in enumerate(alleles):
                    item = {
                        "CHROM": str(rec.CHROM),
                        "POS": int(rec.POS),
                        "ID": rec.ID if rec.ID not in (None, "") else ".",
                        "MARKER": str(marker),
                        "REF": str(rec.REF).upper(),
                        "ALTS": ",".join([str(a).upper() for a in (rec.ALT or [])]),
                        "ALLELE_INDEX": f"a{ix}",
                        "ALLELE": allele,
                    }
                    vec = []
                    for c in clades:
                        f = counts[c][ix] / totals[c] if totals[c] else 0.0
                        item[c] = float(f)
                        vec.append(float(f))
                    rows.append(item)
                    freq_map[(str(rec.CHROM), int(rec.POS), allele)] = np.asarray(vec, dtype=float)
    finally:
        vcf.close()
    return pd.DataFrame(rows), freq_map, membership


def _pca_from_matrix(X: np.ndarray, n_components: int = 2) -> Tuple[np.ndarray, np.ndarray]:
    X = np.asarray(X, dtype=float)
    if X.ndim != 2 or X.shape[0] < 2 or X.shape[1] < 1:
        return np.zeros((X.shape[0] if X.ndim == 2 else 0, n_components)), np.zeros(n_components)
    finite_cols = np.isfinite(X).any(axis=0)
    if not finite_cols.any():
        return np.zeros((X.shape[0], n_components)), np.zeros(n_components)
    X = X[:, finite_cols]
    X = X - np.nanmean(X, axis=0, keepdims=True)
    X = np.nan_to_num(X, nan=0.0)
    keep = X.std(axis=0) > 0
    if keep.any():
        X = X[:, keep]
    else:
        return np.zeros((X.shape[0], n_components)), np.zeros(n_components)
    U, S, _ = np.linalg.svd(X, full_matrices=False)
    coords = U[:, :n_components] * S[:n_components]
    if coords.shape[1] < n_components:
        coords = np.pad(coords, ((0, 0), (0, n_components - coords.shape[1])))
    denom = np.sum(S ** 2)
    explained = (S[:n_components] ** 2 / denom) if denom > 0 else np.zeros(min(n_components, len(S)))
    if len(explained) < n_components:
        explained = np.pad(explained, (0, n_components - len(explained)))
    return coords, explained


def resolve_pca_targets(
    pca: Optional[str],
    *,
    hap_genotype_path: Optional[str] = None,
    vcf_path: Optional[str] = None,
) -> Set[str]:
    """
    Resolve which reference inputs should be used for PCA.

    - pca=None disables PCA
    - pca="auto" runs PCA on whichever reference inputs were provided
    - pca="hap", "vcf", or "both" explicitly request those sources
    """
    if pca in (None, False):
        return set()

    available = set()
    if hap_genotype_path:
        available.add("hap")
    if vcf_path:
        available.add("vcf")
    if not available:
        raise ValueError("--pca requires --hap-genotype and/or --vcf during reference building.")

    mode = str(pca).strip().lower()
    if mode == "auto":
        return available
    if mode == "both":
        missing = [name for name in ("hap", "vcf") if name not in available]
        if missing:
            raise ValueError(
                f"--pca both requires both --hap-genotype and --vcf, but missing: {missing}"
            )
        return {"hap", "vcf"}
    if mode in {"hap", "vcf"}:
        if mode not in available:
            missing_flag = "--hap-genotype" if mode == "hap" else "--vcf"
            raise ValueError(f"--pca {mode} requires {missing_flag}.")
        return {mode}
    raise ValueError("--pca must be one of: auto, hap, vcf, both")


def validate_downsample_proportion(downsample: float) -> float:
    """Validate the PCA feature downsampling proportion."""
    downsample = float(downsample)
    if not (0.0 < downsample <= 1.0):
        raise ValueError("--downsample must be greater than 0 and less than or equal to 1.")
    return downsample


def downsample_feature_matrix(
    X: np.ndarray,
    *,
    proportion: float = 1.0,
    rng_seed: int = 0,
) -> Tuple[np.ndarray, Dict[str, int]]:
    """
    Randomly downsample PCA feature columns for speed.

    Downsampling is deterministic for a given matrix width because a fixed seed is
    used by default.
    """
    X = np.asarray(X, dtype=float)
    n_features = int(X.shape[1]) if X.ndim == 2 else 0
    proportion = validate_downsample_proportion(proportion)
    if X.ndim != 2 or n_features == 0 or proportion >= 1.0:
        return X, {
            "features_before_downsample": n_features,
            "features_after_downsample": n_features,
        }

    keep_n = max(1, int(round(n_features * proportion)))
    keep_n = min(keep_n, n_features)
    if keep_n == n_features:
        return X, {
            "features_before_downsample": n_features,
            "features_after_downsample": n_features,
        }

    rng = np.random.default_rng(rng_seed)
    keep_idx = np.sort(rng.choice(n_features, size=keep_n, replace=False))
    return X[:, keep_idx], {
        "features_before_downsample": n_features,
        "features_after_downsample": keep_n,
    }


def _hap_sample_feature_matrix(hap_gt: pd.DataFrame, samples: Sequence[str]) -> Tuple[np.ndarray, List[str]]:
    feature_keys = []
    for marker, row in hap_gt[list(samples)].iterrows():
        alleles = sorted({a for sample in samples for a in row[sample] if a is not None})
        feature_keys.extend([(str(marker), str(a)) for a in alleles])
    col = {k: i for i, k in enumerate(feature_keys)}
    X = np.zeros((len(samples), len(feature_keys)), dtype=float)
    for i, sample in enumerate(samples):
        for marker, val in hap_gt[sample].items():
            for a in val:
                if a is not None and (str(marker), str(a)) in col:
                    X[i, col[(str(marker), str(a))]] += 1.0
    return X, [f"{m}|{a}" for m, a in feature_keys]


def _vcf_sample_feature_matrix(vcf_path: str, samples: Sequence[str], contigs: Optional[Set[str]] = None) -> np.ndarray:
    if VCF is None:
        raise ImportError("cyvcf2 is required for VCF PCA.")
    vcf = VCF(vcf_path, gts012=True, samples=list(samples))
    rows = []
    try:
        contig_list = [c for c in vcf.seqnames if contigs is None or c in contigs]
        for chrom in contig_list:
            for rec in vcf(f"{chrom}"):
                dosage = []
                for gt in rec.genotypes:
                    vals = [ix for ix in gt[:2] if ix is not None and ix >= 0]
                    dosage.append(float(sum(1 for ix in vals if ix > 0)) if vals else np.nan)
                rows.append(dosage)
    finally:
        vcf.close()
    if not rows:
        return np.zeros((len(samples), 0), dtype=float)
    return np.asarray(rows, dtype=float).T


def write_reference_pca_and_metrics(
    *,
    matrix: np.ndarray,
    samples: Sequence[str],
    membership: pd.DataFrame,
    outdir: str,
    prefix: str = "reference",
) -> Dict[str, str]:
    """
    Write PCA coordinates, a PCA plot, and simple group-differentiation metrics.
    """
    Path(outdir).mkdir(parents=True, exist_ok=True)
    sample_to_clade = dict(zip(membership["sample"], membership["clade"]))
    groups = [sample_to_clade[s] for s in samples]
    coords, explained = _pca_from_matrix(matrix, n_components=2)
    pca_df = pd.DataFrame({
        "sample": list(samples),
        "clade": groups,
        "PC1": coords[:, 0] if coords.size else [],
        "PC2": coords[:, 1] if coords.size else [],
    })
    pca_path = Path(outdir) / f"{prefix}_pca.tsv"
    pca_df.to_csv(pca_path, sep="\t", index=False)

    metrics = {
        "n_samples": int(len(samples)),
        "n_features": int(matrix.shape[1] if matrix.ndim == 2 else 0),
        "pc1_explained_variance": float(explained[0]) if len(explained) else 0.0,
        "pc2_explained_variance": float(explained[1]) if len(explained) > 1 else 0.0,
    }
    centroids = pca_df.groupby("clade")[["PC1", "PC2"]].mean()
    if len(centroids) > 1:
        dists = []
        names = list(centroids.index)
        for i in range(len(names)):
            for j in range(i + 1, len(names)):
                dists.append(float(np.linalg.norm(centroids.loc[names[i]] - centroids.loc[names[j]])))
        within = []
        for _, row in pca_df.iterrows():
            within.append(float(np.linalg.norm(row[["PC1", "PC2"]].to_numpy(dtype=float) - centroids.loc[row["clade"]].to_numpy(dtype=float))))
        metrics["mean_pairwise_centroid_distance_pc"] = float(np.mean(dists)) if dists else 0.0
        metrics["mean_within_clade_distance_pc"] = float(np.mean(within)) if within else 0.0
        metrics["centroid_to_within_distance_ratio_pc"] = metrics["mean_pairwise_centroid_distance_pc"] / max(metrics["mean_within_clade_distance_pc"], 1e-12)
    else:
        metrics["mean_pairwise_centroid_distance_pc"] = 0.0
        metrics["mean_within_clade_distance_pc"] = 0.0
        metrics["centroid_to_within_distance_ratio_pc"] = 0.0

    metrics_path = Path(outdir) / f"{prefix}_differentiation_metrics.json"
    metrics_path.write_text(json.dumps(metrics, indent=2, sort_keys=True) + "\n")
    plot_path = Path(outdir) / f"{prefix}_pca.png"
    try:
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(figsize=(7, 5))
        for clade, sub in pca_df.groupby("clade"):
            ax.scatter(sub["PC1"], sub["PC2"], label=str(clade), s=40)
        ax.set_xlabel(f"PC1 ({metrics['pc1_explained_variance']:.1%})")
        ax.set_ylabel(f"PC2 ({metrics['pc2_explained_variance']:.1%})")
        ax.legend(title="Clade")
        fig.tight_layout()
        fig.savefig(plot_path, dpi=250)
        plt.close(fig)
    except Exception as e:
        plot_path = Path(outdir) / f"{prefix}_pca_plot_error.txt"
        plot_path.write_text(f"Could not create PCA plot: {e}\n")
    return {"pca": str(pca_path), "metrics": str(metrics_path), "plot": str(plot_path)}


def build_reference_files(
    *,
    outdir: str,
    clades: Optional[Sequence[str]],
    membership_path: str,
    hap_genotype_path: Optional[str] = None,
    vcf_path: Optional[str] = None,
    contigs: Optional[Set[str]] = None,
    require_marker: Optional[bool] = None,
    verbose: bool = False,
    pca: Optional[str] = None,
    downsample: float = 1.0,
    prefix: str = "reference",
) -> Dict[str, Optional[str]]:
    """
    Build reference allele-frequency/profile files from reference inputs.
    """

    Path(outdir).mkdir(parents=True, exist_ok=True)
    outputs: Dict[str, Optional[str]] = {}
    membership = load_reference_membership(membership_path)

    if clades:
        clades = _validate_clades(clades)
    else:
        clades = sorted(set(membership["clade"]))
    effective_require_marker = resolve_require_marker(require_marker, vcf_path, hap_genotype_path)
    pca_targets = resolve_pca_targets(
        pca,
        hap_genotype_path=hap_genotype_path,
        vcf_path=vcf_path,
    )
    downsample = validate_downsample_proportion(downsample)

    sample_to_clade = dict(zip(membership["sample"], membership["clade"]))

    if verbose:
        counts = membership_clade_counts(membership, clades)
        print(
            f"Loaded reference membership file with {len(membership)} samples:\n"
            f"{counts.to_string()}"
        )

    if hap_genotype_path:
        hap_freq_df, hap_lookup, membership = calculate_hap_reference_frequencies(
            hap_genotype_path, membership_path, clades, verbose=verbose)
        hap_freq_path = Path(outdir) / f"{prefix}_hap_allele_frequencies.tsv"
        hap_inf_path = Path(outdir) / f"{prefix}_hap_allele_informativeness.tsv"
        hap_freq_df.to_csv(hap_freq_path, sep="\t", index=False)
        if verbose: print(f"Calculated haplotype allele frequencies for clades: {clades}")
        

        hap_allele_informativeness(hap_lookup, clades).to_csv(hap_inf_path, sep="\t", index=False)
        if verbose: print(f"Calculated haplotype allele informativeness.")

        if "hap" in pca_targets:
            if verbose: print(f"Loading hap_genotype matrix for PCA...")
            hap_gt = load_hap_genotypes(hap_genotype_path)
            ref_samples, _ = select_reference_samples(
                hap_gt.columns,
                sample_to_clade,
                clades,
                "hap_genotype columns",
            )
            X, _ = _hap_sample_feature_matrix(hap_gt, ref_samples)
            X, downsample_info = downsample_feature_matrix(X, proportion=downsample)
            if verbose:
                print(
                    f"Running PCA on the hap_genotype reference sample matrix: "
                    f"{X.shape[0]} samples x {X.shape[1]} haplotype-allele features...")
                if downsample != 1:
                    print(
                        f"Downsampled feature matrix to {downsample:g} of features "
                        f"({downsample_info['features_after_downsample']} of "
                        f"{downsample_info['features_before_downsample']} total features retained)."
                    )
            outputs.update(write_reference_pca_and_metrics(
                matrix=X, samples=ref_samples, membership=membership, outdir=outdir, prefix=f"{prefix}_hap"
            ))
        outputs.update({
            "hap_frequencies": str(hap_freq_path),
            "hap_informativeness": str(hap_inf_path),
        })

    if vcf_path:
        variant_df, _, membership = calculate_variant_reference_frequencies(
            vcf_path, membership_path, clades, contigs=contigs, require_marker=effective_require_marker)
        if verbose: print(f"Processed VCF input: {len(variant_df)} variants")

        variant_inf = variant_allele_informativeness(variant_df, clades)
        locus_inf = allele_informativeness_by_locus(variant_df, clades)
        if verbose: print(f"Calculated informativeness metrics for VCF variants")

        variant_path = Path(outdir) / f"{prefix}_variant_profiles.tsv"
        locus_path = Path(outdir) / f"{prefix}_variant_locus_informativeness.tsv"
        variant_inf.to_csv(variant_path, sep="\t", index=False)
        locus_inf.to_csv(locus_path, sep="\t", index=False)
        header_samples = get_vcf_samples(vcf_path)
        ref_samples, matched_counts = select_reference_samples(
            header_samples,
            sample_to_clade,
            clades,
            "VCF header",
        )
        if "vcf" in pca_targets:
            if verbose: print(f"Loading VCF matrix for PCA...")
            if verbose:
                print(
                    f"Found {len(ref_samples)} of {len(membership)} reference membership samples "
                    f"in the VCF file."
                )
            X = _vcf_sample_feature_matrix(vcf_path, ref_samples, contigs=contigs)
            X, downsample_info = downsample_feature_matrix(X, proportion=downsample)
            if verbose:
                if downsample != 1:
                    print(
                        f"Running PCA on the VCF reference sample matrix: "
                        f"{X.shape[0]} samples x {X.shape[1]} variant-dosage features "
                        f"after downsampling from {downsample_info['features_before_downsample']} total "
                        f"features (downsample={downsample:g})..."
                    )
                else:
                    print(
                        f"Running PCA on the VCF reference sample matrix: "
                        f"{X.shape[0]} samples x {X.shape[1]} variant-dosage features..."
                    )
            pca_outputs = write_reference_pca_and_metrics(
                matrix=X, samples=ref_samples, membership=membership, outdir=outdir, prefix=f"{prefix}_vcf")
            outputs.update({
                "vcf_pca": pca_outputs["pca"],
                "vcf_pca_plot": pca_outputs["plot"],
                "vcf_metrics": pca_outputs["metrics"],
            })
        outputs.update({
            "variant_profiles": str(variant_path),
            "variant_locus_informativeness": str(locus_path),
        })
    

    if not hap_genotype_path and not vcf_path:
        raise ValueError("Provide hap_genotype_path, vcf_path, or both to build reference files.")
    return outputs


def get_vcf_samples(vcf_path: str) -> List[str]:
    """Return sample names from a VCF/BCF header without loading all genotypes."""
    if VCF is None:
        raise ImportError("cyvcf2 is required for VCF input. Install it with `pip install cyvcf2`.")
    vcf = VCF(vcf_path, gts012=True)
    try:
        return list(vcf.samples)
    finally:
        vcf.close()


def load_chrom_lengths(path: str, sep: Optional[str] = None) -> pd.DataFrame:
    """
    Load chromosome lengths for plotting.

    Accepts FASTA .fai files or delimited tables with chrom/pos-like columns.
    Returns columns chrom and pos.
    """
    lower = str(path).lower()
    if lower.endswith((".fai", ".fai.gz")):
        df = pd.read_csv(path, sep="\t", header=None, usecols=[0, 1], names=["chrom", "pos"])
    else:
        df = read_table(path, sep=sep)
        lower_to_col = {str(c).strip().lower(): c for c in df.columns}
        chrom_col = next((lower_to_col[a] for a in ("chrom", "chr", "chromosome") if a in lower_to_col), df.columns[0])
        pos_col = next((lower_to_col[a] for a in ("pos", "length", "bp", "end") if a in lower_to_col), df.columns[1])
        df = df[[chrom_col, pos_col]].rename(columns={chrom_col: "chrom", pos_col: "pos"})
    df["chrom"] = df["chrom"].astype(str)
    df["pos"] = pd.to_numeric(df["pos"], errors="coerce")
    return df.dropna(subset=["chrom", "pos"]).copy()


def parse_samples_arg(samples: Optional[Sequence[str]] = None, samples_file: Optional[str] = None) -> Optional[List[str]]:
    """Combine repeated/comma-separated sample CLI values and an optional sample file."""
    out: List[str] = []
    if samples:
        for item in samples:
            out.extend([x.strip() for x in str(item).split(",") if x.strip()])
    if samples_file:
        with open(samples_file) as f:
            for line in f:
                s = line.strip()
                if s and not s.startswith("#"):
                    out.append(s)
    if not out:
        return None
    seen = set()
    unique = []
    for s in out:
        if s not in seen:
            unique.append(s)
            seen.add(s)
    return unique


def choose_samples(
    vcf_samples: Optional[Sequence[str]] = None,
    hap_samples: Optional[Sequence[str]] = None,
    requested_samples: Optional[Sequence[str]] = None,
    all_samples: bool = False,
    sample_source: str = "auto"
) -> List[str]:
    """
    Pick sample names to run.

    If both VCF and hap_genotype inputs are present, auto mode uses the intersection
    when running all samples so combined evidence is used consistently. Explicitly requested
    samples may be present in either input; the available source(s) are used.
    """
    vcf_samples = list(vcf_samples or [])
    hap_samples = list(hap_samples or [])
    vcf_set = set(vcf_samples)
    hap_set = set(hap_samples)
    have_both_input_types = bool(vcf_samples) and bool(hap_samples)

    if requested_samples:
        available = vcf_set | hap_set
        missing = [s for s in requested_samples if s not in available]
        if missing:
            raise ValueError(f"Requested samples not found in any input: {missing[:10]}")
        return list(requested_samples)

    if have_both_input_types:
        if sample_source == "union":
            ordered = vcf_samples + [s for s in hap_samples if s not in vcf_set]
        else:
            common = vcf_set & hap_set
            ordered = [s for s in vcf_samples if s in common] + [s for s in hap_samples if s in common and s not in vcf_set]
    else:
        ordered = vcf_samples or hap_samples

    if not ordered:
        raise ValueError("No samples found in the provided inputs.")
    return ordered if all_samples else [ordered[0]]


def parse_contigs(contigs: Optional[str]) -> Optional[Set[str]]:
    if contigs is None or str(contigs).strip() == "":
        return None
    values: List[str] = []
    for piece in str(contigs).replace(";", ",").split(","):
        piece = piece.strip()
        if piece:
            values.append(piece)
    return set(values) if values else None



# ---------------------------
# Emissions from per-variant profiles
# ---------------------------

def states_from_clades(clades: List[str]) -> List[Tuple[str, str]]:
    return list(combinations_with_replacement(clades, 2))

def build_logB_from_hapgeno(
    test_series, freq_lookup, clades, states, marker_order,
    strict_alleles=None, strict_boost=HMM_DEFAULTS["strict_boost"],      # per-matching hap copy (log-additive)
    e_geno=HMM_DEFAULTS["e_geno"],                                 # genotype error / miscoding rate
    eps=1e-12,                                   # more dynamic range
    # Homozygous tempering: flatten when allele isn't clearly clade-specific
    e_homo=HMM_DEFAULTS["e_homo"],                               # heterozygous undercalling rate
    hom_soften_delta=HMM_DEFAULTS["hom_soften_delta"],                       # min (top - second) diff to avoid flattening
    hom_soften_width=HMM_DEFAULTS["hom_soften_width"],                       # slope of logistic around the threshold
    hom_min_mix=HMM_DEFAULTS["hom_min_mix"],                            # how much to mix with neutral when weak
    hom_neutral=HMM_DEFAULTS["hom_neutral"],                             #neutral target for mixing
    cap_total_boost_per_marker=HMM_DEFAULTS["cap_total_boost_per_marker"],       # cap total boost added to ANY state (in logs)
    debug=None):
    """
    Emission log-likelihoods for unordered diploid states using an error-aware model.
    - For homozygotes, evidence is *softened* unless the allele is clearly clade-diagnostic.
    - Heterozygotes use standard unordered diploid likelihood with genotype error.
    - Strict allele boosts are copy-aware, but capped so they can't dominate weak evidence.
    Returns: logB[T, S] max-centered per marker.
    """
    K, S = len(clades), len(states)
    clade2idx = {c: i for i, c in enumerate(clades)}
    states_arr = np.asarray(states, dtype=object)  # (S, 2)
    ii = np.fromiter((clade2idx[a] for a, _ in states), dtype=int, count=S)
    jj = np.fromiter((clade2idx[b] for _, b in states), dtype=int, count=S)

    # precompute (1-2e) term
    one_minus_2e = (1.0 - 2.0 * e_geno)

    # strict boosts
    use_boosts = strict_alleles is not None
    if use_boosts:
        if strict_boost is None or strict_boost < 1.0:
            raise ValueError("strict_boost must be >= 1.0 when strict_alleles is provided.")
        LOG_BOOST = np.log(strict_boost)
        def diag_clade(marker, allele):
            return strict_alleles.get((marker, allele))
    else:
        LOG_BOOST = 0.0
        diag_clade = lambda *_: None

    def allele_freq(marker, allele):
        # Vector of length K; if missing, return small flat mass (not zeros)
        f = freq_lookup.get((marker, allele))
        if f is None:
            return np.full(K, eps, dtype=float)
        f = np.asarray(f, dtype=float)
        # clip to (eps, 1-eps) for stability
        return np.clip(f, eps, 1.0 - eps)

    def soften_mix_for_homozygote(f_a):
        top = float(np.max(f_a))
        second = float(np.partition(f_a, -2)[-2]) if f_a.size > 1 else 0.0
        sep = max(0.0, top - second)
        x = (sep - hom_soften_delta) / max(1e-9, hom_soften_width)
        w = 1.0 / (1.0 + np.exp(-x))            # 0..1
        return hom_min_mix + (1.0 - hom_min_mix) * w 

    T = len(marker_order)
    logB = np.empty((T, S), dtype=float)

    from collections import defaultdict
    site2alleles = defaultdict(list)
    for mk, al in freq_lookup.keys():
        site2alleles[mk].append(al)

    for t, m in enumerate(marker_order):
        alleles = test_series.get(m, None)

        # fallback: uninformative row if malformed genotype
        if not isinstance(alleles, (list, tuple)) or len(alleles) != 2:
            logB[t] = 0.0
            continue

        a, b = str(alleles[0]), str(alleles[1])
        f_a = allele_freq(m, a)  # shape (K,)
        f_b = allele_freq(m, b)

        # Per-haplotype observation probabilities with error:
        # P(obs a | clade k) = e + (1-2e)*f_a[k]
        pa = e_geno + one_minus_2e * f_a      # per-hap P(obs 'a' | clade k)
        pb = e_geno + one_minus_2e * f_b      # per-hap P(obs 'b' | clade k)

        if a == b:
            # Homozygous case: account for het. undercalling 
            # and temper evidence if allele 'a' is not clearly clade-specific

            # Genotype homozygote likelihood (unordered state i,j):
            p_hom = pa[ii] * pa[jj]

            # Only iterate alleles for THIS marker
            alleles_here = site2alleles.get(m, ())
            other_vec = None
            other_cnt = 0
            for y in alleles_here:
                if y == a:
                    continue
                f_y = allele_freq(m, y)
                py  = e_geno + one_minus_2e * f_y
                other_vec = py if other_vec is None else (other_vec + py)
                other_cnt += 1

            if other_cnt > 0:
                other_vec *= (1.0 / other_cnt)     # mean partner emission over other alleles
                # Unordered het emission if the true genotype were a/other
                p_het = pa[ii] * other_vec[jj] + pa[jj] * other_vec[ii]
                # Mix in dropout: true het miscalled as hom with probability e_homo
                p_dropout = (1.0 - e_homo) * p_hom + e_homo * p_het
            else:
                # No known alternative allele at this site → fall back to base homozygote
                p_dropout = p_hom

            # Informativeness-aware tempering:
            w = soften_mix_for_homozygote(f_a) #weight for neutral distribution
            U_hom = float(hom_neutral) * float(hom_neutral)  # neutral genotype mass for hom
            # Mix towards "neutral" to avoid overconfidence on weakly-diagnostic homozygotes
            p = w * p_dropout + (1.0 - w) * U_hom
            
        else:
            # --- Heterozygous case: standard unordered sum of the two assignments ---
            p = pa[ii] * pb[jj] + pa[jj] * pb[ii]

        # Numeric safety
        p = np.maximum(p, eps)
        lp = np.log(p)

        # --- Strict diagnostic boosts  ---
        if use_boosts:
            reqs = set()
            ra = diag_clade(m, a)
            rb = diag_clade(m, b)
            if ra: reqs.add(str(ra))
            if rb: reqs.add(str(rb))

            if reqs:
                # Count how many diagnostic copies match each state (0,1,2)
                # For homozygous, if both alleles map to same clade, that's 2 copies support.
                # For heterozygous with two distinct diagnostics, can yield up to 2 as well.
                match_count = np.zeros(S, dtype=int)
                for req in reqs:
                    match_count += ((states_arr[:,0] == req) | (states_arr[:,1] == req)).astype(int)

                boost = match_count * LOG_BOOST
                # Cap total boost so boosts can’t overwhelm contradictory frequencies
                if cap_total_boost_per_marker is not None:
                    boost = np.minimum(boost, cap_total_boost_per_marker)
                lp = lp + boost

        # Max-center for numerical stability and better contrast
        lp = lp - np.max(lp)
        logB[t] = lp

    return logB
    


def build_logB_unordered_from_variants_v4(
    variants: List["VCFVariant"],
    freq_map: Dict[Tuple[str,int,str], np.ndarray],
    clades: List[str],
    states: List[Tuple[str,str]],
    marker_order: List[str],
    prof_df: pd.DataFrame, # per-allele informativeness (nats)
    e_geno: float = HMM_DEFAULTS["e_geno"],
    e_homo: float = HMM_DEFAULTS["e_homo"],
    eps: float = 1e-12,
    hom_soften_delta: float = HMM_DEFAULTS["hom_soften_delta"],
    hom_soften_width: float = HMM_DEFAULTS["hom_soften_width"],
    hom_min_mix: float = HMM_DEFAULTS["hom_min_mix"],
    hom_neutral: float = HMM_DEFAULTS["hom_neutral"],
    debug: bool = False
) -> np.ndarray:
    """
    Build sample-specific emission log-likelihoods per marker from variant genotypes,
    weighting each variant site by its allele-specific informativeness drawn from `prof_df`.

    `prof_df` with per-allele specific informativeness in nats (base-e)

    Site weight for marker-level combine:
        S_site = s(a1) + s(a2)  (sum over the individual's two alleles at that site),
        w_site = (S_site ** weight_gamma) / sum_sites ( ... )   (normalized to sum 1).
    If all S_site are 0 or missing, falls back to equal weights.

    Returns
    -------
    logB : np.ndarray, shape (T, S)
        One row per marker in `marker_order`, columns are unordered diploid states.
    """
    # --- setup
    one_minus_2e = (1.0 - 2.0 * e_geno)
    K, S = len(clades), len(states)
    clade2idx = {c: i for i, c in enumerate(clades)}
    ii = np.fromiter((clade2idx[a] for a, _ in states), dtype=int, count=S)
    jj = np.fromiter((clade2idx[b] for _, b in states), dtype=int, count=S)

    # allele frequency accessor
    af_cache: Dict[Tuple[str,int,str], np.ndarray] = {}
    def allele_freq(chrom, pos, allele):
        key = (str(chrom), int(pos), str(allele).upper())
        f = af_cache.get(key)
        if f is None:
            f = freq_map.get(key)
            if f is None:
                f = np.full(K, eps, dtype=float)
            f = np.clip(np.asarray(f, dtype=float), eps, 1.0 - eps)
            af_cache[key] = f
        return f

    # homozygote tempering function
    def soften_mix_for_homozygote(f_a):
        top = float(np.max(f_a))
        second = float(np.partition(f_a, -2)[-2]) if f_a.size > 1 else 0.0
        sep = max(0.0, top - second)
        x = (sep - hom_soften_delta) / max(1e-9, hom_soften_width)
        w = 1.0 / (1.0 + np.exp(-x))
        return hom_min_mix + (1.0 - hom_min_mix) * w

    # variant list grouped by marker for this sample
    by_marker: Dict[str, List["VCFVariant"]] = {}
    for v in variants:
        key = str(v.marker).strip() if v.marker is not None else ""
        if key:
            by_marker.setdefault(key, []).append(v)

    # build lookup table for per-allele "specific" in nats from prof_df
    if not {"CHROM","POS","ALLELE"}.issubset(prof_df.columns):
        raise ValueError("prof_df must have columns CHROM, POS, ALLELE, specific")
    if "specific" not in prof_df.columns:
        raise ValueError("prof_df must have 'specific' (nats)'.")
    else:
        prof_work = prof_df

    prof_work = prof_work.copy()
    prof_work["CHROM"] = prof_work["CHROM"].astype(str)
    prof_work["POS"]   = prof_work["POS"].astype(int)
    prof_work["ALLELE"] = prof_work["ALLELE"].astype(str)

    spec_LUT = pd.Series(
        prof_work["specific"].to_numpy(),
        index=pd.MultiIndex.from_frame(prof_work[["CHROM","POS","ALLELE"]])
    ).to_dict()

    # also gather other alleles present at each site (for hom tempering partner-mean)
    site2alleles = defaultdict(list)
    for (c, p, al) in freq_map.keys():
        site2alleles[(str(c), int(p))].append(str(al).upper())

    # ---- main loop
    T = len(marker_order)
    logB = np.zeros((T, S), dtype=float)

    for t, m in enumerate(marker_order):
        varlist = by_marker.get(m, [])
        if not varlist:
            continue

        site_logps: List[np.ndarray] = []
        site_weights: List[float] = []

        for v in varlist:
            # skip missing GT
            if not hasattr(v, "gt_alleles") or v.gt_alleles is None:
                continue
            if not isinstance(v.gt_alleles, (tuple, list)) or len(v.gt_alleles) != 2:
                continue

            chrom = str(v.chrom); pos = int(v.pos)
            a = str(v.gt_alleles[0]).upper()
            b = str(v.gt_alleles[1]).upper()

            # allele freqs -> genotype emission for unordered states
            f_a = allele_freq(chrom, pos, a)
            f_b = allele_freq(chrom, pos, b)
            pa = e_geno + one_minus_2e * f_a
            pb = e_geno + one_minus_2e * f_b

            if a == b:
                # homozygote: global tempering to account for het. undercalling errors
                # and site-specific tempering for uninformative sites
                p_hom = pa[ii] * pa[jj]

                site_key = (chrom, pos)
                alleles_here = site2alleles.get(site_key, ())
                other_vec = None; other_cnt = 0
                for y in alleles_here:
                    if y == a: continue
                    py = e_geno + one_minus_2e * allele_freq(chrom, pos, y)
                    other_vec = py if other_vec is None else (other_vec + py)
                    other_cnt += 1
                if other_cnt > 0:
                    other_vec *= (1.0 / other_cnt)
                    p_het = pa[ii] * other_vec[jj] + pa[jj] * other_vec[ii]
                    p_dropout = (1.0 - e_homo) * p_hom + e_homo * p_het
                else:
                    p_dropout = p_hom

                w_sep = soften_mix_for_homozygote(f_a)
                U_hom = float(hom_neutral) * float(hom_neutral)
                p = w_sep * p_dropout + (1.0 - w_sep) * U_hom
            else:
                # heterozygote unordered emission
                p = pa[ii] * pb[jj] + pa[jj] * pb[ii]

            p = np.maximum(p, eps)
            site_logps.append(np.log(p))

            # ---- site weight from allele-specific informativeness
            s1 = float(spec_LUT.get((chrom, pos, a), 0.0))
            s2 = float(spec_LUT.get((chrom, pos, b), 0.0))
            S_site = max(0.0, s1 + s2) # additive in nats
            site_weights.append(S_site)

        if not site_logps:
            continue

        w = np.asarray(site_weights, dtype=float)
        if not np.any(w > 0):
            w = np.ones_like(w)
        w = w / w.sum()

        # weighted geometric mean across variant sites in this marker
        lp = np.tensordot(w, np.vstack(site_logps), axes=(0, 0))

        # stabilize
        if np.all(np.isfinite(lp)):
            lp -= np.max(lp)
        logB[t] = lp

    return logB





# Combine haplotype ID and variant-based logB

def combine_logB_mixture(
    logB_fromhaps: np.ndarray,            # (T,S)
    logB_fromvariants: np.ndarray,        # (T,S)
    hap_known_ref: np.ndarray,            # (T,) bool
    hap_score_norm: np.ndarray,           # (T,) from hap_strength_metrics_sample(..)
    var_score_norm: np.ndarray,           # (T,) from sum_In_per_marker_sample(..)
    *,
    floor: float = -60.0,
    # gate weights
    b0: float = 0.0,      # baseline tilt toward hap-ID (positive → more hap)
    bH: float = 1.0,      # sensitivity to hap score
    bV: float = 1.0,      # sensitivity to variant score (subtracts)
    w_min: float = 0.05,
    w_max: float = 0.95,
    verbose: bool = False,
    debug: bool = False) -> Tuple[np.ndarray, np.ndarray]:
    """
    p = wH * softmax(LH) + (1-wH) * softmax(LV),
    with wH = sigmoid(b0 + bH*hap_score_norm - bV*var_score_norm),
    clamped to [w_min, w_max], and forced to 0 when hap is unknown in reference.
    """
    T, S = logB_fromhaps.shape
    assert logB_fromvariants.shape == (T, S)
    assert hap_known_ref.shape == (T,)
    assert hap_score_norm.shape == (T,)
    assert var_score_norm.shape == (T,)

    def _row_norm(L):
        X = np.array(L, copy=True)
        X[~np.isfinite(X)] = floor
        X -= X.max(axis=1, keepdims=True)
        return np.clip(X, floor, 0.0)

    # Normalize logs → probs
    LH = _row_norm(logB_fromhaps)
    LV = _row_norm(logB_fromvariants)
    PH = np.exp(LH); PH /= PH.sum(axis=1, keepdims=True)
    PV = np.exp(LV); PV /= PV.sum(axis=1, keepdims=True)

    # Gate
    x = b0 + bH * np.clip(hap_score_norm, 0, 1) - bV * np.clip(var_score_norm, 0, 1)
    wH_soft = 1.0 / (1.0 + np.exp(-x))
    wH = np.clip(wH_soft, w_min, w_max)
    wH = np.where(hap_known_ref.astype(bool), wH, 0.0)

    # Mix back to log space
    P = wH[:, None] * PH + (1.0 - wH)[:, None] * PV
    L = np.log(np.maximum(P, 1e-300))
    L -= L.max(axis=1, keepdims=True)
    L = np.maximum(L, floor)

    # Print stats
    if verbose==True:
        print(f"average weighting of haplotype ID vs variants: {round(wH.mean(), 2)}")

    return L, wH



def uninformative_marker_breakdown(
    logB, marker_order, variants, freq_map,
    hap_test_series=None,
    hap_freq_lookup=None,
    tol=1e-6, round_decimals=6,
    entropy_frac_thresh=0.97 # “flat” if entropy close to uniform
):
    """
    Classify & summarize per-marker informativeness.
      - variant_missing_pct: % of variant sites in the marker with missing diploid GT
      - hap_missing_genotype: haplotype ID missing for this sample at marker
      - hap_missing_in_reference: hap ID present in sample but not found in reference set
      - flat: “relatively flat” emission (high entropy)
    """

    T, S = logB.shape
    if len(marker_order) != T:
        raise ValueError("marker_df length must match logB rows.")

    # --- Prob-space stats from logB ---
    row_max = logB.max(axis=1, keepdims=True)
    logits  = logB - row_max
    probs   = np.exp(np.clip(logits, -60, 0))
    probs  /= probs.sum(axis=1, keepdims=True)

    with np.errstate(divide='ignore', invalid='ignore'):
        ent = -(probs * np.log(np.clip(probs, 1e-300, 1))).sum(axis=1)
        ent_norm = ent / np.log(S)

    # “Relatively flat” if near-uniform 
    flat = (ent_norm >= float(entropy_frac_thresh)) 

    # --- Group variants by marker (for missingness & reference checks) ---
    by_marker = {}
    for v in variants:
        by_marker.setdefault(str(v.marker), []).append(v)

    allele_not_in_reference = np.zeros(T, dtype=bool)  # any observed allele lacks profile (variant-level)
    missing_genotype        = np.zeros(T, dtype=bool)  # no complete diploid GT at any site (variant-level)
    no_variants             = np.zeros(T, dtype=bool)
    variant_missing_pct     = np.zeros(T, dtype=float) # % missing GT among variant sites in the marker amplicon

    for t, m in enumerate(marker_order):
        varlist = by_marker.get(m, [])
        if not varlist:
            no_variants[t] = True
            missing_genotype[t] = True
            variant_missing_pct[t] = 100.0
            continue

        n_sites = len(varlist)
        n_missing = 0
        has_any_genotyped_variant = False
        missing_in_profiles = False

        for v in varlist:
            a, b = v.gt_alleles
            # Handle missing genotype cases (None, ".", or -1 encodings)
            if (a is None or b is None
                or a in {".", "./.", "-1"} or b in {".", "./.", "-1"}
                or (isinstance(a, (int, float)) and a == -1)
                or (isinstance(b, (int, float)) and b == -1) ):
                    n_missing += 1
                    continue
            has_any_genotyped_variant = True

            key_a = (str(v.chrom), int(v.pos), str(a).upper())
            key_b = (str(v.chrom), int(v.pos), str(b).upper())
            if (key_a not in freq_map) or (key_b not in freq_map):
                missing_in_profiles = True
                

        variant_missing_pct[t] = (n_missing / n_sites) * 100.0

        if not has_any_genotyped_variant:
            missing_genotype[t] = True
        if has_any_genotyped_variant and missing_in_profiles:
            allele_not_in_reference[t] = True

    # --- Haplotype ID diagnostics (sample hap presence & reference presence) ---
    hap_missing_genotype     = np.zeros(T, dtype=bool)  # sample lacks hap ID / ill-formed
    hap_missing_in_reference = np.zeros(T, dtype=bool)  # sample hap present but absent from reference

    if hap_test_series is not None:
        # normalize index lookup once
        ht = hap_test_series
        # check hap allele ID presence in reference set
        ref_lookup = hap_freq_lookup if hap_freq_lookup is not None else {}

        for t, m in enumerate(marker_order):
            a, b = _normalize_hap_value(ht.get(m, None))
            if a is None and b is None:
                hap_missing_genotype[t] = True
                continue

            # If hap present but not in reference set:
            in_ref_a = ((m, a) in ref_lookup)
            in_ref_b = ((m, b) in ref_lookup)
            # Treat as "missing in ref" if neither allele has a profile
            if not (in_ref_a or in_ref_b):
                hap_missing_in_reference[t] = True

    # Check whether both haplotype ID and variants are missing 
    all_variants_missing = ((variant_missing_pct >= 100.0 - tol) |  # every site missing
                            missing_genotype |                      # no complete diploid GT anywhere
                            no_variants   )                       # marker had zero variant records

    hap_and_all_variants_missing = hap_missing_genotype & all_variants_missing

    variant_missing_pct_nonmissing = variant_missing_pct[variant_missing_pct!=100]

    # --- Assemble details DataFrame ---
    details = pd.DataFrame({
        "Marker": marker_order,
        # Variant-level flags
        "variant_missing_pct": np.round(variant_missing_pct, round_decimals),
        "variant_no_variants": no_variants,
        "variant_missing_genotype": missing_genotype,
        "variant_allele_not_in_reference": allele_not_in_reference,
        # Hap-level flags
        "hap_missing_genotype": hap_missing_genotype,
        "hap_missing_in_reference": hap_missing_in_reference,
        "missing_both": hap_and_all_variants_missing,
        # Emission flatness
        "ent": ent,
        "flat": flat  })

    # Define “uninformative” as flat emissions after combining logB
    details["uninformative"] = details["flat"]

    summary = {
        "total_markers": int(T),
        "uninformative_total": int(details["uninformative"].sum()),
        "by_reason": {
            "variant_no_variants": int(no_variants.sum()),
            "variant_missing_genotype": int(missing_genotype.sum()),
            "variant_allele_not_in_reference": int(allele_not_in_reference.sum()),
            "hap_missing_genotype": int(hap_missing_genotype.sum()),
            "hap_missing_in_reference": int(hap_missing_in_reference.sum()),
            "missing_both": int(hap_and_all_variants_missing.sum()),
            "flat": int(flat.sum()),
        },
        "variant_missingness_pct_mean": float(np.nanmean(variant_missing_pct_nonmissing) if variant_missing_pct_nonmissing.size else np.nan),
        "variant_missingness_pct_median": float(np.nanmedian(variant_missing_pct_nonmissing) if variant_missing_pct_nonmissing.size else np.nan),
    }
    return summary, details



# -------------------------------
# Transition probabilities & HMM
# -------------------------------

def transmat_for_gap_unordered(
    d_bp,
    states,
    lam_per_Mb=HMM_DEFAULTS["lam_per_Mb"],
    same_chrom=True,
    trans_temp=HMM_DEFAULTS["trans_temp"],
    floor=1e-300,
):
    S = len(states)

    # Reset at chromosome boundary: no carry-over from previous state
    if not same_chrom:
        A = np.full((S, S), 1.0 / S, dtype=float)
    else:
        d_mb = max(float(d_bp), 0.0) / 1e6
        p_stay_h = np.exp(-lam_per_Mb * d_mb)
        p0 = p_stay_h**2
        p1_total = 2 * (1 - p_stay_h) * p_stay_h
        p2_total = (1 - p_stay_h)**2

        s2i = {s: i for i, s in enumerate(states)}
        clades = sorted({c for (i, j) in states for c in (i, j)})

        A = np.zeros((S, S), dtype=float)
        for k, (i, j) in enumerate(states):
            # stay
            A[k, k] += p0

            # one-haplotype switches
            single_targets = set()
            for x in clades:
                if x != i:
                    single_targets.add(tuple(sorted((x, j))))
                if x != j:
                    single_targets.add(tuple(sorted((i, x))))
            single_targets.discard((i, j))
            if single_targets:
                mass = p1_total / len(single_targets)
                for tgt in single_targets:
                    A[k, s2i[tgt]] += mass

            # two-haplotype switches
            double_targets = {
                tuple(sorted((x, y)))
                for x in clades for y in clades
                if (x, y) != (i, j) and x not in (i, j) and y not in (i, j)
            }
            if double_targets:
                mass2 = p2_total / len(double_targets)
                for tgt in double_targets:
                    A[k, s2i[tgt]] += mass2

        # normalize row-stochastic
        A /= A.sum(axis=1, keepdims=True)

    # ---- Tempering 
    # Temperature in log-space: divide logs by T, then renormalize rows.
    if trans_temp and trans_temp != 1.0:
        logA = np.log(np.clip(A, floor, 1.0))
        logA = logA / float(trans_temp)
        # row-wise log-softmax
        m = np.max(logA, axis=1, keepdims=True)
        A = np.exp(logA - (np.log(np.sum(np.exp(logA - m), axis=1, keepdims=True)) + m))

    # final numeric hygiene
    A = np.clip(A, floor, 1.0)
    A /= A.sum(axis=1, keepdims=True)
    return A


def forward_backward_unordered(
    logB,
    marker_positions,
    states,
    lam_per_Mb=HMM_DEFAULTS["lam_per_Mb"],
    floor=1e-300,
    trans_temp=HMM_DEFAULTS["trans_temp"],
):
    chroms = marker_positions[:,0]
    positions_bp = marker_positions[:,1].astype(float)
    T, S = logB.shape
    alpha = np.zeros((T,S), dtype=float)
    beta  = np.zeros((T,S), dtype=float)
    c = np.zeros(T, dtype=float)

    emis0 = np.exp(logB[0]); emis0 = np.maximum(emis0, 1e-300)
    alpha[0] = emis0 / S
    s0 = alpha[0].sum()
    if not np.isfinite(s0) or s0 < 1e-300:
        alpha[0] = np.full(S, 1.0/S); c[0] = 1.0
    else:
        c[0] = 1.0 / s0; alpha[0] *= c[0]

    for t in range(1, T):
        same = (chroms[t] == chroms[t-1])
        gap = (positions_bp[t] - positions_bp[t-1]) if same else 0.0
        A = transmat_for_gap_unordered(gap, states, lam_per_Mb, same_chrom=same, trans_temp=trans_temp)
        emis = np.exp(logB[t]); emis = np.maximum(emis, 1e-300)
        alpha[t] = emis * (alpha[t-1] @ A)
        s = alpha[t].sum()
        if not np.isfinite(s) or s < 1e-300:
            alpha[t] = np.full(S, 1.0/S); c[t] = 1.0
        else:
            c[t] = 1.0 / s; alpha[t] *= c[t]

    beta[-1] = 1.0
    for t in range(T-2, -1, -1):
        same = (chroms[t+1] == chroms[t])
        gap = (positions_bp[t+1] - positions_bp[t]) if same else 0.0
        A = transmat_for_gap_unordered(gap, states, lam_per_Mb, same_chrom=same, trans_temp=trans_temp)
        emis_next = np.exp(logB[t+1]); emis_next = np.maximum(emis_next, 1e-300)
        tmp = (A @ (emis_next * beta[t+1]))
        if not np.isfinite(tmp).all():
            beta[t] = 1.0
        else:
            beta[t] = tmp
            s = beta[t].sum()
            beta[t] = beta[t] / max(s, floor)      # <- normalize instead of * c[t+1]

    gamma = alpha * beta
    gamma /= np.clip(gamma.sum(axis=1, keepdims=True),floor, None)
    loglik = -np.sum(np.log(c))
    return gamma, loglik

def posterior_guided_viterbi_unordered(
    logB,
    marker_positions,
    states,
    lam_per_Mb=HMM_DEFAULTS["lam_per_Mb"],
    gamma=None,
    tau=HMM_DEFAULTS["tau"],
    emiss_temperature=1.0,
    floor=1e-300,
    trans_temp=HMM_DEFAULTS["trans_temp"],
):
    if gamma is not None and tau > 0.0:
        g = np.clip(gamma, floor, 1.0)
        logG = np.log(g)
        logB = (1.0 - tau) * logB + tau * logG
        logB = logB - logB.max(axis=1, keepdims=True)
    if emiss_temperature and emiss_temperature != 1.0:
        logB = logB / float(emiss_temperature)

    chroms = marker_positions[:,0]
    positions_bp = marker_positions[:,1].astype(float)
    T, S = logB.shape
    delta = np.zeros((T,S), dtype=float)
    psi = np.zeros((T,S), dtype=int)

    delta[0] = (-np.log(S)) + logB[0]
    psi[0] = -1

    for t in range(1, T):
        same = (chroms[t] == chroms[t-1])
        gap = (positions_bp[t] - positions_bp[t-1]) if same else 0.0
        A = transmat_for_gap_unordered(gap, states, lam_per_Mb=lam_per_Mb, same_chrom=same, trans_temp=trans_temp)
        logA = np.log(np.maximum(A, floor))

        for j in range(S):
            prev = delta[t-1] + logA[:, j]
            psi[t,j] = int(np.argmax(prev))
            delta[t,j] = np.max(prev) + logB[t,j]

    path_idx = np.zeros(T, dtype=int)
    path_idx[-1] = int(np.argmax(delta[-1]))
    for t in range(T-2, -1, -1):
        path_idx[t] = psi[t+1, path_idx[t+1]]

    path_states = [states[k] for k in path_idx]
    path_logprob = float(np.max(delta[-1]))
    return path_states, path_logprob



import numpy as np

def orient_min_switch(diploid_path, marker_positions):
    chroms = marker_positions[:, 0]
    T = len(diploid_path)

    # force 64-bit ints for DP to avoid overflow
    dp   = np.zeros((T, 2), dtype=np.int64)
    back = np.zeros((T, 2), dtype=np.int8)

    # Precompute same-chrom flags
    same_chrom = np.zeros(T, dtype=bool)
    same_chrom[1:] = chroms[1:] == chroms[:-1]

    # Big but safe sentinel in int64 space
    BIG = np.int64(2**62)

    # t=0 already zeros
    for t in range(1, T):
        i_prev, j_prev = diploid_path[t-1]
        i_cur,  j_cur  = diploid_path[t]

        # keep everything as strings; only the cost is integers
        for cur in (0, 1):  # 0: (i_cur,j_cur), 1: (j_cur,i_cur)
            h1c, h2c = (i_cur, j_cur) if cur == 0 else (j_cur, i_cur)
            best_val = BIG
            best_prev = 0
            for prev in (0, 1):
                h1p, h2p = (i_prev, j_prev) if prev == 0 else (j_prev, i_prev)
                c = (h1p != h1c) + (h2p != h2c) if same_chrom[t] else 0
                v = dp[t-1, prev] + np.int64(c)
                if v < best_val:
                    best_val = v
                    best_prev = prev
            dp[t, cur]   = best_val
            back[t, cur] = best_prev

    # Backtrack
    o = np.zeros(T, dtype=np.int8)
    o[-1] = 0 if dp[-1, 0] <= dp[-1, 1] else 1
    for t in range(T-2, -1, -1):
        o[t] = back[t+1, int(o[t+1])]

    # Build oriented haplotype label streams (still clade names)
    hap1, hap2 = [], []
    for t, (i, j) in enumerate(diploid_path):
        if o[t] == 0:
            hap1.append(i); hap2.append(j)
        else:
            hap1.append(j); hap2.append(i)
    return hap1, hap2



def orient_min_switch_v0(diploid_path, marker_positions):
    chroms = marker_positions[:,0]
    T = len(diploid_path)
    dp = np.zeros((T,2), dtype=np.int64)
    back = np.zeros((T,2), dtype=np.int8)
    for t in range(1, T):
        i_prev, j_prev = diploid_path[t-1]
        i_cur,  j_cur  = diploid_path[t]
        same_chrom = (chroms[t] == chroms[t-1])
        for cur in (0,1):
            h1_cur, h2_cur = (i_cur, j_cur) if cur==0 else (j_cur, i_cur)
            best = (10**9, 0)
            for prev in (0,1):
                h1_prev, h2_prev = (i_prev, j_prev) if prev==0 else (j_prev, i_prev)
                cost = (h1_prev != h1_cur) + (h2_prev != h2_cur) if same_chrom else 0
                val = dp[t-1, prev] + cost
                if val < best[0]:
                    best = (val, prev)
            dp[t,cur], back[t,cur] = best
    o = np.zeros(T, dtype=int)
    o[-1] = 0 if dp[-1,0] <= dp[-1,1] else 1
    for t in range(T-2, -1, -1):
        o[t] = back[t+1, o[t+1]]
    hap1, hap2 = [], []
    for t, (i,j) in enumerate(diploid_path):
        if o[t]==0:
            hap1.append(i); hap2.append(j)
        else:
            hap1.append(j); hap2.append(i)
    return hap1, hap2

def soften_posterior(gamma, T=3, floor=1e-300):
    g = np.clip(gamma, floor, 1.0)
    gT = g ** (1.0 / T)
    return gT / gT.sum(axis=1, keepdims=True)

# ---------------------------
# Driver
# ---------------------------

@dataclass
class RunParams:
    lam_per_Mb: float = HMM_DEFAULTS["lam_per_Mb"]
    strict_boost: float = HMM_DEFAULTS["strict_boost"]
    cap_total_boost_per_marker: Optional[float] = HMM_DEFAULTS["cap_total_boost_per_marker"]
    trans_temp: float = HMM_DEFAULTS["trans_temp"]
    alpha: float = HMM_DEFAULTS["alpha"]          # smoothing blend; 0 disables smoothing
    windowsize: int = HMM_DEFAULTS["windowsize"]         # smoothing window; 0 disables smoothing
    e_geno: float = HMM_DEFAULTS["e_geno"]
    e_homo: float = HMM_DEFAULTS["e_homo"]
    b0: float = HMM_DEFAULTS["b0"]
    certainty_softener: float = HMM_DEFAULTS["certainty_softener"]
    hom_soften_delta: float = HMM_DEFAULTS["hom_soften_delta"]
    hom_soften_width: float = HMM_DEFAULTS["hom_soften_width"]
    hom_min_mix: float = HMM_DEFAULTS["hom_min_mix"]
    hom_neutral: float = HMM_DEFAULTS["hom_neutral"]
    use_min_pos_for_marker: bool = True
    tau: float = HMM_DEFAULTS["tau"]            # posterior-guided Viterbi blend

def smooth_logB_by_context(logB, marker_df, halfwin=3, alpha=0.5, floor=-20.0):
    T,S = logB.shape
    if len(marker_df) != T:
        raise ValueError("marker_df length must match logB rows.")
    if "Chrom_numeric" not in marker_df.columns:
        raise KeyError("marker_df must have Chrom_numeric.")
    chroms = marker_df["Chrom_numeric"].to_numpy()

    # clean inputs
    X = np.array(logB, copy=True)
    X[~np.isfinite(X)] = floor

    # chromosome breakpoints
    cuts = np.r_[0, 1 + np.where(chroms[1:] != chroms[:-1])[0], T]

    sm = np.empty_like(X)
    base_hw = max(0, int(halfwin))
    
    for a,b in zip(cuts[:-1], cuts[1:]):
        seg = X[a:b]; L = b - a
        n = seg.shape[0]
        if n == 0:
            continue
            
        # Ensure effective half-window does not exceed segment size
        hw_eff = min(base_hw, max(0, (n - 1) // 2))
        if hw_eff == 0:
            # no smoothing possible; copy through
            sm[a:b, :] = seg
            continue
            
        # build Gaussian-like kernel with odd length <= n
        x = np.arange(-hw_eff, hw_eff + 1, dtype=float)
        denom = max(hw_eff, 1e-9)
        k = np.exp(-0.5 * (x / denom) ** 2)
        k /= k.sum()

        # convolve each state; 'same' returns length n
        for s in range(S):
            sm[a:b, s] = np.convolve(seg[:,s], k, mode="same")

    #blend with original
    blended = (1 - alpha) * X + alpha * sm
    row_max = blended.max(axis=1, keepdims=True)
    blended = np.maximum(blended - row_max, floor)
    return blended







def run_unordered_hmm_from_vcf(
    variants: Optional[List[VCFVariant]] = None,
    marker_df: Optional[pd.DataFrame] = None,
    clades: Optional[List[str]] = None,
    states: Optional[List[Tuple[str, str]]] = None,
    prof_df: Optional[pd.DataFrame] = None,
    freq_map: Optional[Dict[Tuple[str, int, str], np.ndarray]] = None,
    haplotypes: Optional[pd.Series] = None,
    hap_inf: Optional[pd.DataFrame] = None,
    freq_lookup: Optional[Dict[Tuple[str, str], np.ndarray]] = None,
    strict_alleles: Optional[Dict[Tuple[str, str], str]] = None,
    mus_hap_alleles: Optional[pd.DataFrame] = None,
    NONMus_hap_alleles: Optional[pd.DataFrame] = None,
    sample: Optional[str]=None,
    lam_per_Mb: float=HMM_DEFAULTS["lam_per_Mb"],
    strict_boost=HMM_DEFAULTS["strict_boost"],
    cap_total_boost_per_marker=HMM_DEFAULTS["cap_total_boost_per_marker"],
    trans_temp: float=HMM_DEFAULTS["trans_temp"],
    alpha: float=HMM_DEFAULTS["alpha"],
    windowsize: int=HMM_DEFAULTS["windowsize"],
    e_geno: float=HMM_DEFAULTS["e_geno"],
    e_homo: float=HMM_DEFAULTS["e_homo"],
    b0: float=HMM_DEFAULTS["b0"], #base preference for haplotype IDs (+) or variants (-)
    certainty_softener: float=HMM_DEFAULTS["certainty_softener"],
    hom_soften_delta: float=HMM_DEFAULTS["hom_soften_delta"],
    hom_soften_width: float=HMM_DEFAULTS["hom_soften_width"],
    hom_min_mix: float=HMM_DEFAULTS["hom_min_mix"],
    hom_neutral: float=HMM_DEFAULTS["hom_neutral"],
    tau: float=HMM_DEFAULTS["tau"], 
    verbose: bool=False,
    debug: bool=False):

    variants = list(variants or [])
    freq_map = freq_map or {}
    if marker_df is None:
        if variants:
            _, marker_df = build_variant_pos_df(variants)
        else:
            raise ValueError("marker_df or marker_positions_path is required when no VCF variants are available.")
    marker_df = normalize_marker_positions(marker_df)
    clades = _validate_clades(clades or DEFAULT_CLADES)
    states = states or states_from_clades(clades)

    has_variant_input = prof_df is not None and freq_map is not None
    has_hap_input = haplotypes is not None and freq_lookup is not None
    if not has_variant_input and not has_hap_input:
        raise ValueError("Provide VCF variant inputs, hap_genotype inputs, or both.")
    
    # 1) Extract positions/order
    marker_order = marker_df["marker"].tolist()
    marker_positions = marker_df[["Chrom_numeric","pos"]].to_numpy(dtype=int)

    if debug == True:
        print("marker_positions shape: ", marker_positions.shape)

    # 2) Check for chromosome 20 only when Mus chr20 diagnostic haplotype alleles are provided
    chr20_adjustments_enabled = (
        has_hap_input
        and mus_hap_alleles is not None
        and {"Marker", "Allele"}.issubset(mus_hap_alleles.columns)
    )
    if chr20_adjustments_enabled:
        chr20presence, homozygous, mus_matches = (
            check_for_chr_20(mus_hap_alleles, NONMus_hap_alleles, haplotypes) )
    else:
        chr20presence, homozygous, mus_matches = False, False, pd.DataFrame()

    if verbose == True and chr20_adjustments_enabled:
        print("chr20presence:", chr20presence)
        print("homozygous:", homozygous)
        if homozygous == "hemizygous":
            print("possible chr.7/chr.20 hemizygosity")


    # 3) Calculate emissions from the available input type(s)
    if verbose == True:
        print("Building emissions")

    logB_fromvariants = None
    logB_fromhaps = None
    if has_variant_input:
        logB_fromvariants = build_logB_unordered_from_variants_v4(
            variants=variants, freq_map=freq_map, clades=clades,
            states=states, marker_order=marker_order, prof_df=prof_df,
            e_geno=e_geno, e_homo=e_homo, hom_soften_delta=hom_soften_delta,
            hom_soften_width=hom_soften_width,
            hom_min_mix=hom_min_mix, hom_neutral=hom_neutral,
            debug=debug)

    if has_hap_input:
        logB_fromhaps = build_logB_from_hapgeno(test_series = haplotypes,
            strict_boost = strict_boost, strict_alleles=strict_alleles,
            freq_lookup=freq_lookup, clades=clades, states=states,
            marker_order=marker_order, eps=1e-9, e_geno=e_geno, e_homo=e_homo,
            hom_soften_delta=hom_soften_delta, hom_soften_width=hom_soften_width,
            hom_min_mix=hom_min_mix, hom_neutral=hom_neutral,
            cap_total_boost_per_marker=cap_total_boost_per_marker, debug=debug)

    if debug == True:
        if logB_fromvariants is not None:
            print("variant-based logB shape:", logB_fromvariants.shape)
        if logB_fromhaps is not None:
            print("haplotype-based logB shape:", logB_fromhaps.shape)

    # Combine haplotype and variant-based emissions when both are available.
    # Otherwise, use the single available source directly.
    if has_hap_input and has_variant_input:
        if hap_inf is None:
            raise ValueError("hap_inf/hap_informativeness_path is required when combining hap_genotype and VCF inputs.")
        hap_known_ref = hap_known_from_reference(test_series = haplotypes,
                    freq_lookup=freq_lookup, marker_order=marker_order)

        hap_scores = hap_strength_metrics_sample(
            haplotypes=haplotypes, hap_inf=hap_inf,
            marker_order=marker_order, clades=clades)
        if debug == True:
            print("hap informativeness scores calculated")

        var_scores = sum_In_per_marker_sample(variants=variants, marker_order=marker_order, prof_df=prof_df)
        if debug == True:
            print("variant informativeness scores calculated")

        logB_combined, logB_weights = combine_logB_mixture(
            logB_fromhaps, logB_fromvariants,
            hap_known_ref=hap_known_ref,
            hap_score_norm=hap_scores,
            var_score_norm=var_scores,
            b0=b0, debug=debug, verbose=verbose)
        if debug == True:
            print("logBs combined")
    elif has_hap_input:
        logB_combined = _row_norm(logB_fromhaps)
        logB_weights = np.ones(len(marker_order), dtype=float)
        if verbose:
            print("Using hap_genotype emissions only")
    else:
        logB_combined = _row_norm(logB_fromvariants)
        logB_weights = np.zeros(len(marker_order), dtype=float)
        if verbose:
            print("Using VCF emissions only")


    # 4) Check marker missingness/informativeness
    summary, details = uninformative_marker_breakdown(
        logB=logB_combined, marker_order=marker_order, 
        variants=variants, freq_map=freq_map,
        hap_test_series=haplotypes if has_hap_input else None,
        hap_freq_lookup=freq_lookup if has_hap_input else None,
        entropy_frac_thresh=0.95)

    if has_hap_input and has_variant_input:
        missing_markers = details["missing_both"]
    elif has_hap_input:
        missing_markers = details["hap_missing_genotype"]
    else:
        missing_markers = details["variant_missing_genotype"] | details["variant_no_variants"]

    if verbose==True or debug==True:
        print(f"Uninformative: {summary['uninformative_total']}/{summary['total_markers']}")
        if has_hap_input:
            print(
            f"\t missing hap genotype: {summary['by_reason']['hap_missing_genotype']}\n"
            f"\t haplotype ID not in reference set: {summary['by_reason']['hap_missing_in_reference']}"
            )
        if has_variant_input:
            print(
            f"\t markers with all missing variants: {summary['by_reason']['variant_missing_genotype']}"
            )
        if has_hap_input and has_variant_input:
            print(
            f"\t markers missing both hap allele ID and VCF variants: {summary['by_reason']['missing_both']}"
            )

    # Print warning if many markers are missing
    fraction_uninformative =  summary['uninformative_total'] / summary['total_markers']
    fraction_missing = float(np.mean(np.asarray(missing_markers, dtype=bool)))

    if fraction_missing > 0.5:
        if verbose ==True:
            print("Warning: Many missing markers. Results may be unreliable.")
        message = f"Warning: {round(fraction_missing * 100)}% markers missing. Results may be unreliable."
    elif summary['total_markers'] < 1000:
        if verbose ==True:
            print("Warning: Few markers available. Results may be unreliable.")
        message = f"Warning: Only {summary['total_markers']} markers available. Results may be unreliable."
    elif fraction_uninformative > 0.7:
        if verbose ==True:
            print("Warning: Many uninformative markers. Results may be unreliable.")
        message = f"Warning: {round(fraction_uninformative * 100)}% of markers are uninformative. Results may be unreliable."
    else: 
        message = None
    

    # 5) Smooth
    if verbose == True:
        print("Smoothing logB")
    if alpha and alpha > 0 and windowsize and windowsize > 1:
        logB_ctx = smooth_logB_by_context(logB_combined, marker_df, halfwin=windowsize//2, alpha=alpha)
    else:
        logB_ctx = logB_combined


    # 6) Forward-backward
    if verbose == True:
        print("Starting forward-backward algorithm")
    gamma, loglik = forward_backward_unordered(
        logB_ctx, marker_positions=marker_positions,
        states=states, lam_per_Mb=lam_per_Mb, trans_temp=trans_temp)

    # Certainty softening step:
    if certainty_softener and certainty_softener != 1.0:
        gamma = soften_posterior(gamma, T=certainty_softener)


    # 7) Viterbi (posterior-guided if tau>0)
    if verbose == True:
        print("Running Viterbi")
    diploid_path, v_ll = posterior_guided_viterbi_unordered(
        logB=logB_ctx, marker_positions=marker_positions, states=states,
        lam_per_Mb=lam_per_Mb, gamma=gamma, tau=tau, trans_temp=trans_temp)

    s2i = {s:k for k,s in enumerate(states)}
    path_idx = np.array([s2i[s] for s in diploid_path], dtype=int)
    dipl_post = gamma[np.arange(gamma.shape[0]), path_idx].tolist()

    
    # 8) Calculate margin certainty for graph
    def margin_certainty_called(gamma, path_idx, floor=1e-300):
        g = np.clip(gamma, floor, 1.0)
        T, S = g.shape
        called_p = g[np.arange(T), path_idx]
        return np.clip((called_p - 1.0/S) / (1.0 - 1.0/S), 0.0, 1.0)
    dipl_cert = margin_certainty_called(gamma, path_idx).tolist()

    hap1, hap2 = orient_min_switch(diploid_path, marker_positions)

    # 9) Generate final dataframes
    if verbose == True:
        print("Generating results")
    posteriors_df = pd.DataFrame(gamma, index=marker_order, columns=[f'{{{i},{j}}}' for (i,j) in states]).reset_index().rename(columns={'index':'Marker'})
    calls_unordered = [f'{{{i},{j}}}' for (i,j) in diploid_path]
    per_marker_df = pd.DataFrame({
        "Marker": marker_order,
        "UnorderedCall": calls_unordered,
        "DiplPostProb": dipl_post,
        "DiplCertainty": dipl_cert,
        "Hap1Clade": hap1,
        "Hap2Clade": hap2,
        "Chrom": marker_df["chrom"].tolist(),
        "Pos": marker_df["pos"].tolist(),
        "MissingMarkers": missing_markers
    })

    return {
        "loglik": float(v_ll),
        "posteriors_df": posteriors_df,
        "per_marker_df": per_marker_df,
        "logB_combined": logB_combined,
        "logB_weights": logB_weights,
        "chr20_adjustments_enabled": chr20_adjustments_enabled,
        "chr20presence": chr20presence,
        "chr20zygosity": homozygous,
        #"uninformative_details": details,
        #"missing_summary": summary,
        "missingness_rate": round(fraction_missing,2),
        "WarningMessage": message }






def _safe_sample_name(sample: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", str(sample))


def _make_output_dirs(outdir: str, make_plots: bool = False) -> Dict[str, Path]:
    out = Path(outdir)
    paths = {"outdir": out, "posteriors": out / "posteriors"}
    paths["posteriors"].mkdir(parents=True, exist_ok=True)
    if make_plots:
        paths["figures"] = out / "figures"
        paths["figures"].mkdir(parents=True, exist_ok=True)
    else:
        out.mkdir(parents=True, exist_ok=True)
    return paths


def _build_params_dict(kwargs: Dict[str, Any], outdir: str, clades: Sequence[str]) -> Dict[str, Any]:
    params = dict(kwargs)
    params["outdir"] = str(outdir)
    params["Kclust"] = len(clades)
    params["clades"] = list(clades)
    params.setdefault("refv", "NA")
    params.setdefault("modelv", "4")
    params.setdefault("vcf_refv", "NA")
    return params


def _write_run_log(outdir: str, params: Dict[str, Any], inputs: Dict[str, Any]) -> Path:
    timestamp = time.strftime("%Y-%m-%d_%H-%M-%S")
    logfile = Path(outdir) / f"run_{timestamp}.log"
    payload = {
        "timestamp": timestamp,
        "params": params,
        "inputs": {k: str(v) if v is not None else None for k, v in inputs.items()},
    }
    logfile.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    return logfile


def run_lai_hmm_for_sample(
    sample: str,
    *,
    marker_df: pd.DataFrame,
    clades: Sequence[str],
    states: Optional[List[Tuple[str, str]]] = None,
    variants: Optional[List[VCFVariant]] = None,
    haplotypes: Optional[pd.Series] = None,
    prof_df: Optional[pd.DataFrame] = None,
    freq_map: Optional[Dict[Tuple[str, int, str], np.ndarray]] = None,
    hap_inf: Optional[pd.DataFrame] = None,
    freq_lookup: Optional[Dict[Tuple[str, str], np.ndarray]] = None,
    strict_alleles: Optional[Dict[Tuple[str, str], str]] = None,
    mus_hap_alleles: Optional[pd.DataFrame] = None,
    nonmus_hap_alleles: Optional[pd.DataFrame] = None,
    outdir: Optional[str] = None,
    params: Optional[Dict[str, Any]] = None,
    chromlengths: Optional[pd.DataFrame] = None,
    save_outputs: bool = True,
    make_plots: bool = True,
    show_plots: bool = False,
    output_date: Optional[str] = None,
    verbose: bool = True,
    debug: bool = False,
    **hmm_kwargs
) -> Dict[str, Any]:
    """
    Run one sample from already loaded inputs.

    Returns a dictionary containing the raw HMM result, clade percentages, and any
    output file paths that were written.
    """
    output_date = output_date or time.strftime("%Y-%m-%d")
    states = states or states_from_clades(list(clades))
    params = dict(params or _build_params_dict(hmm_kwargs, outdir or ".", clades))
    params.setdefault("clades", list(clades))
    hmm_results = run_unordered_hmm_from_vcf(
        sample=sample,
        variants=variants,
        marker_df=marker_df,
        clades=list(clades),
        states=states,
        prof_df=prof_df,
        freq_map=freq_map,
        strict_alleles=strict_alleles,
        hap_inf=hap_inf,
        haplotypes=haplotypes,
        freq_lookup=freq_lookup,
        mus_hap_alleles=mus_hap_alleles,
        NONMus_hap_alleles=nonmus_hap_alleles,
        verbose=verbose,
        debug=debug,
        **hmm_kwargs,
    )
    params["chr20_adjustments_enabled"] = bool(hmm_results.get("chr20_adjustments_enabled"))
    params["chr20presence"] = bool(hmm_results.get("chr20presence"))

    per_marker_results_df = pd.merge(
        hmm_results["per_marker_df"],
        hmm_results["posteriors_df"],
        on="Marker",
        how="left",
    )
    hap_long, clade_percentages = prepare_to_plot(hmm_results)
    clade_row = {c: float(clade_percentages.get(c, 0.0)) for c in clades}
    clade_row["missingness"] = float(hmm_results["missingness_rate"])
    clade_row["warning"] = hmm_results.get("WarningMessage")

    written: Dict[str, str] = {}
    if save_outputs:
        if outdir is None:
            raise ValueError("outdir is required when save_outputs=True")
        _make_output_dirs(outdir, make_plots=make_plots)
        safe = _safe_sample_name(sample)
        out_csv = Path(outdir) / "posteriors" / f"{safe}_{output_date}.csv"
        per_marker_results_df.to_csv(out_csv, index=False)
        written["per_marker_posteriors"] = str(out_csv)

    if make_plots:
        if chromlengths is None:
            raise ValueError("chrom_lengths_path/chromlengths is required when make_plots=True")
        chromosome_painting(
            hap_long,
            clade_percentages,
            sample,
            params,
            chromlengths,
            message=hmm_results.get("WarningMessage"),
            showplot=show_plots,
            savefig=save_outputs,
        )

    return {
        "sample": sample,
        "status": "ok",
        "hmm_result": hmm_results,
        "per_marker_posteriors_df": per_marker_results_df,
        "clade_percentages": clade_row,
        "written": written,
    }


def run_lai_hmm(
    *,
    vcf_path: Optional[str] = None,
    hap_genotype_path: Optional[str] = None,
    marker_positions_path: Optional[str] = None,
    variant_profiles_path: Optional[str] = None,
    hap_freq_lookup_path: Optional[str] = None,
    reference_membership_path: Optional[str] = None,
    build_reference: bool = False,
    reference_outdir: Optional[str] = None,
    strict_alleles_path: Optional[str] = None,
    hap_informativeness_path: Optional[str] = None,
    mus_hap_alleles_path: Optional[str] = None,
    nonmus_hap_alleles_path: Optional[str] = None,
    chrom_lengths_path: Optional[str] = None,
    outdir: str = "results",
    samples: Optional[Sequence[str]] = None,
    sample: Optional[Any] = None,
    samples_file: Optional[str] = None,
    all_samples: bool = False,
    sample_source: str = "auto",
    threads: int = 1,
    clades: Sequence[str] = DEFAULT_CLADES,
    contigs: Optional[Set[str]] = None,
    require_marker: Optional[bool] = None,
    make_plots: bool = True,
    show_plots: bool = False,
    save_outputs: bool = True,
    write_log: bool = True,
    reference_pca: Optional[str] = None,
    reference_pca_downsample: float = 1.0,
    verbose: bool = True,
    debug: bool = False,
    **hmm_kwargs
) -> Dict[str, Any]:
    """
    High-level wrapper for command-line and Python use.

    Provide vcf_path, hap_genotype_path, or both. When both are present and a
    sample exists in both inputs, the emissions are combined; if a requested
    sample exists in only one input, that available source is used.
    """
    if not vcf_path and not hap_genotype_path:
        raise ValueError("Provide at least one genotype input: vcf_path and/or hap_genotype_path.")
    clades = _validate_clades(clades)
    effective_require_marker = resolve_require_marker(require_marker, vcf_path, hap_genotype_path)
    if "kappa" in hmm_kwargs:
        raise TypeError("The kappa parameter has been removed; please rerun without kappa.")
    reference_paths = None
    need_reference = (
        build_reference
        or (vcf_path and not variant_profiles_path)
        or (hap_genotype_path and not hap_freq_lookup_path)
        or (hap_genotype_path and vcf_path and not hap_informativeness_path)
    )
    if need_reference:
        if not reference_membership_path:
            missing = []
            if vcf_path and not variant_profiles_path:
                missing.append("variant_profiles_path")
            if hap_genotype_path and not hap_freq_lookup_path:
                missing.append("hap_freq_lookup_path")
            if hap_genotype_path and vcf_path and not hap_informativeness_path:
                missing.append("hap_informativeness_path")
            raise ValueError(
                "Missing precomputed reference files "
                f"({', '.join(missing) or 'requested build_reference'}). "
                "Provide them directly or provide reference_membership_path to build them."
            )
        reference_paths = build_reference_files(
            outdir=reference_outdir or str(Path(outdir) / "reference"),
            clades=clades,
            membership_path=reference_membership_path,
            hap_genotype_path=hap_genotype_path if hap_genotype_path else None,
            vcf_path=vcf_path if vcf_path else None,
            contigs=contigs,
            require_marker=effective_require_marker,
            pca=reference_pca,
            downsample=reference_pca_downsample,
            verbose=verbose,
        )
        variant_profiles_path = variant_profiles_path or reference_paths.get("variant_profiles")
        hap_freq_lookup_path = (
            hap_freq_lookup_path
            or reference_paths.get("hap_frequencies")
            or reference_paths.get("hap_frequency_lookup")
        )
        hap_informativeness_path = hap_informativeness_path or reference_paths.get("hap_informativeness")

    if vcf_path and not variant_profiles_path:
        raise ValueError("variant_profiles_path is required when vcf_path is provided.")
    if hap_genotype_path and not hap_freq_lookup_path:
        raise ValueError("hap_freq_lookup_path is required when hap_genotype_path is provided.")

    if sample_source not in {"auto", "intersection", "union"}:
        raise ValueError("sample_source must be one of: auto, intersection, union")
    if threads < 1:
        raise ValueError("threads must be >= 1")
    for key, value in HMM_DEFAULTS.items():
        hmm_kwargs.setdefault(key, value)

    sample_values: List[str] = []
    if sample is not None:
        if isinstance(sample, str):
            sample_values.append(sample)
        else:
            sample_values.extend([str(s) for s in sample])
    if samples is not None:
        if isinstance(samples, str):
            sample_values.append(samples)
        else:
            sample_values.extend([str(s) for s in samples])
    requested_samples = parse_samples_arg(sample_values, samples_file)
    run_all_unrequested = all_samples or requested_samples is None
    output_date = time.strftime("%Y-%m-%d")
    states = states_from_clades(clades)
    out_paths = _make_output_dirs(outdir, make_plots=make_plots)

    params = _build_params_dict(
        {
            "switch_rate": hmm_kwargs.get("lam_per_Mb", HMM_DEFAULTS["lam_per_Mb"]),
            "strict_boost": hmm_kwargs.get("strict_boost", HMM_DEFAULTS["strict_boost"]),
            "trans_temp": hmm_kwargs.get("trans_temp", HMM_DEFAULTS["trans_temp"]),
            "certainty_softener": hmm_kwargs.get("certainty_softener", HMM_DEFAULTS["certainty_softener"]),
            "e_geno": hmm_kwargs.get("e_geno", HMM_DEFAULTS["e_geno"]),
            "e_homo": hmm_kwargs.get("e_homo", HMM_DEFAULTS["e_homo"]),
            "b0": hmm_kwargs.get("b0", HMM_DEFAULTS["b0"]),
            "alpha": hmm_kwargs.get("alpha", HMM_DEFAULTS["alpha"]),
            "windowsize": hmm_kwargs.get("windowsize", HMM_DEFAULTS["windowsize"]),
            "tau": hmm_kwargs.get("tau", HMM_DEFAULTS["tau"]),
            "hom_soften_delta": hmm_kwargs.get("hom_soften_delta", HMM_DEFAULTS["hom_soften_delta"]),
            "hom_soften_width": hmm_kwargs.get("hom_soften_width", HMM_DEFAULTS["hom_soften_width"]),
            "hom_min_mix": hmm_kwargs.get("hom_min_mix", HMM_DEFAULTS["hom_min_mix"]),
            "hom_neutral": hmm_kwargs.get("hom_neutral", HMM_DEFAULTS["hom_neutral"]),
            "verbose": verbose,
        },
        outdir,
        clades,
    )

    if verbose:
        print("Loading inputs")

    marker_df = load_marker_positions(marker_positions_path) if marker_positions_path else None
    prof_df = None
    freq_map = None
    chunk = None
    vcf_header_samples: List[str] = []
    hap_gt = None
    hap_samples: List[str] = []
    hap_inf = None
    freq_lookup = None
    strict_alleles = None
    mus_hap_alleles = None
    nonmus_hap_alleles = None
    chromlengths = None

    if vcf_path:
        vcf_header_samples = get_vcf_samples(vcf_path)
    if hap_genotype_path:
        hap_gt = load_hap_genotypes(hap_genotype_path)
        hap_samples = list(hap_gt.columns)

    targets = choose_samples(
        vcf_samples=vcf_header_samples,
        hap_samples=hap_samples,
        requested_samples=requested_samples,
        all_samples=run_all_unrequested,
        sample_source=sample_source,
    )
    if verbose:
        print(f"Running HMM on {len(targets)} sample(s)")

    if vcf_path:
        vcf_sample_set = set(vcf_header_samples)
        load_samples = [s for s in targets if s in vcf_sample_set]
        if load_samples:
            prof_df, freq_map = load_variant_profiles(variant_profiles_path, clades=clades)
            if "specific" not in prof_df.columns:
                prof_df = variant_allele_informativeness(prof_df, clades)
            chunk = load_vcf_chunk_matrix(
                vcf_path,
                contig_whitelist=contigs,
                sample_names=load_samples,
                require_marker=effective_require_marker,
                debug=debug,
            )
            if marker_df is None:
                marker_df = marker_df_from_chunk(chunk)

    if hap_genotype_path:
        freq_lookup = load_hap_frequency_lookup(hap_freq_lookup_path, clades=clades)
        if hap_informativeness_path:
            hap_inf = read_table(hap_informativeness_path)
        elif freq_lookup is not None:
            hap_inf = hap_allele_informativeness(freq_lookup, clades)
        if strict_alleles_path:
            strict_alleles = load_pickle(strict_alleles_path)

    if marker_df is None:
        raise ValueError("marker_positions_path is required for hap_genotype-only runs.")
    marker_df = normalize_marker_positions(marker_df)

    if mus_hap_alleles_path:
        mus_hap_alleles = read_table(mus_hap_alleles_path)
    if nonmus_hap_alleles_path:
        nonmus_hap_alleles = read_table(nonmus_hap_alleles_path)
    if chrom_lengths_path:
        chromlengths = load_chrom_lengths(chrom_lengths_path)

    if make_plots and chromlengths is None:
        raise ValueError("chrom_lengths_path is required when make_plots=True.")

    if write_log and save_outputs:
        logfile = _write_run_log(
            outdir,
            params,
            {
                "vcf_path": vcf_path,
                "hap_genotype_path": hap_genotype_path,
                "marker_positions_path": marker_positions_path,
                "variant_profiles_path": variant_profiles_path,
                "hap_freq_lookup_path": hap_freq_lookup_path,
                "reference_membership_path": reference_membership_path,
                "reference_outputs": reference_paths,
                "strict_alleles_path": strict_alleles_path,
                "hap_informativeness_path": hap_informativeness_path,
                "mus_hap_alleles_path": mus_hap_alleles_path,
                "nonmus_hap_alleles_path": nonmus_hap_alleles_path,
                "chrom_lengths_path": chrom_lengths_path,
            },
        )
    else:
        logfile = None

    def _inputs_for_sample(s: str) -> Tuple[List[VCFVariant], Optional[pd.Series]]:
        sample_variants: List[VCFVariant] = []
        sample_haps: Optional[pd.Series] = None
        if chunk is not None and s in chunk.sample_to_idx:
            sample_variants = chunk.get_sample_variants(s)
        if hap_gt is not None and s in hap_gt.columns:
            sample_haps = hap_gt[s]
        if chunk is None and sample_haps is None:
            raise ValueError(f"{s} is not present in the hap_genotype input.")
        if chunk is not None and s not in chunk.sample_to_idx and sample_haps is None:
            raise ValueError(f"{s} is not present in the loaded VCF samples or hap_genotype input.")
        return sample_variants, sample_haps

    def _run_target(s: str) -> Dict[str, Any]:
        try:
            sample_variants, sample_haps = _inputs_for_sample(s)
            has_vcf_source = chunk is not None and s in chunk.sample_to_idx
            has_hap_source = sample_haps is not None
            return run_lai_hmm_for_sample(
                s,
                marker_df=marker_df,
                clades=clades,
                states=states,
                variants=sample_variants if has_vcf_source else None,
                haplotypes=sample_haps if has_hap_source else None,
                prof_df=prof_df if has_vcf_source else None,
                freq_map=freq_map if has_vcf_source else None,
                hap_inf=hap_inf,
                freq_lookup=freq_lookup if has_hap_source else None,
                strict_alleles=strict_alleles if has_hap_source else None,
                mus_hap_alleles=mus_hap_alleles,
                nonmus_hap_alleles=nonmus_hap_alleles,
                outdir=outdir,
                params=params,
                chromlengths=chromlengths,
                save_outputs=save_outputs,
                make_plots=make_plots,
                show_plots=show_plots,
                output_date=output_date,
                verbose=verbose,
                debug=debug,
                **hmm_kwargs,
            )
        except Exception as e:
            return {
                "sample": s,
                "status": "error",
                "error": f"{e.__class__.__name__}: {e}",
            }

    results: List[Dict[str, Any]] = []
    if threads == 1 or len(targets) == 1:
        for i, s in enumerate(targets, start=1):
            if verbose:
                print(f"[{i}/{len(targets)}] {s}")
            results.append(_run_target(s))
    else:
        with ThreadPoolExecutor(max_workers=threads) as ex:
            future_to_sample = {ex.submit(_run_target, s): s for s in targets}
            for i, fut in enumerate(as_completed(future_to_sample), start=1):
                s = future_to_sample[fut]
                if verbose:
                    print(f"[{i}/{len(targets)}] finished {s}")
                results.append(fut.result())

    summary_rows = []
    for res in results:
        if res.get("status") == "ok":
            row = {"sample": res["sample"], **res["clade_percentages"]}
        else:
            row = {"sample": res["sample"], "error": res.get("error")}
        summary_rows.append(row)
    summary_df = pd.DataFrame(summary_rows)
    summary_csv = None
    if save_outputs:
        summary_csv = out_paths["outdir"] / f"clade_percentage_summary_{output_date}.csv"
        summary_df.to_csv(summary_csv, index=False)

    return {
        "samples": targets,
        "results": results,
        "summary_df": summary_df,
        "summary_csv": str(summary_csv) if summary_csv else None,
        "logfile": str(logfile) if logfile else None,
        "outdir": str(out_paths["outdir"]),
        "reference_outputs": reference_paths,
    }



def prepare_to_plot(hmm_result):
    
    clade_marker_pos = hmm_result['per_marker_df']

    # Pivot to long format
    clade_marker_pos_renamed = clade_marker_pos.rename(
        columns={'Hap1Clade': 'Clade1', 'Hap2Clade': 'Clade2'})
    
    hap_long = pd.wide_to_long(clade_marker_pos_renamed, 
                              stubnames=['Clade'], 
                              i=['Marker', 'Chrom', 'Pos', 'MissingMarkers'], 
                              j='Haplotype', 
                              suffix='\\d+').reset_index()
    
    # Rename Haplotype values to 'Hap1' and 'Hap2'
    hap_long['Haplotype'] = 'Hap' + hap_long['Haplotype'].astype(str)
    hap_long['Chrom'] = hap_long['Chrom'].astype(str)

    assert hap_long['MissingMarkers'].dtype == bool

    
    df = hap_long.copy()
    # Ensure proper dtypes
    df['Chrom'] = df['Chrom'].astype(str)
    df['Pos']   = pd.to_numeric(df['Pos'], errors='coerce')
    df = df.dropna(subset=['Pos'])


    # Work haplotype-by-chromosome, sorted by position
    df = df.sort_values(['Haplotype', 'Chrom', 'Pos'])

    # Compute next position and next chromosome within (Haplotype, Chrom)
    df['next_Pos']   = df.groupby(['Haplotype', 'Chrom'])['Pos'].shift(-1)
    df['interval_bp'] = df['next_Pos'] - df['Pos']

    # Keep only valid within-chromosome intervals (drop last marker per chrom and any non-positive gaps)
    intervals = df.dropna(subset=['next_Pos']).copy()
    intervals = intervals[intervals['interval_bp'] > 0]

    # All clades: aggregate length * certainty per clade
    intervals['len_x_cert'] = intervals['interval_bp'] * intervals['DiplPostProb']
    agg = intervals.groupby('Clade', as_index=False)['len_x_cert'].sum()

    # Total diploid covered length = sum over both haplotypes, all chromosomes
    total_len_weighted = intervals['len_x_cert'].sum()
    
    agg['Percent'] = round(100.0 * agg['len_x_cert'] / total_len_weighted, 1)
    
    clade_percentages = dict(zip(agg['Clade'], agg['Percent']))


    chr20_adjustments_enabled = bool(hmm_result.get("chr20_adjustments_enabled"))
    if chr20_adjustments_enabled and hmm_result["chr20presence"] == True:
        chr20_start = CHR20_START_BP

        chr20 = hap_long[
            (hap_long["Chrom"] == "chr07")
            & (hap_long["Pos"] > chr20_start)
            & (hap_long["Clade"] == "Mus") 
            #& (hap_long["Marker"].isin(mus_hap_alleles["Marker"]))
        ].copy()
        
        chr20["Chrom"] = "20"
        chr20["Pos"] = chr20["Pos"] - chr20_start

        hap_long_nochr20 = hap_long[ ~(
        (hap_long["Chrom"] == "chr07")
        & (hap_long["Pos"] > chr20_start)
        & (hap_long["Clade"] == "Mus") )
        ].copy()
        
        hap_long_new = pd.concat([hap_long_nochr20, chr20], axis=0)

        return hap_long_new, clade_percentages
        
    else:
        return hap_long, clade_percentages
    


## Define plotting function

def chromosome_painting(hap_long, clade_percentages, indv, params, chromlengths, markerticks=True, message=None, annotation=None, savefig=True, showplot=True):
    import numpy as np
    import pandas as pd
    import matplotlib.pyplot as plt
    import matplotlib.patches as mpatches
    import matplotlib.lines as mlines
    from matplotlib.patches import Rectangle
    from matplotlib.ticker import FuncFormatter

    from datetime import datetime
    today = datetime.now().strftime("%m-%d-%y")

    p = params

    # {chrom(str): length_bp(float)}
    chrom_len_map = {
    str(r.chrom): float(r.pos)
    for r in chromlengths[['chrom', 'pos']].itertuples(index=False) }

    clade_fullnames = {
    "EA": "East Asian",
    "Mus": "Muscadinia",
    "NA1": "N. American 1",
    "NA2": "N. American 2",
    "Vv": "Vitis vinifera"}

    # set dtypes 
    df = hap_long.copy()
    # Chrom as string label
    df['Chrom'] = df['Chrom'].astype(str)
    # Positions as exact basepairs (int64) to avoid float equality issues
    df['Pos'] = pd.to_numeric(df['Pos'], errors='coerce').astype('Int64')  # pandas nullable int
    df = df.dropna(subset=['Pos']).copy()
    df['Pos'] = df['Pos'].astype(np.int64)

    requested_clades = [str(c).strip() for c in p.get("clades", []) if str(c).strip()]
    observed_clades = list(dict.fromkeys(str(c).strip() for c in df['Clade'].dropna() if str(c).strip()))
    palette_clades = requested_clades or observed_clades
    clade_colors = get_clade_colors(palette_clades)
    missing_palette_clades = [c for c in observed_clades if c not in clade_colors]
    for clade, color in get_clade_colors(missing_palette_clades).items():
        clade_colors.setdefault(clade, color)
    chr20_adjustment_active = bool(p.get("chr20_adjustments_enabled")) and bool(p.get("chr20presence"))

    # Sort chromosomes numerically
    def _chrom_key(c):
        s = str(c).replace("chr", "").replace("CHR", "")
        try:    return (0, int(s))
        except: return (1, s)
    chromosomes = sorted(df['Chrom'].unique(), key=_chrom_key)

    # Get max chromosome lengths to define plot boundaries:
    chrom_max_lengths = []
    plottable = []
    for chrom in chromosomes:
        chrom_key = str(chrom)
        chrom_data = df.loc[df['Chrom'] == chrom_key]
        if chrom_data.empty:
            continue
        pos = np.unique(chrom_data['Pos'].to_numpy())
        if pos.size == 0:
            continue

        # Keep chr07 plotted to the chr20 split point only when Mus chr20 adjustment is active.
        if chr20_adjustment_active and (chrom == "chr07") and (max(pos) < CHR20_START_BP):
            L = float(CHR20_START_BP)
        else:
            L = max(float(chrom_len_map.get(chrom_key, float(pos[-1]))), max(pos))
        
        chrom_max_lengths.append(L)
        plottable.append(chrom_key)


    # Guard against empty data
    if (not plottable) or (len(hap_long["MissingMarkers"].value_counts())==1):
        fig, ax = plt.subplots(figsize=(8, 3))
        ax.axis("off")
        ax.text(0.5, 0.5, f"No markers to plot for {indv}", ha="center", va="center", fontsize=16)
        if savefig:
            fig.savefig(f"{p['outdir']}/figures/{indv}_{today}.png", dpi=300, bbox_inches="tight")
        if showplot:
            plt.show()
        plt.close(fig)
        return

    maxL = max(chrom_max_lengths)   # the global x-axis maximum

    fig, ax = plt.subplots(figsize=(12, max(3, len(chromosomes) * 0.9)))
    ax.set_xlim(0, maxL)

    # invert to make chrom 1 appear at the top
    ax.invert_yaxis()

    ax.set_facecolor("white")                    # white background
    ax.grid(False, which='both', axis='both')    # no x and y grids


    for i, chrom in enumerate(chromosomes):
        chrom_key = str(chrom)
        chrom_data = df.loc[df['Chrom'] == chrom_key].copy()

        # Sort per-haplotype by position
        chrom_data = chrom_data.sort_values(['Haplotype', 'Pos'], kind='mergesort')

        # Master position vector (sorted, unique) exactly from hap_long
        pos = np.unique(chrom_data['Pos'].to_numpy())
        if pos.size == 0:
            # nothing to draw for this chromosome
            continue

        # Chromosome length (bp)
        L = chrom_max_lengths[i]

        # Edges by midpoints between marker positions
        mids = (pos[:-1].astype(np.float64) + pos[1:].astype(np.float64)) / 2.0
        edges = np.empty(pos.size + 1, dtype=float)
        edges[0]    = 0.0
        edges[1:-1] = mids
        edges[-1]   = L

        left_bounds, right_bounds = edges[:-1], edges[1:]

        #put an opaque background band behind each haplotype
        lane_h = 0.18  # height of each hap "track" in data y-units

        # Draw both haplotype tracks using the same bounds
        for hap, offset in (('Hap1', +0.15), ('Hap2', -0.15)):
            h = chrom_data.loc[chrom_data['Haplotype'] == hap, ['Pos','Clade','DiplCertainty']].copy()
            h['Pos'] = pd.to_numeric(h['Pos'], errors='coerce').dropna().astype(np.int64)
            h['DiplCertainty'] = pd.to_numeric(h['DiplCertainty'], errors='coerce').fillna(0.0)
            
            # exact mapping from integer bp positions
            pos_to_clade = dict(zip(h['Pos'].to_numpy(), h['Clade'].to_numpy()))
            pos_to_cert  = dict(zip(h['Pos'].to_numpy(), h['DiplCertainty'].astype(float).to_numpy()))

            # rectangle geometry for this hap's lane
            lane_h = 0.2
            y0 = i + offset - lane_h/2

            for j, ppos in enumerate(pos):
                clade = pos_to_clade.get(ppos, None)
                raw_cert = pos_to_cert.get(ppos, 0.0)
                cert = float(raw_cert)
                if not np.isfinite(cert):
                    cert = 0.0

                # clamp alpha to [0,1] to avoid matplotlib warnings
                alpha = float(np.clip(cert, 0.0, 1.0))
                if alpha <= 0.0:
                    continue

                L = float(left_bounds[j]); R = float(right_bounds[j]); W = R - L
                if W <= 0:
                    continue


                ax.add_patch(Rectangle(
                    (L, y0), W, lane_h,
                    facecolor=clade_colors.get(clade, 'gray'), edgecolor='none',
                    linewidth=0, antialiased=False, zorder=3, alpha=alpha))

        
        if markerticks == True:
            # Get unique tick positions for this chromosome where markers are NOT missing
            ticks = (
                df.loc[(df['Chrom'] == chrom_key) & (~df['MissingMarkers']),
                       ['Marker', 'Pos']]
                  .drop_duplicates(subset=['Marker'])  # one tick per marker
                  .sort_values('Pos'))
            
            y0, y1 = i + 0.50, i + 0.35
            ax.vlines(ticks['Pos'].to_numpy(), y0, y1, color='k', linewidth=0.8, zorder=6)
                    

        # Chromosome label at left
        plt.text(-2_500_000, i, f'Chr {i+1}', va='center', ha='right', fontsize=20)
        


    # Aesthetics
    pad = 0.05 * maxL   # 2% of the chromosome length padding on either side
    ax.set_xlim(-pad, maxL + pad) # set x-axis limit to max chrom length + padding

    # Format x-axis ticks in Mbp (bp / 1e6)
    ax.xaxis.set_major_formatter(FuncFormatter(lambda x, pos: f"{x/1e6:.0f}"))
    ax.set_xlabel("Physical Position (Mbp)", fontsize=20)
    
    ax.set_yticks([])
    ax.tick_params(axis='x', labelsize=20)
    ax.set_title(f"Clade Assignment for {indv}", fontsize=20)
    fig.tight_layout()

    # Legend
    unique_clades = [c for c in requested_clades if c in observed_clades]
    unique_clades.extend([c for c in observed_clades if c not in unique_clades])
    legend_handles = [
        mpatches.Patch(
            color=clade_colors.get(clade, 'gray'),
            label=f"{clade_fullnames.get(clade, clade)} - {clade_percentages.get(clade, 0):.1f}%"
        )
        for clade in unique_clades
    ]
    if markerticks == True:
        legend_handles.append(
            mlines.Line2D([], [], color='k', linestyle='None', marker='|', markersize=18,
                          markeredgewidth=2, label='rhAmpSeq markers'))
    ax.legend(handles=legend_handles, title="Clade (certainty = opacity)",
              bbox_to_anchor=(0.97, 0.98), loc='upper center', fontsize=20, title_fontsize=20,
             framealpha=1)

    # Print warning message or other annotations on the plot
    if annotation != None and message != None:
        plt.suptitle(f"{annotation}\n {message}", y=0.97)
    if annotation != None and message == None:
        plt.suptitle(annotation, y=0.96, fontsize=14)
    if annotation == None and message != None:
        plt.suptitle(f"{message}", y=0.96, fontsize=14)

    # Print version number on the plot
    version_label = f"version: {p.get('refv', 'NA')}_v{p.get('modelv', 'NA')}_{p.get('vcf_refv', 'NA')}"
    fig.text(0.01, 0.01, # bottom left corner
             version_label, ha="left", va="bottom", fontsize=8)
    
    # Save
    if savefig == True:
        fig.savefig(f"{p['outdir']}/figures/{indv}_{today}.png", dpi=400, bbox_inches='tight')
            
#save with parameters:
#lmda{p['switch_rate']}_K{p['Kclust']}_Q{p['Qmax']}"    #f"_min{p['COMMON_FREQUENCY']}_max{p['MAX_FREQUENCY']}_sb{p['strict_boost']}.png"

    if showplot == True:
        plt.show()

    plt.close()






def build_arg_parser():
    import argparse

    parser = argparse.ArgumentParser(
        description="Run the LAI HMM from VCF variants, hap_genotype calls, or both.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  Build reference files only:
    python LAI_HMM_v0.4.py --step build-reference --hap-genotype example_data/hap_genotype_refset_example.gz --vcf example_data/example_refset.vcf.gz --reference-membership example_data/example_reference_membership.tsv --clades EA,Mus,NA,Vv --reference-outdir demo_reference --pca both --downsample 0.5

  VCF only:
    python LAI_HMM_v0.4.py --vcf example_data/example_samples.vcf.gz --variant-profiles example_data/reference_variant_profiles.tsv --marker-positions example_data/marker_positions.csv --no-plots --sample SAMPLE1

  hap_genotype only:
    python LAI_HMM_v0.4.py --hap-genotype example_data/hap_genotype_samples_example --hap-frequencies example_data/reference_hap_allele_frequency_lookup.pkl --marker-positions example_data/marker_positions.csv --chrom-lengths example_data/chrom_lengths.fai --all-samples

  combined VCF + hap_genotype, four worker threads:
    python LAI_HMM_v0.4.py --vcf example_data/example_samples.vcf.gz --hap-genotype example_data/hap_genotype_samples_example --variant-profiles example_data/reference_variant_profiles.tsv --hap-frequencies example_data/reference_hap_allele_frequency_lookup.pkl --hap-informativeness example_data/reference_hap_allele_informativeness.tsv --marker-positions example_data/marker_positions.csv --chrom-lengths example_data/chrom_lengths.fai --all-samples --threads 4
""",
    )

    inputs = parser.add_argument_group("inputs")
    inputs.add_argument("--vcf", help="Sample or cohort VCF/BCF. Indexed bgzipped VCF is recommended.")
    inputs.add_argument("--hap-genotype", "--hap_genotype", dest="hap_genotype", help="hap_genotype matrix with marker rows and sample columns.")
    inputs.add_argument("--marker-positions", "--marker_positions", dest="marker_positions", help="Marker position table with marker, chrom, pos columns.")
    inputs.add_argument("--variant-profiles", "--profiles", dest="variant_profiles", help="Clade-specific variant allele frequency/profile TSV.")
    inputs.add_argument("--hap-frequencies", "--hap_freq_lookup", dest="hap_freq_lookup", help="Clade-specific haplotype allele frequencies as either a pickle lookup or a delimited text table with marker, allele_id, and clade columns.")
    inputs.add_argument("--reference-membership", "--reference_membership", dest="reference_membership", help="Reference sample-to-clade table with sample/IID and clade/group columns.")
    inputs.add_argument("--reference-outdir", "--reference_outdir", dest="reference_outdir", help="Directory for generated reference files. Defaults to OUTDIR/reference.")
    inputs.add_argument("--step", choices=["all", "build-reference", "run-hmm", "variant-informativeness", "hap-informativeness"], default="all", help="Run the whole workflow or one reusable step.")
    inputs.add_argument("--strict-alleles", "--strict_alleles", dest="strict_alleles", help="Optional pickle of strict diagnostic haplotype alleles.")
    inputs.add_argument("--hap-informativeness", "--hap_inf", dest="hap_informativeness", help="Optional haplotype allele informativeness TSV. Required for combined VCF + hap_genotype weighting.")
    inputs.add_argument("--mus-hap-alleles", dest="mus_hap_alleles", help="Optional Mus chr20 diagnostic hap allele CSV/TSV.")
    inputs.add_argument("--nonmus-hap-alleles", dest="nonmus_hap_alleles", help="Optional non-Mus chr7 diagnostic hap allele CSV/TSV.")
    inputs.add_argument("--chrom-lengths", "--chrom_lengths", dest="chrom_lengths", help="Optional chromosome length table or FASTA .fai for plotting.")
    inputs.add_argument(
        "--pca",
        nargs="?",
        const="auto",
        choices=["auto", "hap", "vcf", "both"],
        default=None,
        help=(
            "During reference-building, run PCA on reference samples and save PCA coordinates, "
            "metrics, and plots. Use --pca alone to run on whichever of --hap-genotype and/or "
            "--vcf were provided, or specify one of: hap, vcf, both."
        ),
    )
    inputs.add_argument(
        "--downsample",
        type=float,
        default=1.0,
        help=(
            "During PCA for reference-building, randomly retain this proportion of feature columns "
            "before PCA. Default 1 uses the full matrix; 0.5 keeps about half of the features."
        ),
    )

    samples = parser.add_argument_group("samples and output")
    samples.add_argument("--sample", action="append", help="Sample name to run. Can be repeated or comma-separated.")
    samples.add_argument("--samples-file", "--samples_file", dest="samples_file", help="Text file with one sample name per line.")
    samples.add_argument("--all-samples", "--all_samples", dest="all_samples", action="store_true", help="Run all selected samples. This is the default when no --sample or --samples-file is provided.")
    samples.add_argument("--sample-source", choices=["auto", "intersection", "union"], default="auto", help="When both VCF and hap_genotype inputs are present, auto/intersection runs shared samples; union runs samples present in either input.")
    samples.add_argument("--outdir", default="results", help="Output directory.")
    samples.add_argument("--threads", type=int, default=1, help="Number of worker threads for per-sample parallel runs.")
    samples.add_argument("--plots", dest="plots", action="store_true", default=True, help="Save chromosome painting figures. Enabled by default and requires --chrom-lengths.")
    samples.add_argument("--no-plots", "--no_plots", dest="plots", action="store_false", help="Skip chromosome painting figures.")
    samples.add_argument("--show-plots", "--show_plots", dest="show_plots", action="store_true", help="Display plots interactively while saving them.")
    samples.add_argument("--no-save", dest="save_outputs", action="store_false", help="Run without writing CSV/plot/log outputs.")
    samples.set_defaults(save_outputs=True)

    model = parser.add_argument_group("model parameters")
    model.add_argument("--clades", default="EA,Mus,NA1,NA2,Vv", help="Comma-separated group/clade/population/species names in the order used by reference profiles.")
    model.add_argument("--contigs", help="Comma-separated VCF contigs to load, for example chr01,chr02.")
    model.add_argument("--allow-missing-marker", action="store_true", help="Legacy compatibility flag for VCF-only workflows. If the selected VCF records do not use INFO/MARKER, unlabeled variants fall back to CHROM:POS site labels; otherwise unlabeled records are skipped. Combined VCF + hap_genotype runs still require INFO/MARKER.")
    model.add_argument("--lam-per-Mb", "--lam_per_Mb", dest="lam_per_Mb", type=float, default=HMM_DEFAULTS["lam_per_Mb"])
    model.add_argument("--strict-boost", "--strict_boost", dest="strict_boost", type=float, default=HMM_DEFAULTS["strict_boost"])
    model.add_argument("--cap-total-boost-per-marker", "--cap_total_boost_per_marker", dest="cap_total_boost_per_marker", type=float, default=HMM_DEFAULTS["cap_total_boost_per_marker"])
    model.add_argument("--trans-temp", "--trans_temp", dest="trans_temp", type=float, default=HMM_DEFAULTS["trans_temp"])
    model.add_argument("--alpha", type=float, default=HMM_DEFAULTS["alpha"], help="Context smoothing blend. Default 0 disables smoothing.")
    model.add_argument("--windowsize", type=int, default=HMM_DEFAULTS["windowsize"], help="Context smoothing window. Default 0 disables smoothing.")
    model.add_argument("--e-geno", "--e_geno", dest="e_geno", type=float, default=HMM_DEFAULTS["e_geno"])
    model.add_argument("--e-homo", "--e_homo", dest="e_homo", type=float, default=HMM_DEFAULTS["e_homo"], help="Homozygote dropout softening rate; mixes in possible undercalled heterozygotes before informativeness-based homozygote flattening.")
    model.add_argument("--b0", type=float, default=HMM_DEFAULTS["b0"])
    model.add_argument("--certainty-softener", "--certainty_softener", dest="certainty_softener", type=float, default=HMM_DEFAULTS["certainty_softener"], help="Posterior softening temperature. Default 1 disables softening.")
    model.add_argument("--hom-soften-delta", "--hom_soften_delta", dest="hom_soften_delta", type=float, default=HMM_DEFAULTS["hom_soften_delta"])
    model.add_argument("--hom-soften-width", "--hom_soften_width", dest="hom_soften_width", type=float, default=HMM_DEFAULTS["hom_soften_width"])
    model.add_argument("--hom-min-mix", "--hom_min_mix", dest="hom_min_mix", type=float, default=HMM_DEFAULTS["hom_min_mix"], help="Minimum weight retained on the dropout-adjusted homozygote emission during informativeness-based softening.")
    model.add_argument("--hom-neutral", "--hom_neutral", dest="hom_neutral", type=float, default=HMM_DEFAULTS["hom_neutral"])
    model.add_argument("--tau", type=float, default=HMM_DEFAULTS["tau"])

    debug = parser.add_argument_group("logging")
    debug.add_argument("--verbose", dest="verbose", action="store_true", default=True, help="Print progress messages. Enabled by default.")
    debug.add_argument("--quiet", dest="verbose", action="store_false", help="Suppress progress messages.")
    debug.add_argument("--debug", action="store_true")
    return parser


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = build_arg_parser()
    args = parser.parse_args(argv)

    clades = _validate_clades([c.strip() for c in args.clades.split(",") if c.strip()])
    reference_outdir = args.reference_outdir or str(Path(args.outdir) / "reference")
    try:
        downsample = validate_downsample_proportion(args.downsample)
    except ValueError as e:
        parser.error(str(e))

    if args.step == "build-reference":
        if not args.reference_membership:
            parser.error("--reference-membership is required for --step build-reference")
        outputs = build_reference_files(
            outdir=reference_outdir,
            clades=clades,
            membership_path=args.reference_membership,
            hap_genotype_path=args.hap_genotype,
            vcf_path=args.vcf,
            verbose=args.verbose,
            pca=args.pca,
            downsample=downsample,
            contigs=parse_contigs(args.contigs),
            require_marker=False if args.allow_missing_marker else None,
        )
        print(json.dumps(outputs, indent=2, sort_keys=True))
        return 0

    if args.step == "variant-informativeness":
        if not args.variant_profiles:
            parser.error("--variant-profiles is required for --step variant-informativeness")
        df = read_table(args.variant_profiles)
        out = variant_allele_informativeness(df, clades)
        Path(args.outdir).mkdir(parents=True, exist_ok=True)
        out_path = Path(args.outdir) / "variant_profiles_with_informativeness.tsv"
        out.to_csv(out_path, sep="\t", index=False)
        locus = allele_informativeness_by_locus(df, clades)
        locus_path = Path(args.outdir) / "variant_locus_informativeness.tsv"
        locus.to_csv(locus_path, sep="\t", index=False)
        print(f"Wrote: {out_path}")
        print(f"Wrote: {locus_path}")
        return 0

    if args.step == "hap-informativeness":
        if not args.hap_freq_lookup:
            parser.error("--hap-frequencies is required for --step hap-informativeness")
        freq_lookup = load_hap_frequency_lookup(args.hap_freq_lookup, clades=clades)
        out = hap_allele_informativeness(freq_lookup, clades)
        Path(args.outdir).mkdir(parents=True, exist_ok=True)
        out_path = Path(args.outdir) / "hap_allele_informativeness.tsv"
        out.to_csv(out_path, sep="\t", index=False)
        print(f"Wrote: {out_path}")
        return 0

    result = run_lai_hmm(
        vcf_path=args.vcf,
        hap_genotype_path=args.hap_genotype,
        marker_positions_path=args.marker_positions,
        variant_profiles_path=args.variant_profiles,
        hap_freq_lookup_path=args.hap_freq_lookup,
        reference_membership_path=args.reference_membership,
        build_reference=(args.step == "all" and bool(args.reference_membership)),
        reference_outdir=args.reference_outdir,
        strict_alleles_path=args.strict_alleles,
        hap_informativeness_path=args.hap_informativeness,
        mus_hap_alleles_path=args.mus_hap_alleles,
        nonmus_hap_alleles_path=args.nonmus_hap_alleles,
        chrom_lengths_path=args.chrom_lengths,
        outdir=args.outdir,
        samples=args.sample,
        samples_file=args.samples_file,
        all_samples=args.all_samples,
        sample_source=args.sample_source,
        threads=args.threads,
        clades=clades,
        contigs=parse_contigs(args.contigs),
        require_marker=False if args.allow_missing_marker else None,
        make_plots=args.plots,
        show_plots=args.show_plots,
        save_outputs=args.save_outputs,
        reference_pca=args.pca,
        reference_pca_downsample=downsample,
        verbose=args.verbose,
        debug=args.debug,
        lam_per_Mb=args.lam_per_Mb,
        strict_boost=args.strict_boost,
        cap_total_boost_per_marker=args.cap_total_boost_per_marker,
        trans_temp=args.trans_temp,
        alpha=args.alpha,
        windowsize=args.windowsize,
        e_geno=args.e_geno,
        e_homo=args.e_homo,
        b0=args.b0,
        certainty_softener=args.certainty_softener,
        hom_soften_delta=args.hom_soften_delta,
        hom_soften_width=args.hom_soften_width,
        hom_min_mix=args.hom_min_mix,
        hom_neutral=args.hom_neutral,
        tau=args.tau,
    )

    errors = [r for r in result["results"] if r.get("status") != "ok"]
    if result.get("summary_csv"):
        print(f"Wrote summary: {result['summary_csv']}")
    print(f"Finished {len(result['results']) - len(errors)}/{len(result['results'])} sample(s)")
    if errors:
        for err in errors[:10]:
            print(f"ERROR {err['sample']}: {err.get('error')}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
