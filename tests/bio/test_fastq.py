"""
Tests for edms.bio.fastq

fastq.py is ~6200 lines, the majority of which are either:
  (a) thin wrappers around external CLI bioinformatics tools (umi_tools, cutadapt,
      bowtie2, samtools, fgbio) that live in a separate 'umi_tools' conda env not
      available here, or
  (b) matplotlib/seaborn plotting functions built on top of (c).
Real test effort is focused on (c): the pure-Python/pandas/Biopython parsing,
quality-score handling, sequence alignment, and data-wrangling functions that do the
actual scientific computation. Three real bugs were found and fixed directly in
fastq.py during this work (split_fastqs_by_expected_errors()'s unsupported
write_gzip kwarg, trim_filter()'s 5' trim off-by-one, and extract_umis()'s
.extract_umi/.extract_umis log-path typo); the tests below assert the corrected
behavior.

Covers:
- fuzzy_substring_search(): Levenshtein-distance substring search
- _iter_fastq_files() / _natural_sorted_paths(): file discovery & natural sorting
- revcom_fastqs() / unzip_fastqs() / comb_fastqs(): real file I/O via Bio.SeqIO/gzip
  (no external binaries involved)
- split_fastqs_by_expected_errors(): expected-error quality computation
- trim_filter(): quality-based read filtering, trimming, and N-masking
- find_AA_edits(), trim(), format_alignment(), find_indel(): sequence-alignment
  helpers behind genotype()
- genotype(), outcomes(), abundances(): edit-calling and tidy summarization
- edit_change(): amino-acid property/conservation classification
- _signature_phred_summary(): per-signature Phred quality summaries
- extract_umis() / trim_motifs(): subprocess command construction (mocked; the
  underlying binaries are not installed in this env)
"""
import gzip
import os

import numpy as np
import pandas as pd
import pytest
from Bio.Seq import Seq
from Bio.Align import PairwiseAligner
from Bio.SeqRecord import SeqRecord

from edms.bio import fastq
from edms.bio.signature import signature_from_alignment


# =========================================================================== #
# fuzzy_substring_search()
# =========================================================================== #
def test_fuzzy_substring_search_exact_match():
    df = fastq.fuzzy_substring_search(text="AAACGTAAA", pattern="CGT", max_distance=0)
    assert len(df) == 1
    row = df.iloc[0]
    assert row["window"] == "CGT"
    assert row["distance"] == 0
    assert row["start_i"] == 3
    assert row["end_i"] == 6


def test_fuzzy_substring_search_within_distance():
    df = fastq.fuzzy_substring_search(text="AAACGAAAA", pattern="CGT", max_distance=1)
    assert len(df) == 1
    assert df.iloc[0]["window"] == "CGA"
    assert df.iloc[0]["distance"] == 1


def test_fuzzy_substring_search_no_match_above_distance():
    df = fastq.fuzzy_substring_search(text="AAAAAAAAA", pattern="CGT", max_distance=0)
    assert len(df) == 0
    assert list(df.columns) == ["pattern", "window", "distance", "start_i", "end_i"]


def test_fuzzy_substring_search_multiple_hits():
    df = fastq.fuzzy_substring_search(text="CGTAACGT", pattern="CGT", max_distance=0)
    assert sorted(df["start_i"]) == [0, 5]


# =========================================================================== #
# _iter_fastq_files() / _natural_sorted_paths()
# =========================================================================== #
def test_iter_fastq_files_finds_fastq_and_gz_recursively(tmp_path):
    (tmp_path / "sub").mkdir()
    (tmp_path / "a.fastq").write_text("@x\nA\n+\nI\n")
    (tmp_path / "sub" / "b.fastq.gz").write_bytes(b"")
    (tmp_path / "not_fastq.txt").write_text("ignore me")

    found = sorted(fastq._iter_fastq_files(str(tmp_path)))
    basenames = sorted(os.path.basename(p) for p in found)
    assert basenames == ["a.fastq", "b.fastq.gz"]


def test_iter_fastq_files_empty_dir(tmp_path):
    assert list(fastq._iter_fastq_files(str(tmp_path))) == []


def test_natural_sorted_paths_basename_orders_numerically():
    paths = ["/x/file10.fastq", "/x/file2.fastq", "/x/file1.fastq"]
    result = fastq._natural_sorted_paths(paths)
    assert [os.path.basename(p) for p in result] == ["file1.fastq", "file2.fastq", "file10.fastq"]


def test_natural_sorted_paths_key_root_uses_relpath():
    root = "/root"
    paths = ["/root/sub/file2.fastq", "/root/sub/file10.fastq"]
    result = fastq._natural_sorted_paths(paths, key_root=root)
    assert result == ["/root/sub/file2.fastq", "/root/sub/file10.fastq"]


# =========================================================================== #
# revcom_fastqs() / unzip_fastqs() / comb_fastqs()
# (real Bio.SeqIO + gzip file I/O -- no external binaries)
# =========================================================================== #
FASTQ_CONTENT = "@r1\nACGT\n+\nIIII\n@r2\nTTTT\n+\nIIII\n"


@pytest.fixture
def fastq_in_dir(tmp_path):
    in_dir = tmp_path / "in"
    in_dir.mkdir()
    (in_dir / "a.fastq").write_text(FASTQ_CONTENT)
    with gzip.open(in_dir / "b.fastq.gz", "wt") as f:
        f.write(FASTQ_CONTENT)
    return in_dir


def test_revcom_fastqs_uncompressed(fastq_in_dir, tmp_path):
    out_dir = tmp_path / "revcom"
    fastq.revcom_fastqs(in_dir=str(fastq_in_dir), out_dir=str(out_dir))

    content = (out_dir / "a.fastq").read_text()
    lines = content.strip().split("\n")
    # r1 = ACGT (palindromic revcomp), r2 = TTTT -> revcomp AAAA
    assert lines[1] == "ACGT"
    assert lines[5] == "AAAA"


def test_revcom_fastqs_gzipped(fastq_in_dir, tmp_path):
    out_dir = tmp_path / "revcom"
    fastq.revcom_fastqs(in_dir=str(fastq_in_dir), out_dir=str(out_dir))
    with gzip.open(out_dir / "b.fastq.gz", "rt") as f:
        lines = f.read().strip().split("\n")
    assert lines[5] == "AAAA"


def test_unzip_fastqs_decompresses_only_gz_files(fastq_in_dir, tmp_path):
    out_dir = tmp_path / "unzip"
    fastq.unzip_fastqs(in_dir=str(fastq_in_dir), out_dir=str(out_dir))
    # Only the .fastq.gz input produces an unzipped output (unzip_fastqs ignores
    # already-uncompressed .fastq files under its own logic).
    assert os.listdir(out_dir) == ["b.fastq"]
    assert (out_dir / "b.fastq").read_text() == FASTQ_CONTENT


def test_comb_fastqs_non_recursive_concatenates_in_natural_order(fastq_in_dir, tmp_path):
    out_dir = tmp_path / "comb"
    fastq.comb_fastqs(in_dir=str(fastq_in_dir), out_dir=str(out_dir), out_file="combined.fastq")
    combined = (out_dir / "combined.fastq").read_text()
    assert combined == FASTQ_CONTENT + FASTQ_CONTENT  # a.fastq then b.fastq.gz (natural sort)


def test_comb_fastqs_gzipped_output(fastq_in_dir, tmp_path):
    out_dir = tmp_path / "comb"
    fastq.comb_fastqs(in_dir=str(fastq_in_dir), out_dir=str(out_dir), out_file="combined.fastq.gz")
    with gzip.open(out_dir / "combined.fastq.gz", "rt") as f:
        combined = f.read()
    assert combined == FASTQ_CONTENT + FASTQ_CONTENT


def test_comb_fastqs_invalid_suffix_raises(fastq_in_dir, tmp_path):
    with pytest.raises(ValueError, match="fastq"):
        fastq.comb_fastqs(in_dir=str(fastq_in_dir), out_dir=str(tmp_path / "comb"), out_file="combined.txt")


def test_comb_fastqs_recursive_one_output_per_subdir(tmp_path):
    in_dir = tmp_path / "in"
    (in_dir / "groupA").mkdir(parents=True)
    (in_dir / "groupB").mkdir(parents=True)
    (in_dir / "groupA" / "x.fastq").write_text(FASTQ_CONTENT)
    (in_dir / "groupB" / "y.fastq").write_text(FASTQ_CONTENT)

    out_dir = tmp_path / "comb_recursive"
    fastq.comb_fastqs(in_dir=str(in_dir), out_dir=str(out_dir), recursive=True, recursive_zip=False)

    assert sorted(os.listdir(out_dir)) == ["groupA.fastq", "groupB.fastq"]
    assert (out_dir / "groupA.fastq").read_text() == FASTQ_CONTENT


# =========================================================================== #
# split_fastqs_by_expected_errors()
# CONFIRMED BUG: the nested helper is called with a `write_gzip` keyword argument
# that its signature does not accept, so the public function raises a TypeError
# on every call regardless of input -- the function is completely unusable.
# fastq.py ~line 560-566 (call site) vs ~line 440-445 (def site).
# =========================================================================== #
def test_split_fastqs_by_expected_errors_is_actually_callable(tmp_path):
    (tmp_path / "sample1.fastq").write_text("@r1\nACGT\n+\nIIII\n")
    out = fastq.split_fastqs_by_expected_errors(indir=str(tmp_path))
    assert out["total_reads"].iloc[0] == 1


def test_split_fastqs_by_expected_errors_splits_reads_by_threshold(tmp_path):
    # r1: all Q40 -> expected_errors ~= 4e-4 (< 1.0 threshold) -> ee0 bucket
    # r2: all Q0  -> expected_errors = 4.0   (>= 1.0 threshold) -> eeplus bucket
    (tmp_path / "sample1.fastq").write_text(
        "@r1\nACGT\n+\nIIII\n"
        "@r2\nACGT\n+\n!!!!\n"
    )
    out = fastq.split_fastqs_by_expected_errors(indir=str(tmp_path), threshold=1.0)
    row = out.iloc[0]
    assert row["total_reads"] == 2
    assert row["reads_expected_errors_lt_threshold"] == 1
    assert row["reads_expected_errors_ge_threshold"] == 1

    sample_dir = tmp_path / "split_fastqs_by_expected_errors" / "sample1"
    with gzip.open(sample_dir / "sample1_EE_lt_1.0.fastq.gz", "rt") as f:
        assert f.read().splitlines()[0] == "@r1"
    with gzip.open(sample_dir / "sample1_EE_ge_1.0.fastq.gz", "rt") as f:
        assert f.read().splitlines()[0] == "@r2"


def test_split_fastqs_by_expected_errors_write_gzip_false_writes_plain_fastq(tmp_path):
    (tmp_path / "sample1.fastq").write_text("@r1\nACGT\n+\nIIII\n")
    fastq.split_fastqs_by_expected_errors(indir=str(tmp_path), write_gzip=False)
    sample_dir = tmp_path / "split_fastqs_by_expected_errors" / "sample1"
    assert (sample_dir / "sample1_EE_lt_1.0.fastq").is_file()
    assert not (sample_dir / "sample1_EE_lt_1.0.fastq.gz").exists()


def test_split_fastqs_by_expected_errors_no_files_raises_filenotfound(tmp_path):
    empty_dir = tmp_path / "empty"
    empty_dir.mkdir()
    with pytest.raises(FileNotFoundError):
        fastq.split_fastqs_by_expected_errors(indir=str(empty_dir))


# =========================================================================== #
# trim_filter()
# The forward (5') quality-trim scan (fastq.py ~L2378-2380) records i+1 -- the
# position just past the last low-quality base -- as the slice start, mirroring
# the reversed (3') scan (~L2381-2383), which records the first low-quality
# index found scanning backward as the (exclusive) slice end. A single
# low-quality base at either end of a read is trimmed symmetrically.
# =========================================================================== #
def _mk_record(seq, quals, rid="r1"):
    r = SeqRecord(Seq(seq), id=rid)
    r.letter_annotations["phred_quality"] = quals
    return r


def test_trim_filter_passthrough_when_all_thresholds_met():
    r = _mk_record("ACGTACGTAC", [40] * 10)
    rid, seq, seqN, quals, alls, avgs, trims, masks = fastq.trim_filter(
        r, qall=10, qavg=30, qtrim=0, qmask=0, alls=0, avgs=0, trims=0, masks=0,
    )
    assert rid == "r1"
    assert str(seq) == "ACGTACGTAC"
    assert alls == 0 and avgs == 0 and trims == 0 and masks == 0


def test_trim_filter_drops_read_below_qall_threshold():
    r = _mk_record("ACGTACGTAC", [40, 40, 40, 5, 40, 40, 40, 40, 40, 40])
    rid, seq, seqN, quals, alls, avgs, trims, masks = fastq.trim_filter(
        r, qall=10, qavg=30, qtrim=0, qmask=0, alls=0, avgs=0, trims=0, masks=0,
    )
    assert rid is None and seq is None
    assert alls == 1 and avgs == 0


def test_trim_filter_drops_read_below_qavg_threshold():
    r = _mk_record("ACGTACGTAC", [15] * 10)
    rid, seq, seqN, quals, alls, avgs, trims, masks = fastq.trim_filter(
        r, qall=10, qavg=30, qtrim=0, qmask=0, alls=0, avgs=0, trims=0, masks=0,
    )
    assert rid is None
    assert alls == 0 and avgs == 1


def test_trim_filter_masks_low_quality_bases_to_N():
    quals = [40, 40, 5, 40, 40, 40, 40, 40, 40, 40]
    r = _mk_record("ACGTACGTAC", quals)
    rid, seq, seqN, out_quals, alls, avgs, trims, masks = fastq.trim_filter(
        r, qall=1, qavg=1, qtrim=0, qmask=10, alls=0, avgs=0, trims=0, masks=1,
    )
    assert str(seq) == "ACGTACGTAC"   # unmasked sequence untouched
    assert str(seqN) == "ACNTACGTAC"  # low-quality base #2 masked to N
    assert masks == 2                 # counter incremented from initial value of 1


def test_trim_filter_3prime_trims_single_trailing_low_quality_base():
    quals = [40, 40, 40, 40, 40, 40, 40, 40, 40, 5]
    r = _mk_record("ACGTACGTAC", quals)
    rid, seq, seqN, out_quals, alls, avgs, trims, masks = fastq.trim_filter(
        r, qall=1, qavg=1, qtrim=10, qmask=0, alls=0, avgs=0, trims=0, masks=0,
    )
    assert str(seq) == "ACGTACGTA"  # trailing low-quality base correctly dropped
    assert trims == 1


def test_trim_filter_5prime_trims_single_leading_low_quality_base():
    quals = [5, 40, 40, 40, 40, 40, 40, 40, 40, 40]
    r = _mk_record("ACGTACGTAC", quals)
    rid, seq, seqN, out_quals, alls, avgs, trims, masks = fastq.trim_filter(
        r, qall=1, qavg=1, qtrim=10, qmask=0, alls=0, avgs=0, trims=0, masks=0,
    )
    assert str(seq) == "CGTACGTAC"  # expect the single bad leading base dropped
    assert trims == 1


def test_trim_filter_qtrim_disabled_skips_trimming():
    quals = [5, 5, 40, 40, 40, 40, 40, 40, 5, 5]
    r = _mk_record("ACGTACGTAC", quals)
    rid, seq, seqN, out_quals, alls, avgs, trims, masks = fastq.trim_filter(
        r, qall=1, qavg=1, qtrim=0, qmask=0, alls=0, avgs=0, trims=0, masks=0,
    )
    assert str(seq) == "ACGTACGTAC"
    assert trims == 0


# =========================================================================== #
# find_AA_edits(), trim(), format_alignment(), find_indel()
# =========================================================================== #
def test_find_AA_edits_single_substitution():
    assert fastq.find_AA_edits(wt="MAK", res=1, seq="MAN") == "K3N"


def test_find_AA_edits_no_difference_returns_unknown():
    # find_AA_edits() is only ever called by genotype() after it has already
    # confirmed wt != seq (the wt==seq case is special-cased to 'WT' upstream), so
    # 'Unknown' here just documents find_AA_edits()'s own standalone behavior.
    assert fastq.find_AA_edits(wt="MAK", res=1, seq="MAK") == "Unknown"


def test_find_AA_edits_multiple_substitutions_joined_with_comma():
    assert fastq.find_AA_edits(wt="MAKG", res=1, seq="MANA") == "K3N, G4A"


def test_find_AA_edits_respects_res_offset():
    assert fastq.find_AA_edits(wt="MAK", res=101, seq="MAN") == "K103N"


def test_trim_truncates_to_multiple_of_three():
    assert fastq.trim("ATGGCTA") == "ATGGCT"
    assert str(fastq.trim(Seq("ATGGCTA"))) == "ATGGCT"


def test_trim_already_multiple_of_three_unchanged():
    assert fastq.trim("ATGGCTAAA") == "ATGGCTAAA"


def test_format_alignment_middle_string():
    mid = fastq.format_alignment("ATG-CT", "ATGACT")
    assert mid == "|||-||"


def test_format_alignment_return_full_text():
    full = fastq.format_alignment("AT", "AC", return_alignment=True)
    assert full == "AT\n|.\nAC"


def test_find_indel_in_frame_deletion():
    wt = "ATGGCTAAA"      # M A K
    mut = "ATGAAA"        # M K   (codon 'GCT'/A deleted)
    edit, category = fastq.find_indel(wt=wt, mut=mut, res=1)
    assert category == "In-frame Deletion"
    assert edit == "AK2K"


def test_find_indel_in_frame_insertion():
    wt = "ATGGCTAAA"          # M A K
    mut = "ATGGCTGCTAAA"      # M A A K (extra 'GCT'/A codon inserted)
    edit, category = fastq.find_indel(wt=wt, mut=mut, res=1)
    assert category == "In-frame Insertion"
    assert edit == "M1MA"


def test_find_indel_frameshift():
    wt = "ATGGCTAAA"     # 9 nt
    mut = "ATGGCAAA"     # 8 nt -> not a multiple of 3
    edit, category = fastq.find_indel(wt=wt, mut=mut, res=1)
    assert category == "Frameshift Indel"
    assert edit == "K3fs"


# =========================================================================== #
# genotype(), outcomes(), abundances()
# =========================================================================== #
@pytest.fixture
def wt_seq():
    return "ATGGCTAAA"  # M A K


def test_genotype_calls_wt_reads_as_WT(wt_seq):
    fastqs = {"s1": pd.DataFrame({"nuc": [wt_seq], "prot": [str(Seq(wt_seq).translate())]})}
    out = fastq.genotype(fastqs=fastqs, res=1, wt=wt_seq, save=False)
    assert out["s1"].iloc[0]["Edit"] == "WT"
    assert out["s1"].iloc[0]["Category"] == "WT"


def test_genotype_calls_substitution(wt_seq):
    mut_nuc = "ATGGCTAAC"  # M A N  (K -> N)
    fastqs = {"s1": pd.DataFrame({"nuc": [mut_nuc], "prot": [str(Seq(mut_nuc).translate())]})}
    out = fastq.genotype(fastqs=fastqs, res=1, wt=wt_seq, save=False)
    assert out["s1"].iloc[0]["Edit"] == "K3N"
    assert out["s1"].iloc[0]["Category"] == "Substitution"


def test_genotype_calls_indel_when_length_differs(wt_seq):
    del_nuc = "ATGAAA"  # deletion relative to wt_seq
    fastqs = {"s1": pd.DataFrame({"nuc": [del_nuc], "prot": [str(Seq(fastq.trim(del_nuc)).translate())]})}
    out = fastq.genotype(fastqs=fastqs, res=1, wt=wt_seq, save=False)
    assert out["s1"].iloc[0]["Edit"] == "AK2K"
    assert out["s1"].iloc[0]["Category"] == "In-frame Deletion"


def test_genotype_empty_translation_calls_unknown_flanks(wt_seq):
    fastqs = {"s1": pd.DataFrame({"nuc": [wt_seq], "prot": [""]})}
    out = fastq.genotype(fastqs=fastqs, res=1, wt=wt_seq, save=False)
    assert out["s1"].iloc[0]["Edit"] == "Unknown"
    assert out["s1"].iloc[0]["Category"] == "Flanks"


def test_genotype_rejects_out_of_frame_wt():
    with pytest.raises(ValueError, match="not in-frame"):
        fastq.genotype(fastqs={}, res=1, wt="ATGGC", save=False)


def test_outcomes_counts_and_fractions(wt_seq):
    fastqs = {"s1": pd.DataFrame({
        "nuc": [wt_seq, "ATGGCTAAC", "ATGAAA", ""],
        "prot": [
            str(Seq(wt_seq).translate()),
            str(Seq("ATGGCTAAC").translate()),
            str(Seq(fastq.trim("ATGAAA")).translate()),
            "",
        ],
    })}
    genotyped = fastq.genotype(fastqs=fastqs, res=1, wt=wt_seq, save=False)
    out = fastq.outcomes(fastqs=genotyped)

    assert set(out["fastq_file"]) == {"s1"}
    assert set(out["count"]) == {1}
    assert np.allclose(out["fraction"], 0.25)
    assert set(out["Edit"]) == {"WT", "K3N", "AK2K", "Unknown"}


def test_abundances_isolates_single_desired_edit():
    df = pd.DataFrame({
        "fastq_file": ["s1"] * 3,
        "Edit": ["K3N", "WT", "A5G"],
        "count": [5, 10, 3],
        "fraction": [0.28, 0.56, 0.16],
    })
    out = fastq.abundances(df=df, desired_edits=["K3N"])
    assert list(out["Edit"]) == ["K3N"]


def test_abundances_combinations_includes_paired_edit_strings():
    df = pd.DataFrame({
        "fastq_file": ["s1"] * 4,
        "Edit": ["K3N", "A5G", "K3N, A5G", "WT"],
        "count": [5, 3, 2, 10],
        "fraction": [0.25, 0.15, 0.10, 0.50],
    })
    out = fastq.abundances(df=df, desired_edits=["K3N", "A5G"], combinations=2)
    assert set(out["Edit"]) == {"K3N", "A5G", "K3N, A5G"}


def test_abundances_reads_from_file_path(tmp_path):
    csv_pt = tmp_path / "outcomes.csv"
    pd.DataFrame({
        "fastq_file": ["s1"],
        "Edit": ["K3N"],
        "count": [5],
        "fraction": [1.0],
    }).to_csv(csv_pt, index=False)
    out = fastq.abundances(df=str(csv_pt), desired_edits=["K3N"])
    assert list(out["Edit"]) == ["K3N"]


# =========================================================================== #
# edit_change()
# =========================================================================== #
def test_edit_change_substitution_classification():
    df = pd.DataFrame({"Edit": ["K3N"]})
    out = fastq.edit_change(df, edit="Edit")
    row = out.iloc[0]
    assert row["Before"] == "K" and row["After"] == "N" and row["AA Number"] == 3
    assert row["Initial AA"] == "K" and row["Resulting AA"] == "N"
    assert row["Change"] == "Polar"  # K (basic) -> N (polar)


def test_edit_change_excludes_WT_and_multi_edit_rows():
    df = pd.DataFrame({"Edit": ["K3N", "WT", "K3N, A5G", "Not WT"]})
    out = fastq.edit_change(df, edit="Edit")
    assert list(out["Edit"]) == ["K3N"]


def test_edit_change_insertion_classification():
    df = pd.DataFrame({"Edit": ["A5AG"]})
    out = fastq.edit_change(df, edit="Edit")
    row = out.iloc[0]
    assert row["Initial AA"] == "A" and row["Resulting AA"] == "insG"


def test_edit_change_deletion_classification():
    df = pd.DataFrame({"Edit": ["GA5A"]})
    out = fastq.edit_change(df, edit="Edit")
    row = out.iloc[0]
    assert row["Initial AA"] == "G" and row["Resulting AA"] == "del"


def test_edit_change_includes_aa_properties_when_requested():
    df = pd.DataFrame({"Edit": ["K3N"]})
    out = fastq.edit_change(df, edit="Edit", aa_properties=["hydrophobicity"])
    assert "Before_AA_Properties" in out.columns
    assert "hydrophobicity" in out.iloc[0]["Before_AA_Properties"]


# =========================================================================== #
# _signature_phred_summary()
# =========================================================================== #
def _global_aligner():
    aligner = PairwiseAligner()
    aligner.mode = "global"
    aligner.match_score = 2
    aligner.mismatch_score = -1
    aligner.open_gap_score = -10
    aligner.extend_gap_score = -0.1
    return aligner


def test_signature_phred_summary_snv_reports_read_base_quality():
    ref = "ACGTACGTAC"
    query = "ACGTGCGTAC"  # SNV at 0-based pos 4: A->G
    aligner = _global_aligner()
    aln = aligner.align(ref, query)[0]
    sig = signature_from_alignment(aln, ref_seq=ref, query_seq=query)

    read_quals = [30] * len(query)
    read_quals[4] = 10
    result = fastq._signature_phred_summary(sig, aln, read_quals)
    snv_min, snv_med, snv_mean, snv_max, read_min, read_med, read_mean, read_max = result
    assert snv_min == 10 and snv_max == 10
    assert read_min == 10 and read_max == 30


def test_signature_phred_summary_no_variant_returns_nan():
    ref = "ACGTACGTAC"
    aligner = _global_aligner()
    aln = aligner.align(ref, ref)[0]
    sig = signature_from_alignment(aln, ref_seq=ref, query_seq=ref)
    result = fastq._signature_phred_summary(sig, aln, [30] * len(ref))
    assert all(np.isnan(v) for v in result[:4])


def test_signature_phred_summary_insertion_uses_inserted_base_quality():
    ref = "ACGTACGTAC"
    query = "ACGTAACGTAC"  # inserted extra 'A' after 0-based ref pos 4
    aligner = _global_aligner()
    aln = aligner.align(ref, query)[0]
    sig = signature_from_alignment(aln, ref_seq=ref, query_seq=query)
    assert len(sig.indels) == 1 and sig.indels[0].ins != ""

    read_quals = [30] * len(query)
    result = fastq._signature_phred_summary(sig, aln, read_quals)
    assert not np.isnan(result[0])  # an inserted-base quality was found


def test_signature_phred_summary_deletion_uses_flanking_base_quality():
    ref = "ACGTACGTAC"
    query = "ACGTCGTAC"  # 1-nt deletion relative to ref
    aligner = _global_aligner()
    aln = aligner.align(ref, query)[0]
    sig = signature_from_alignment(aln, ref_seq=ref, query_seq=query)
    assert len(sig.indels) == 1 and sig.indels[0].dellen > 0

    read_quals = [30] * len(query)
    result = fastq._signature_phred_summary(sig, aln, read_quals)
    assert not np.isnan(result[0])


# =========================================================================== #
# extract_umis() / trim_motifs()
# Subprocess command construction only -- umi_tools/cutadapt are not installed in
# this env (they live in the separate 'umi_tools' conda env per the README), so
# subprocess.run is mocked and we assert on the constructed command string(s).
# =========================================================================== #
class _FakeCompletedProcess:
    def __init__(self):
        self.stdout = ""
        self.stderr = ""


def test_extract_umis_constructs_umi_tools_command(tmp_path, mocker):
    fastq_dir = tmp_path / "fastqs"
    fastq_dir.mkdir()
    (fastq_dir / "sample1.fastq").write_text(FASTQ_CONTENT)
    out_dir = tmp_path / "out"

    mock_run = mocker.patch("edms.bio.fastq.subprocess.run", return_value=_FakeCompletedProcess())
    fastq.extract_umis(fastq_dir=str(fastq_dir), out_dir=str(out_dir), bc_pattern="NNNN", env="myenv")

    assert mock_run.call_count == 1
    command = mock_run.call_args[0][0]
    assert "conda run -n myenv umi_tools extract" in command
    assert "--bc-pattern=NNNN" in command
    assert f"--stdin={os.path.join(str(fastq_dir), 'sample1.fastq')}" in command
    assert mock_run.call_args.kwargs["shell"] is True


def test_extract_umis_log_path_matches_created_directory(tmp_path, mocker):
    """extract_umis() creates the logging subdirectory
    os.path.join(out_dir, '.extract_umis') (plural) via mkdir(), and the
    constructed umi_tools --log path must reference that same directory."""
    fastq_dir = tmp_path / "fastqs"
    fastq_dir.mkdir()
    (fastq_dir / "sample1.fastq").write_text(FASTQ_CONTENT)
    out_dir = tmp_path / "out"

    created_dirs = []
    real_mkdir = fastq.mkdir

    def _tracking_mkdir(dir, *a, **kw):
        created_dirs.append(os.path.abspath(str(dir)))
        return real_mkdir(dir, *a, **kw)

    mocker.patch("edms.bio.fastq.mkdir", side_effect=_tracking_mkdir)
    mock_run = mocker.patch("edms.bio.fastq.subprocess.run", return_value=_FakeCompletedProcess())
    fastq.extract_umis(fastq_dir=str(fastq_dir), out_dir=str(out_dir), bc_pattern="NNNN")

    command = mock_run.call_args[0][0]
    # Extract the directory referenced by --log=<dir>/<file>.log
    log_arg = [tok for tok in command.split() if tok.startswith("--log=")][0]
    log_dir = os.path.dirname(log_arg[len("--log="):])
    assert log_dir in created_dirs


def test_trim_motifs_derives_motif5_motif3_from_target_sequence(tmp_path, mocker):
    fastq_dir = tmp_path / "fastqs"
    fastq_dir.mkdir()
    (fastq_dir / "sample1.fastq").write_text(FASTQ_CONTENT)

    in_file = pd.DataFrame({
        "target_name": ["t1"],
        "target_sequence": ["AAAAAAAAAAGGGGGGGGGG(ATCG)CCCCCCCCCCTTTTTTTTTT"],
    })

    mock_run = mocker.patch("edms.bio.fastq.subprocess.run", return_value=_FakeCompletedProcess())
    fastq.trim_motifs(fastq_dir=str(fastq_dir), out_dir=str(tmp_path / "out"), in_file=in_file, motif_length=10)

    assert mock_run.call_count == 2
    first_command = mock_run.call_args_list[0][0][0]
    second_command = mock_run.call_args_list[1][0][0]
    assert "-g GGGGGGGGGG" in first_command
    assert "-a CCCCCCCCCC" in second_command


def test_trim_motifs_requires_motif_source(tmp_path):
    fastq_dir = tmp_path / "fastqs"
    fastq_dir.mkdir()
    with pytest.raises(ValueError, match="in_file"):
        fastq.trim_motifs(fastq_dir=str(fastq_dir), out_dir=str(tmp_path / "out"))
