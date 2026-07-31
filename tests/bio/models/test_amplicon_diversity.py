"""
Tests for edms.bio.models.amplicon_diversity — expected per-cycle nucleotide
diversity (Shannon entropy) for staggered forward primers.

Pure-math helpers are cross-checked against hand-derived expectations;
main() gets a light smoke test.
"""
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pytest

from edms.bio.models import amplicon_diversity as ad


@pytest.fixture(autouse=True)
def _close_figures():
    yield
    plt.close("all")


def test_clean_seq_uppercases_and_strips_non_acgt():
    assert ad.clean_seq("acgtACGTnN-xyz123") == "ACGTACGT"


def test_base_index():
    assert [ad.base_index(b) for b in "ACGT"] == [0, 1, 2, 3]


def test_per_cycle_base_freqs_no_stagger_matches_sequence():
    # With S=0 there is only stagger s=0. At cycle c=0, the code's `c <= s`
    # branch (0<=0) is taken -> uniform 0.25 for every base (this is the
    # model's convention for the first cycle under zero stagger). For c>=1,
    # pos=c-0=c reads the base at that position directly -> one-hot.
    seq = "ACGT"
    freqs = ad.per_cycle_base_freqs(seq, S=0, max_cycles=4)
    expected = np.array([
        [0.25, 0.25, 0.25, 0.25],
        [0, 1, 0, 0],
        [0, 0, 1, 0],
        [0, 0, 0, 1],
    ])
    assert np.allclose(freqs, expected)


def test_per_cycle_base_freqs_pools_over_stagger():
    # S=1 pools stagger s=0 and s=1 (2 options), each weighted 1/(S+1)=0.5
    seq = "AC"
    freqs = ad.per_cycle_base_freqs(seq, S=1, max_cycles=2)
    # cycle c=0: s=0 -> c<=s (0<=0) True -> uniform 0.25 each;
    #            s=1 -> c<=s (0<=1) True -> uniform 0.25 each
    # average of two uniform distributions is still uniform
    assert np.allclose(freqs[0], [0.25, 0.25, 0.25, 0.25])
    # cycle c=1: s=0 -> c<=s (1<=0) False -> pos=1-0=1 -> seq[1]='C' -> one-hot C
    #            s=1 -> c<=s (1<=1) True -> uniform 0.25 each
    # average: C gets (1+0.25)/2=0.625, others get 0.25/2=0.125
    expected_c1 = np.array([0.125, 0.625, 0.125, 0.125])
    assert np.allclose(freqs[1], expected_c1)


def test_shannon_entropy_bits_uniform_is_two_bits():
    p = np.array([0.25, 0.25, 0.25, 0.25])
    assert ad.shannon_entropy_bits(p) == pytest.approx(2.0, abs=1e-6)


def test_shannon_entropy_bits_certain_is_zero_bits():
    p = np.array([1.0, 0.0, 0.0, 0.0])
    assert ad.shannon_entropy_bits(p) == pytest.approx(0.0, abs=1e-6)


def test_shannon_entropy_bits_matches_manual_log2_sum():
    p = np.array([0.5, 0.25, 0.125, 0.125])
    expected = -np.sum(p * np.log2(p))
    assert ad.shannon_entropy_bits(p) == pytest.approx(expected, abs=1e-6)


def test_main_smoke(tmp_path, monkeypatch):
    out_csv = tmp_path / "out.csv"
    monkeypatch.setattr(sys, "argv", [
        "prog", "--seq", "ACGTACGTACGTACGTACGT",
        "--S_values", "0", "2",
        "--out_csv", str(out_csv),
    ])
    ad.main()
    assert out_csv.exists()
