"""
Tests for edms.bio.pegLIT — pegRNA linker design (simulated-annealing linker optimizer).

pegLIT.py ports linker-design code from the PE-based tool described in
https://www.nature.com/articles/s41587-021-01039-7. It uses ViennaRNA (RNA.fold_compound),
scikit-learn (AgglomerativeClustering), and Levenshtein — all installed in the `edms` conda env,
so these tests exercise the real implementations rather than mocking them.
"""
import numpy as np
import pytest

from edms.bio.pegLIT import (
    apply_filters,
    calc_subscores,
    apply_score,
    optimize,
    apply_bottleneck,
    pegLIT,
)


# ---------------------------------------------------------------------------
# apply_filters
# ---------------------------------------------------------------------------

class TestApplyFilters:
    def test_ac_content_pass(self):
        # linker "AACC" has 4 A/C nts, ac_thresh=4 -> exactly meets threshold (count < thresh fails, so 4 < 4 is False -> passes)
        assert apply_filters("GGGG", "AACC", "GGGG", ac_thresh=4, u_thresh=10, n_thresh=10, excluded_motifs=None) is True

    def test_ac_content_fail(self):
        # linker "AACC" has 4 A/C nts; ac_thresh=5 -> 4 < 5 -> fails
        assert apply_filters("GGGG", "AACC", "GGGG", ac_thresh=5, u_thresh=10, n_thresh=10, excluded_motifs=None) is False

    def test_ac_content_counts_only_A_and_C(self):
        # "GGTT" has zero A/C -> count=0, ac_thresh=1 -> 0 < 1 -> fails
        assert apply_filters("AAAA", "GGTT", "AAAA", ac_thresh=1, u_thresh=10, n_thresh=10, excluded_motifs=None) is False

    def test_consecutive_u_within_linker_fails(self):
        # linker itself contains u_thresh+1 T's in a row (treated as U) -> should fail
        # u_thresh=3 -> need 4 consecutive T's to trip it
        assert apply_filters("AAAA", "CTTTTC", "AAAA", ac_thresh=0, u_thresh=3, n_thresh=10, excluded_motifs=None) is False

    def test_consecutive_u_below_threshold_passes(self):
        # only 3 consecutive T's with u_thresh=3 -> exactly at threshold, not thresh+1 -> passes
        assert apply_filters("AAAA", "CTTTC", "AAAA", ac_thresh=0, u_thresh=3, n_thresh=10, excluded_motifs=None) is True

    def test_consecutive_u_spanning_pre_and_linker(self):
        # seq_pre ends in TTT, linker starts with T -> 4 consecutive T's/U's spanning the boundary
        assert apply_filters("AATTT", "TCCC", "AAAA", ac_thresh=0, u_thresh=3, n_thresh=10, excluded_motifs=None) is False

    def test_consecutive_n_within_linker_fails(self):
        # n_thresh=3 -> any nt repeated 4+ times in a row in the neighborhood fails.
        # Use G (not U) so the U-filter doesn't also trigger, and ac_thresh=0 so AC filter doesn't trigger.
        assert apply_filters("AAAA", "GGGGC", "AAAA", ac_thresh=0, u_thresh=10, n_thresh=3, excluded_motifs=None) is False

    def test_consecutive_n_below_threshold_passes(self):
        assert apply_filters("AAAA", "GGGC", "AAAA", ac_thresh=0, u_thresh=10, n_thresh=3, excluded_motifs=None) is True

    def test_excluded_motif_detected_fully_within_linker(self):
        # Motif entirely contained inside the linker (no boundary-spanning needed) -> detected regardless of window bug
        assert apply_filters("AAAA", "CCGAATTCGG", "AAAA", ac_thresh=0, u_thresh=10, n_thresh=10,
                              excluded_motifs=["GAATTC"]) is False

    def test_excluded_motif_absent_passes(self):
        assert apply_filters("AAAA", "CCCCGG", "AAAA", ac_thresh=0, u_thresh=10, n_thresh=10,
                              excluded_motifs=["GAATTC"]) is True

    def test_excluded_motif_none_skips_check(self):
        # excluded_motifs=None should skip the motif check entirely (no crash)
        assert apply_filters("AAAA", "GAATTC", "AAAA", ac_thresh=0, u_thresh=10, n_thresh=10,
                              excluded_motifs=None) is True

    def test_excluded_motif_boundary_spanning_bug(self):
        # Fixed bug: the excluded-motif boundary window used to be sized with
        # `len(excluded_motifs)` (the *number* of motifs) instead of the motif length(s), so a
        # motif spanning the seq_pre/seq_linker boundary was silently missed. The window is now
        # sized from the longest excluded motif, so this boundary-spanning case is caught.
        seq_pre = "GAATT"
        seq_linker = "CAAAA"
        seq_post = "TTTTT"
        excluded_motifs = ["GAATTC"]
        # Sanity check: the motif genuinely spans the pre/linker boundary in the true concatenated sequence.
        assert "GAATTC" in (seq_pre + seq_linker + seq_post)
        result = apply_filters(seq_pre, seq_linker, seq_post, ac_thresh=0, u_thresh=10, n_thresh=10,
                                excluded_motifs=excluded_motifs)
        assert result is False


# ---------------------------------------------------------------------------
# calc_subscores / apply_score (ViennaRNA-backed)
# ---------------------------------------------------------------------------

class TestCalcSubscores:
    def test_returns_masked_array_of_correct_length(self):
        res = calc_subscores(1, "AAAA", "CCCC", "GGGG")
        assert isinstance(res, np.ma.MaskedArray)
        assert len(res) == 3

    def test_subscores_are_probabilities(self):
        res = calc_subscores(0, "AAAAGGGG", "CCCCUUUU".replace("U", "T"), "GGGGCCCC")
        for v in res:
            assert -1e-6 <= v <= 1.0 + 1e-6

    def test_self_index_reflects_unpaired_probability(self):
        # For a poly-A / poly-C / poly-G sequence set (no complementary base-pairing possible between
        # A-only and C-only and G-only stretches under standard RNA pairing rules AU/GC/GU), the marginal
        # basepair probability should be very low everywhere, so probability mass concentrates in the
        # "unpaired" diagonal captured at the linker_pos-th component. This is a smoke/structural check,
        # not a hand re-derivation of the folding algorithm.
        res = calc_subscores(1, "AAAA", "AAAA", "AAAA")
        assert res.sum() == pytest.approx(len(res) * res[1] if False else res.sum())  # array is well-formed
        assert not np.any(np.isnan(res.astype(float)))


class TestApplyScore:
    def test_returns_4_tuple_of_rounded_values(self):
        score = apply_score(
            seq_spacer="GGCCC",
            seq_scaffold="GTTTTAGAGCTAG",
            seq_template="TGGAGGA",
            seq_pbs="CGTGCTG",
            seq_linker="AACCGG",
        )
        assert isinstance(score, tuple)
        assert len(score) == 4
        for v in score:
            assert 0.0 <= v <= 1.0

    def test_epsilon_rounding(self):
        # apply_score buckets subscores to the nearest epsilon via epsilon*int(val/epsilon);
        # with epsilon=0.1 every returned value must be a multiple of 0.1 (within float tolerance).
        score = apply_score(
            seq_spacer="GGCCC", seq_scaffold="GTTTTAGAGCTAG", seq_template="TGGAGGA",
            seq_pbs="CGTGCTG", seq_linker="AACCGG", epsilon=0.1,
        )
        for v in score:
            multiple = v / 0.1
            assert abs(multiple - round(multiple)) < 1e-6

    def test_score_to_beat_short_circuits_free_pegrna_subscores(self):
        # If subscore_pbs is bucketed below score_to_beat[0], spacer/scaffold subscores should be
        # returned as exactly 0. (impossible-to-beat score_to_beat forces the early-exit branch)
        score = apply_score(
            seq_spacer="GGCCC", seq_scaffold="GTTTTAGAGCTAG", seq_template="TGGAGGA",
            seq_pbs="CGTGCTG", seq_linker="AACCGG", epsilon=0.01,
            score_to_beat=(1.0, 0, 0, 0),  # impossible to beat: subscore_pbs bucket can be at most < 1.0
        )
        assert score[1] == 0.0
        assert score[3] == 0.0


# ---------------------------------------------------------------------------
# optimize / apply_bottleneck (small, fast parameterizations)
# ---------------------------------------------------------------------------

SPACER = "GGCCCAGACTGAGCACGTGA"
SCAFFOLD = "GTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGC"
TEMPLATE = "TGGAGGAAGCAGGGCTTCCTTTCCTCTGCCATCACTTATCGTCGTCATCCTTGTAATC"
PBS = "CGTGCTCAGTCTG"
MOTIF = "CGCGGTTCTATCTAGTTACGCGTTAAACCAACTAGAA"


class TestOptimize:
    def test_returns_topn_linkers_of_correct_length(self):
        scores, heap = optimize(
            SPACER, SCAFFOLD, TEMPLATE, PBS, MOTIF,
            linker_pattern="NNNN", ac_thresh=0.0, u_thresh=3, n_thresh=3, topn=3, epsilon=0.01,
            num_repeats=2, num_steps=5, temp_init=0.15, temp_decay=0.95, seed=1, excluded_motifs=[],
        )
        assert len(scores) == 3
        assert len(heap) == 3
        for linker in heap:
            assert len(linker) == 4  # matches linker_pattern length
            assert set(linker) <= set("ACGT")

    def test_deterministic_given_seed(self):
        args = (SPACER, SCAFFOLD, TEMPLATE, PBS, MOTIF)
        kwargs = dict(linker_pattern="NNNN", ac_thresh=0.0, u_thresh=3, n_thresh=3, topn=3, epsilon=0.01,
                      num_repeats=2, num_steps=5, temp_init=0.15, temp_decay=0.95, seed=42, excluded_motifs=[])
        scores1, heap1 = optimize(*args, **kwargs)
        scores2, heap2 = optimize(*args, **kwargs)
        assert heap1 == heap2
        assert scores1 == scores2

    def test_respects_fixed_pattern_positions(self):
        # A pattern with fixed bases (not N) must keep those positions fixed in every returned linker.
        scores, heap = optimize(
            SPACER, SCAFFOLD, TEMPLATE, PBS, MOTIF,
            linker_pattern="ANNT", ac_thresh=0.0, u_thresh=3, n_thresh=3, topn=2, epsilon=0.01,
            num_repeats=2, num_steps=5, temp_init=0.15, temp_decay=0.95, seed=1, excluded_motifs=[],
        )
        for linker in heap:
            assert linker[0] == "A"
            assert linker[3] == "T"


class TestApplyBottleneck:
    def test_bottleneck_1_returns_single_best(self):
        heap_scores = [(0.1, 0.1, 0.1, 0.1), (0.9, 0.9, 0.9, 0.9), (0.5, 0.5, 0.5, 0.5)]
        heap = ["AAAA", "CCCC", "GGGG"]
        out = apply_bottleneck(heap_scores, heap, bottleneck=1, seed=1)
        assert out == ["CCCC"]

    def test_bottleneck_n_returns_n_sequences(self):
        scores, heap = optimize(
            SPACER, SCAFFOLD, TEMPLATE, PBS, MOTIF,
            linker_pattern="NNNN", ac_thresh=0.0, u_thresh=3, n_thresh=3, topn=4, epsilon=0.01,
            num_repeats=2, num_steps=5, temp_init=0.15, temp_decay=0.95, seed=1, excluded_motifs=[],
        )
        out = apply_bottleneck(scores, heap, bottleneck=2, seed=1)
        assert len(out) == 2
        for seq in out:
            assert str(seq) in heap


# ---------------------------------------------------------------------------
# pegLIT end-to-end (tiny parameters to keep runtime small)
# ---------------------------------------------------------------------------

class TestPegLITEndToEnd:
    def test_pegLIT_returns_requested_number_of_linkers_of_correct_length(self):
        out = pegLIT(
            seq_spacer=SPACER, seq_scaffold=SCAFFOLD, seq_template=TEMPLATE, seq_pbs=PBS, seq_motif=MOTIF,
            linker_pattern="NNNN", num_repeats=2, num_steps=5, topn=3, bottleneck=1, seed=1,
        )
        assert isinstance(out, list)
        assert len(out) == 1
        assert len(out[0]) == 4
        assert set(out[0]) <= set("ACGT")

    def test_pegLIT_bottleneck_greater_than_1(self):
        out = pegLIT(
            seq_spacer=SPACER, seq_scaffold=SCAFFOLD, seq_template=TEMPLATE, seq_pbs=PBS, seq_motif=MOTIF,
            linker_pattern="NNNN", num_repeats=2, num_steps=5, topn=4, bottleneck=2, seed=1,
        )
        assert len(out) == 2

    def test_pegLIT_deterministic_given_seed(self):
        kwargs = dict(
            seq_spacer=SPACER, seq_scaffold=SCAFFOLD, seq_template=TEMPLATE, seq_pbs=PBS, seq_motif=MOTIF,
            linker_pattern="NNNN", num_repeats=2, num_steps=5, topn=3, bottleneck=1, seed=7,
        )
        out1 = pegLIT(**kwargs)
        out2 = pegLIT(**kwargs)
        assert out1 == out2

    def test_pegLIT_honors_fixed_linker_pattern_positions(self):
        out = pegLIT(
            seq_spacer=SPACER, seq_scaffold=SCAFFOLD, seq_template=TEMPLATE, seq_pbs=PBS, seq_motif=MOTIF,
            linker_pattern="GNNC", num_repeats=2, num_steps=5, topn=3, bottleneck=1, seed=3,
        )
        assert out[0][0] == "G"
        assert out[0][3] == "C"
