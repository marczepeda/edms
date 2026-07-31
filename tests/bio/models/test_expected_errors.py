"""
Tests for edms.bio.models.expected_errors — expected errors (EE) vs average
Phred Q for multiple read lengths.

Pure-math helpers are cross-checked against hand-derived formulas;
build_parser()/main() get light smoke tests.
"""
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pytest

from edms.bio.models import expected_errors as ee


@pytest.fixture(autouse=True)
def _close_figures():
    yield
    plt.close("all")


def test_p_error_from_Q_matches_phred_formula():
    assert ee.p_error_from_Q(10) == pytest.approx(0.1)
    assert ee.p_error_from_Q(20) == pytest.approx(0.01)
    assert ee.p_error_from_Q(30) == pytest.approx(0.001)
    assert ee.p_error_from_Q(0) == pytest.approx(1.0)


def test_expected_errors_from_avgQ():
    assert ee.expected_errors_from_avgQ(30, 100) == pytest.approx(100 * 10 ** (-30 / 10))
    assert ee.expected_errors_from_avgQ(30, 100) == pytest.approx(0.1)


def test_expected_errors_from_percycle_Q():
    Qs = [30, 30, 20]
    expected = sum(10 ** (-q / 10) for q in Qs)
    assert ee.expected_errors_from_percycle_Q(Qs) == pytest.approx(expected)
    assert ee.expected_errors_from_percycle_Q(Qs) == pytest.approx(0.012)


def test_expected_errors_from_percycle_Q_accepts_array():
    result = ee.expected_errors_from_percycle_Q(np.array([40.0, 40.0]))
    assert result == pytest.approx(2 * 10 ** (-4))


def test_build_parser_defaults():
    parser = ee.build_parser()
    args = parser.parse_args([])
    assert args.qmin == 10
    assert args.qmax == 45
    assert args.cycles == [100, 200, 300, 400, 500, 600]


def test_build_parser_overrides():
    parser = ee.build_parser()
    args = parser.parse_args(["--qmin", "5", "--cycles", "150", "250"])
    assert args.qmin == 5
    assert args.cycles == [150, 250]


def test_main_smoke(tmp_path, monkeypatch):
    out_csv = tmp_path / "ee.csv"
    monkeypatch.setattr(sys, "argv", [
        "prog", "--no-show", "--csv-out", str(out_csv), "--cycles", "100", "200",
    ])
    ee.main()
    assert out_csv.exists()
