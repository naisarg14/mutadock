"""Tests for the ΔΔG protocol helpers in ``mutadock.mutation.predict_ddG``.

PyRosetta is stubbed by this package's conftest, so only the pure-Python
summary/scaling logic is exercised here; the real repack/minimize protocol is
verified end-to-end elsewhere.
"""

import statistics

import pytest

from mutadock.mutation import predict_ddG


def test_scale_constant_is_sane():
    assert isinstance(predict_ddG.REU_TO_KCAL_SCALE, float)
    # The cited ref2015 cartesian_ddg range is roughly 0.29–0.40.
    assert 0.2 < predict_ddG.REU_TO_KCAL_SCALE < 0.5


def test_summarize_single_replicate():
    s = predict_ddG.summarize_ddg([-5.0], scale=0.34)
    assert s["ddG_value"] == pytest.approx(-5.0)
    assert s["ddG_sd"] == 0.0  # SD undefined for n=1 → reported as 0
    assert s["ddG_kcal"] == pytest.approx(-1.7)
    assert s["n_replicates"] == 1


def test_summarize_multiple_replicates():
    values = [-4.0, -6.0, -5.0]
    s = predict_ddG.summarize_ddg(values, scale=0.5)
    assert s["ddG_value"] == pytest.approx(-5.0)
    assert s["ddG_sd"] == pytest.approx(statistics.stdev(values))
    assert s["ddG_kcal"] == pytest.approx(-2.5)
    assert s["n_replicates"] == 3


def test_summarize_uses_default_scale():
    s = predict_ddG.summarize_ddg([-10.0])
    assert s["ddG_kcal"] == pytest.approx(-10.0 * predict_ddG.REU_TO_KCAL_SCALE)


# --------------------------------------------------------------------------- #
# protocol_label + resolve_ddg_params (the fast/min/cartesian mapping)         #
# --------------------------------------------------------------------------- #
def test_protocol_label():
    assert (
        predict_ddG.protocol_label(backbone_minimization=False, cartesian=False)
        == "fast"
    )
    assert (
        predict_ddG.protocol_label(backbone_minimization=True, cartesian=False) == "min"
    )
    assert (
        predict_ddG.protocol_label(backbone_minimization=True, cartesian=True)
        == "cartesian"
    )
    # cartesian implies minimization regardless of the flag
    assert (
        predict_ddG.protocol_label(backbone_minimization=False, cartesian=True)
        == "cartesian"
    )


def _resolve(argv):
    import argparse

    from mutadock.mutation.helpers import add_ddg_protocol_args, resolve_ddg_params

    parser = argparse.ArgumentParser()
    add_ddg_protocol_args(parser)
    return resolve_ddg_params(parser.parse_args(argv))


def test_default_protocol_is_min():
    p = _resolve([])
    assert p["backbone_minimization"] is True
    assert p["cartesian"] is False
    assert p["replicates"] == 1
    assert "reu_to_kcal" not in p  # omitted unless --reu-to-kcal given


def test_protocol_fast_disables_minimization():
    p = _resolve(["--protocol", "fast"])
    assert p["backbone_minimization"] is False
    assert p["cartesian"] is False


def test_protocol_cartesian_enables_both():
    p = _resolve(
        ["--protocol", "cartesian", "--replicates", "3", "--reu-to-kcal", "0.29"]
    )
    assert p["backbone_minimization"] is True
    assert p["cartesian"] is True
    assert p["replicates"] == 3
    assert p["reu_to_kcal"] == pytest.approx(0.29)
