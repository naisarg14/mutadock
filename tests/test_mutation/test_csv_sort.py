"""
Tests for mutadock.mutation.csv_sort
-------------------------------------
Pure-logic tests (pandas only) for sort_csv().

These cover:
  * ascending vs descending order is actually honoured (bug 2.1),
  * sorting by column 0 is possible (not swallowed by a falsy check),
  * a string column index like "1" is coerced to int,
  * a call with no column specified raises instead of blocking on input().

Run from the project root:
    pytest tests/test_mutation/test_csv_sort.py
"""

import pandas as pd
import pytest

from mutadock.mutation.csv_sort import sort_csv
from mutadock.mutation.exceptions import MutationError

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _write_csv(tmp_path, rows, columns):
    """Write *rows* to a temp CSV and return its path (without the suffix trimmed)."""
    df = pd.DataFrame(rows, columns=columns)
    path = tmp_path / "input.csv"
    df.to_csv(path, index=False)
    return str(path)


def _read_csv(path):
    return pd.read_csv(path)


# ---------------------------------------------------------------------------
# Bug 2.1 — sort direction must be honoured
# ---------------------------------------------------------------------------


def test_descending_reverses_ascending(tmp_path):
    """Descending sort must produce the reverse row order of ascending."""
    # Column 0 ("sr") is dropped + renumbered by sort_csv, so sort by "value".
    rows = [
        [1, "alpha", 30],
        [2, "beta", 10],
        [3, "gamma", 20],
    ]
    in_file = _write_csv(tmp_path, rows, ["sr", "name", "value"])

    asc_out = str(tmp_path / "asc.csv")
    desc_out = str(tmp_path / "desc.csv")

    sort_csv(in_file=in_file, out_file=asc_out, col_num=2, order=True)
    sort_csv(in_file=in_file, out_file=desc_out, col_num=2, order=False)

    asc = _read_csv(asc_out)
    desc = _read_csv(desc_out)

    assert list(asc["value"]) == [10, 20, 30]
    assert list(desc["value"]) == [30, 20, 10]
    # Names should follow the same reversal.
    assert list(asc["name"]) == list(reversed(list(desc["name"])))


def test_ascending_is_default_order(tmp_path):
    """Default order (order not passed) sorts ascending."""
    rows = [
        [1, "a", 3],
        [2, "b", 1],
        [3, "c", 2],
    ]
    in_file = _write_csv(tmp_path, rows, ["sr", "name", "value"])
    out = str(tmp_path / "out.csv")

    sort_csv(in_file=in_file, out_file=out, col_num=2)

    assert list(_read_csv(out)["value"]) == [1, 2, 3]


# ---------------------------------------------------------------------------
# Bug 2.2 (1) — sorting by column 0 must work (0 is a valid index)
# ---------------------------------------------------------------------------


def test_sort_by_column_zero(tmp_path):
    """col_num=0 must sort by the first column, not trigger the fallback."""
    # sr is deliberately out of order so the sort is observable.
    rows = [
        [3, "alpha", 30],
        [1, "beta", 10],
        [2, "gamma", 20],
    ]
    in_file = _write_csv(tmp_path, rows, ["sr", "name", "value"])
    out = str(tmp_path / "out.csv")

    sort_csv(in_file=in_file, out_file=out, col_num=0, order=True)

    # Rows ordered by original sr (1, 2, 3) -> beta, gamma, alpha.
    assert list(_read_csv(out)["name"]) == ["beta", "gamma", "alpha"]


# ---------------------------------------------------------------------------
# Bug 2.2 (2) — string column index must be coerced to int
# ---------------------------------------------------------------------------


def test_string_column_index_is_coerced(tmp_path):
    """A string index like "1" must be accepted (coerced to int) and used."""
    rows = [
        [1, "gamma", 20],
        [2, "alpha", 30],
        [3, "beta", 10],
    ]
    in_file = _write_csv(tmp_path, rows, ["sr", "name", "value"])
    out = str(tmp_path / "out.csv")

    # "1" -> column index 1 == "name"
    sort_csv(in_file=in_file, out_file=out, col_num="1", order=True)

    assert list(_read_csv(out)["name"]) == ["alpha", "beta", "gamma"]


# ---------------------------------------------------------------------------
# Bug 2.2 (3) — no column specified must raise, not block on input()
# ---------------------------------------------------------------------------


def test_no_column_raises(tmp_path):
    """With neither col_num nor col_name, a MutationError is raised (no input())."""
    rows = [[1, "a", 1], [2, "b", 2]]
    in_file = _write_csv(tmp_path, rows, ["sr", "name", "value"])
    out = str(tmp_path / "out.csv")

    with pytest.raises(MutationError):
        sort_csv(in_file=in_file, out_file=out, col_num=None, col_name=None)


# ---------------------------------------------------------------------------
# col_name resolution still works
# ---------------------------------------------------------------------------


def test_sort_by_column_name_case_insensitive(tmp_path):
    """col_name resolves case-insensitively to the right column."""
    rows = [
        [1, "gamma", 20],
        [2, "alpha", 30],
        [3, "beta", 10],
    ]
    in_file = _write_csv(tmp_path, rows, ["sr", "name", "value"])
    out = str(tmp_path / "out.csv")

    sort_csv(in_file=in_file, out_file=out, col_name="VALUE", order=True)

    assert list(_read_csv(out)["value"]) == [10, 20, 30]
