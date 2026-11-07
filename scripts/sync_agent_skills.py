#!/usr/bin/env python3
"""Mirror ``.agents/skills`` (canonical, read by Codex) into ``.claude/skills``.

Claude Code discovers skills under ``.claude/skills/<name>/SKILL.md`` and does
not read ``.agents/``, so the repository carries both trees.  ``.agents/skills``
is the source of truth; this script regenerates the mirror.

Plain file copies are used deliberately: the repository is developed on Windows
and WSL, where git checks symlinks out as ordinary text files.

Usage::

    python scripts/sync_agent_skills.py            # regenerate the mirror
    python scripts/sync_agent_skills.py --check    # exit 1 if the trees differ
"""

from __future__ import annotations

import argparse
import filecmp
import shutil
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
SOURCE = REPO_ROOT / ".agents" / "skills"
MIRROR = REPO_ROOT / ".claude" / "skills"


def _relative_files(root: Path) -> set[Path]:
    if not root.is_dir():
        return set()
    return {
        p.relative_to(root)
        for p in root.rglob("*")
        if p.is_file() and "__pycache__" not in p.parts
    }


def check() -> list[str]:
    """Return a list of human-readable differences between the two trees."""
    src_files = _relative_files(SOURCE)
    dst_files = _relative_files(MIRROR)

    problems = [
        f"missing from .claude/skills: {p}" for p in sorted(src_files - dst_files)
    ]
    problems += [f"stale in .claude/skills: {p}" for p in sorted(dst_files - src_files)]
    problems += [
        f"differs: {p}"
        for p in sorted(src_files & dst_files)
        if not filecmp.cmp(SOURCE / p, MIRROR / p, shallow=False)
    ]
    return problems


def sync() -> int:
    if not SOURCE.is_dir():
        print(f"error: {SOURCE} does not exist", file=sys.stderr)
        return 1
    if MIRROR.exists():
        shutil.rmtree(MIRROR)
    shutil.copytree(SOURCE, MIRROR, ignore=shutil.ignore_patterns("__pycache__"))
    count = len(_relative_files(MIRROR))
    print(f"synced {count} file(s) to {MIRROR.relative_to(REPO_ROOT)}")
    return 0


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument(
        "--check",
        action="store_true",
        help="Verify the mirror matches without writing anything (exit 1 on drift).",
    )
    args = parser.parse_args(argv)

    if not args.check:
        return sync()

    problems = check()
    if problems:
        print(".claude/skills is out of date with .agents/skills:", file=sys.stderr)
        for p in problems:
            print(f"  - {p}", file=sys.stderr)
        print("\nRun: python scripts/sync_agent_skills.py", file=sys.stderr)
        return 1
    print(".claude/skills matches .agents/skills")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
