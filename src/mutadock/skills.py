"""Install MUTADOCK's bundled agent skills into a coding agent's skills folder.

MUTADOCK ships a set of task-focused ``SKILL.md`` guides so that Claude Code,
Codex and other agents that follow the Agent Skills convention already know the
CLI, the ΔΔG units and protocols, the cost model of a full scan, and the common
input failures.

The two agents look in different places, so this command writes to either or
both:

* Claude Code -- ``~/.claude/skills/<name>/`` (or ``<project>/.claude/skills/``)
* Codex       -- ``~/.agents/skills/<name>/`` (or ``<project>/.agents/skills/``)

Only directories named ``mutadock-*`` are ever written or replaced, so a shared
skills folder containing other people's skills is safe.
"""

from __future__ import annotations

import argparse
import shutil
import sys
from pathlib import Path

#: Every bundled skill directory starts with this prefix.
SKILL_PREFIX = "mutadock-"

#: Where each supported agent discovers skills, relative to a scope root.
AGENT_DIRS = {
    "claude": Path(".claude") / "skills",
    "codex": Path(".agents") / "skills",
}


def skills_source() -> Path:
    """Return the directory holding the bundled ``mutadock-*`` skills.

    In an installed wheel the skills live in ``mutadock/agent_skills`` (placed
    there by hatchling's ``force-include``).  In a source checkout -- a clone, an
    editable install, or a pytest run using ``pythonpath = ["src"]`` -- that
    directory does not exist, so fall back to the canonical ``.agents/skills``
    tree at the repository root.

    Raises:
        FileNotFoundError: If neither location exists.
    """
    packaged = Path(__file__).resolve().parent / "agent_skills"
    if packaged.is_dir():
        return packaged

    # src/mutadock/skills.py -> src/mutadock -> src -> <repo root>
    repo_root = Path(__file__).resolve().parent.parent.parent
    source_tree = repo_root / ".agents" / "skills"
    if source_tree.is_dir():
        return source_tree

    raise FileNotFoundError(
        "Could not locate the bundled MUTADOCK skills. Looked in "
        f"'{packaged}' and '{source_tree}'."
    )


def list_skills() -> list[Path]:
    """Return the bundled skill directories, sorted by name."""
    source = skills_source()
    return sorted(
        p
        for p in source.iterdir()
        if p.is_dir() and p.name.startswith(SKILL_PREFIX) and (p / "SKILL.md").is_file()
    )


def _describe(skill_dir: Path) -> str:
    """Return the ``description:`` line of a skill, truncated for display."""
    try:
        lines = (skill_dir / "SKILL.md").read_text(encoding="utf-8").splitlines()
    except OSError:
        return ""
    for line in lines[:20]:
        if line.startswith("description:"):
            text = line.split(":", 1)[1].strip()
            return text if len(text) <= 100 else text[:97] + "..."
    return ""


def _same_path(a: Path, b: Path) -> bool:
    """True when two paths refer to the same location on disk.

    ``Path.samefile`` is authoritative but needs both paths to exist, so fall
    back to comparing fully resolved paths.
    """
    try:
        if a.exists() and b.exists():
            return a.samefile(b)
    except OSError:
        pass
    return a.resolve() == b.resolve()


def install(
    targets: list[str],
    scope_root: Path,
    force: bool = False,
) -> int:
    """Copy every bundled skill into each target agent's skills directory.

    Args:
        targets: Agent keys from :data:`AGENT_DIRS` (``"claude"``, ``"codex"``).
        scope_root: Directory the agent paths are resolved against -- the user's
            home for ``--user``, or a project directory for ``--project``.
        force: Replace an existing ``mutadock-*`` directory instead of skipping.

    Returns:
        Process exit code: 0 on success, 1 if anything was skipped.
    """
    skills = list_skills()
    if not skills:
        print("error: no bundled skills found", file=sys.stderr)
        return 1

    skipped = False
    source_root = skills_source()
    for target in targets:
        dest_root = scope_root / AGENT_DIRS[target]
        if _same_path(dest_root, source_root):
            print(
                f"\n{target}: {dest_root}\n"
                "  nothing to do - this is where the bundled skills already live."
            )
            continue
        dest_root.mkdir(parents=True, exist_ok=True)
        print(f"\n{target}: {dest_root}")

        for skill in skills:
            dest = dest_root / skill.name

            # Installing into the directory the skills are read from would
            # delete the source and then copy from the hole it just made. This
            # happens for real: `--project .` inside a checkout resolves
            # `.agents/skills` to the canonical tree itself.
            if _same_path(dest, skill):
                print(f"  skipped  {skill.name} (source and destination are the same)")
                continue

            if dest.exists() and not force:
                print(f"  skipped  {skill.name} (already exists; use --force)")
                skipped = True
                continue
            if dest.exists():
                shutil.rmtree(dest)
            shutil.copytree(skill, dest, ignore=shutil.ignore_patterns("__pycache__"))
            print(f"  installed {skill.name}")

    if skipped:
        print(
            "\nSome skills already existed and were left alone. Re-run with --force to replace them."
        )
        return 1
    print(
        "\nDone. Restart your agent (or run /reload-skills in Claude Code) to pick them up."
    )
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="md_skills",
        description="Install MUTADOCK's agent skills for Claude Code and Codex.",
    )
    sub = parser.add_subparsers(dest="command")

    sub.add_parser("list", help="List the bundled skills.")
    sub.add_parser("path", help="Print the directory holding the bundled skills.")

    inst = sub.add_parser(
        "install", help="Copy the skills into an agent's skills folder."
    )
    inst.add_argument(
        "--claude",
        action="store_true",
        help="Install for Claude Code (.claude/skills).",
    )
    inst.add_argument(
        "--codex",
        action="store_true",
        help="Install for Codex (.agents/skills).",
    )
    scope = inst.add_mutually_exclusive_group()
    scope.add_argument(
        "--user",
        action="store_true",
        help="Install into the home directory, for every project (the default).",
    )
    scope.add_argument(
        "--project",
        metavar="DIR",
        help="Install into DIR, for that project only.",
    )
    inst.add_argument(
        "--force",
        action="store_true",
        help="Replace existing mutadock-* skill directories.",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    """CLI entry point for ``md_skills``."""
    parser = build_parser()
    args = parser.parse_args(argv)

    try:
        if args.command == "path":
            print(skills_source())
            return 0

        if args.command == "list" or args.command is None:
            for skill in list_skills():
                print(f"{skill.name:24} {_describe(skill)}")
            if args.command is None:
                print("\nRun 'md_skills install --claude --codex' to install them.")
            return 0

        targets = [t for t in ("claude", "codex") if getattr(args, t)]
        if not targets:
            targets = ["claude", "codex"]

        scope_root = (
            Path(args.project).expanduser().resolve() if args.project else Path.home()
        )
        if args.project and not scope_root.is_dir():
            print(
                f"error: project directory '{scope_root}' does not exist",
                file=sys.stderr,
            )
            return 1

        return install(targets, scope_root, force=args.force)
    except FileNotFoundError as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
