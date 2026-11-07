"""Consistency checks for the bundled agent skills.

These keep the skills from silently going stale: the frontmatter has to follow
the Agent Skills convention, every skill has to be routable from ``INDEX.md``,
the Claude Code mirror has to match the canonical Codex tree, and no skill may
tell an agent to run a command that does not exist.

Pure stdlib -- no PyRosetta, no Vina, so these run everywhere.
"""

from __future__ import annotations

import re
import subprocess
import sys
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parent.parent.parent
SKILLS_DIR = REPO_ROOT / ".agents" / "skills"
MIRROR_DIR = REPO_ROOT / ".claude" / "skills"
INDEX = SKILLS_DIR / "INDEX.md"
PYPROJECT = REPO_ROOT / "pyproject.toml"
SYNC_SCRIPT = REPO_ROOT / "scripts" / "sync_agent_skills.py"

#: Per the Agent Skills spec; Claude Code truncates at 1536, Codex is stricter.
MAX_DESCRIPTION_CHARS = 1024
#: Both runtimes recommend keeping a SKILL.md short enough to load cheaply.
MAX_SKILL_LINES = 500


def skill_dirs() -> list[Path]:
    if not SKILLS_DIR.is_dir():
        return []
    return sorted(p for p in SKILLS_DIR.iterdir() if p.is_dir())


def parse_frontmatter(skill_md: Path) -> dict[str, str]:
    """Parse the simple ``key: value`` YAML frontmatter block of a SKILL.md."""
    text = skill_md.read_text(encoding="utf-8")
    match = re.match(r"^---\n(.*?)\n---\n", text, re.DOTALL)
    if not match:
        return {}
    fields: dict[str, str] = {}
    for line in match.group(1).splitlines():
        if ":" in line and not line.startswith((" ", "\t", "-")):
            key, value = line.split(":", 1)
            fields[key.strip()] = value.strip()
    return fields


SKILLS = skill_dirs()
SKILL_IDS = [p.name for p in SKILLS]


def test_skills_directory_is_populated():
    assert SKILLS, f"no skill directories found under {SKILLS_DIR}"


@pytest.mark.parametrize("skill", SKILLS, ids=SKILL_IDS)
def test_skill_has_skill_md(skill: Path):
    assert (skill / "SKILL.md").is_file(), f"{skill.name} has no SKILL.md"


@pytest.mark.parametrize("skill", SKILLS, ids=SKILL_IDS)
def test_frontmatter_name_matches_directory(skill: Path):
    fields = parse_frontmatter(skill / "SKILL.md")
    assert fields, f"{skill.name}/SKILL.md has no YAML frontmatter block"
    assert fields.get("name") == skill.name, (
        f"{skill.name}/SKILL.md declares name '{fields.get('name')}'; the Agent "
        f"Skills spec requires it to equal the directory name"
    )


@pytest.mark.parametrize("skill", SKILLS, ids=SKILL_IDS)
def test_frontmatter_description_is_present_and_portable(skill: Path):
    description = parse_frontmatter(skill / "SKILL.md").get("description", "")
    assert description, f"{skill.name}/SKILL.md has no description"
    assert len(description) <= MAX_DESCRIPTION_CHARS, (
        f"{skill.name} description is {len(description)} chars; keep it under "
        f"{MAX_DESCRIPTION_CHARS} so it is not truncated by either runtime"
    )


@pytest.mark.parametrize("skill", SKILLS, ids=SKILL_IDS)
def test_skill_md_is_not_too_long(skill: Path):
    lines = (skill / "SKILL.md").read_text(encoding="utf-8").splitlines()
    assert len(lines) <= MAX_SKILL_LINES, (
        f"{skill.name}/SKILL.md is {len(lines)} lines; move detail into "
        f"references/ and keep the entry point under {MAX_SKILL_LINES}"
    )


def test_every_skill_is_listed_in_index():
    index_text = INDEX.read_text(encoding="utf-8")
    missing = [s.name for s in SKILLS if s.name not in index_text]
    assert not missing, f"not routable from INDEX.md: {missing}"


def test_index_lists_no_unknown_skills():
    index_text = INDEX.read_text(encoding="utf-8")
    known = {s.name for s in SKILLS}
    referenced = set(re.findall(r"`(mutadock-[a-z0-9-]+)`", index_text))
    assert not (
        referenced - known
    ), f"INDEX.md routes to skills that do not exist: {sorted(referenced - known)}"


def test_claude_mirror_matches_canonical_tree():
    result = subprocess.run(
        [sys.executable, str(SYNC_SCRIPT), "--check"],
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0, (
        ".claude/skills has drifted from .agents/skills.\n"
        f"{result.stdout}\n{result.stderr}"
    )


def console_script_names() -> set[str]:
    """Return the ``md_*`` commands declared in ``[project.scripts]``."""
    text = PYPROJECT.read_text(encoding="utf-8")
    block = re.search(r"\[project\.scripts\]\n(.*?)(?=\n\[|\Z)", text, re.DOTALL)
    assert block, "could not find [project.scripts] in pyproject.toml"
    return set(re.findall(r"^(md_[a-z_]+)\s*=", block.group(1), re.MULTILINE))


@pytest.mark.parametrize("skill", SKILLS, ids=SKILL_IDS)
def test_skill_only_names_real_commands(skill: Path):
    """Catch a renamed, removed, or invented ``md_*`` command in a skill."""
    declared = console_script_names()
    for path in sorted(skill.rglob("*.md")):
        text = path.read_text(encoding="utf-8")
        used = set(re.findall(r"\bmd_[a-z_]+", text))
        # A shorthand stem such as `md_ddg_*` or `md_ddg_single|double|triple`
        # scans as the bare prefix `md_ddg_`; accept a token that prefixes at
        # least one real command, so only invented names fail.
        unknown = {
            name
            for name in used - declared
            if not any(cmd.startswith(name) for cmd in declared)
        }
        assert not unknown, (
            f"{path.relative_to(REPO_ROOT)} references command(s) that are not "
            f"in [project.scripts]: {sorted(unknown)}"
        )


def test_bundled_skills_are_discoverable_from_the_package():
    """``md_skills`` must find the skills from a source checkout too."""
    sys.path.insert(0, str(REPO_ROOT / "src"))
    try:
        from mutadock.skills import list_skills, skills_source
    finally:
        sys.path.pop(0)

    assert skills_source().is_dir()
    assert {p.name for p in list_skills()} == {s.name for s in SKILLS}


def _load_skills_module():
    sys.path.insert(0, str(REPO_ROOT / "src"))
    try:
        import mutadock.skills as mod

        return mod
    finally:
        sys.path.pop(0)


def test_install_never_deletes_its_own_source(tmp_path, monkeypatch, capsys):
    """Installing into the directory the skills are read from must be a no-op.

    ``md_skills install --project .`` inside a checkout resolves the Codex
    target to the canonical ``.agents/skills`` tree. Without a source ==
    destination guard the installer removes each source directory and then
    copies from the hole it just made, permanently destroying the skills.
    """
    mod = _load_skills_module()

    source = tmp_path / ".agents" / "skills"
    source.mkdir(parents=True)
    for name in ("mutadock-alpha", "mutadock-beta"):
        (source / name).mkdir()
        (source / name / "SKILL.md").write_text(
            f"---\nname: {name}\ndescription: test\n---\n", encoding="utf-8"
        )

    monkeypatch.setattr(mod, "skills_source", lambda: source)

    assert mod.install(["codex"], tmp_path, force=True) == 0

    survivors = sorted(p.name for p in source.iterdir())
    assert survivors == [
        "mutadock-alpha",
        "mutadock-beta",
    ], "install() destroyed its own source directory"
    for name in survivors:
        assert (source / name / "SKILL.md").is_file()
    assert "nothing to do" in capsys.readouterr().out


def test_install_copies_into_a_separate_directory(tmp_path, monkeypatch):
    """The ordinary case must still copy both agents' trees."""
    mod = _load_skills_module()

    source = tmp_path / "bundled"
    source.mkdir()
    (source / "mutadock-alpha").mkdir()
    (source / "mutadock-alpha" / "SKILL.md").write_text(
        "---\nname: mutadock-alpha\ndescription: test\n---\n", encoding="utf-8"
    )
    monkeypatch.setattr(mod, "skills_source", lambda: source)

    target = tmp_path / "project"
    target.mkdir()
    assert mod.install(["claude", "codex"], target, force=True) == 0
    assert (target / ".claude" / "skills" / "mutadock-alpha" / "SKILL.md").is_file()
    assert (target / ".agents" / "skills" / "mutadock-alpha" / "SKILL.md").is_file()
    assert (source / "mutadock-alpha" / "SKILL.md").is_file()
